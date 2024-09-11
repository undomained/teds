// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class to hold partially or fully calibrated level 1 measurements. 
// The data level can range from L1A to L1B depending on the calibration

#include "l1_measurement.h"
#include "l1.h"
#include "time.h"

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

namespace tango {

L1Measurement::L1Measurement(const std::string filename, const int image_start, const int image_end, const std::string config)
{

    // When we create a L1Measurement object need to create it from file
    // So l1_filename should be a member
    l1_filename = filename;
    const netCDF::NcFile nc { l1_filename, netCDF::NcFile::read };
    setLevel( nc);

    // TODO do we want to obtain image_start and image_end from settings?
    alt_beg = static_cast<size_t>(image_start);
    int n_alt = nc.getDim("along_track").getSize();
    if (image_end == fill::i){
        alt_end = n_alt;
    } else {
        alt_end = static_cast<size_t>(image_end) + 1;
    }
    n_images = alt_end - alt_beg ;

    if (n_images > (n_alt - alt_beg)){
        spdlog::warn("Given image_end exceeds number of images in input");
        n_images = n_alt - alt_beg; // set last image to analyze to last image possible
    }

    l1_measurement.resize(n_images);
    image_time.resize(n_images);
    binning_table.resize(n_images);
    nr_coadditions.resize(n_images);
    exposure_time.resize(n_images);

    read(l1_filename, config);
//    copyGeometry(config);
}

L1& L1Measurement::operator[](const int img){
    // acces specific image (L1 object)
    return l1_measurement[img];
}

int L1Measurement::size(){
    // Return the size i.e number of images (L1 objects)
    return l1_measurement.size();
}

L1& L1Measurement::front(){
    // Return the first image (L1 object)
    return l1_measurement.front();
}

L1& L1Measurement::back(){
    // Return the last image (L1 object)
    return l1_measurement.back();
}

void L1Measurement::read(const std::string& filename,
    const std::string& config)
{
    l1_filename = filename;
    const netCDF::NcFile nc { l1_filename, netCDF::NcFile::read };
    // possibly given file name is not equal to filename at construction.
    // need to set level
    setLevel(nc);

    readMetaData(nc, config);
    if (l1_level == "SGM"){
        readSceneData(nc);
    } else {
        readScienceData(nc);
    }
}

// Basically Edwards writeL1product fct from io.cpp
void L1Measurement::write(const std::string& filename,
    const std::string& level,
    const std::string& config, 
    const int argc, const char* const argv[]) 
{
    spdlog::info("write L1: {}", filename);
    
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };

    writeGlobalAttributes(nc, level, config,argc, argv);
    writeMetaData(nc);
    
    if (level == "SGM"){
        // Need to check if proctable has only ISRF in it.
        writeObservationData(nc);
    } else {
        writeScienceData(nc);
//        writeGeolocationData(nc);
    }

}

void L1Measurement::setLevel(const netCDF::NcFile& nc){
    // Determine data level
    std::string level {};
    try {
        nc.getAtt("product_type").getValues(level);
    } catch (const netCDF::exceptions::NcBadId&) {
        spdlog::error("Can't find product_type variable in input data, aborting");
        std::exit(EXIT_FAILURE);
    }
    spdlog::info("Input data calibration level: {}", level);
    l1_level = level;
}

void L1Measurement::writeGlobalAttributes(netCDF::NcFile& nc, const std::string& level, const std::string& config, const int argc, const char* const argv[]) {
    // Edwards splitup functions copied from io.cpp
    // TODO kunnen we niet settings gebruiken? Dat is toch al geladen config??????? Dan niet YAML::load nodig (?????)
    spdlog::info("Writing Global Attributes");
    const std::string instrument { YAML::Load(config)["instrument"].as<std::string>() };
    nc.putAtt("Conventions", "CF-1.11");
    nc.putAtt("title", "Tango " + instrument + " " + level + " data");
    nc.putAtt("product_type", level);
    nc.putAtt("project", "TANGO");
    nc.putAtt("instrument", instrument);
    nc.putAtt("creator_name", "SRON/Earth Science");
    nc.putAtt("creator_url", "https://www.sron.nl/missions-earth");
    nc.putAtt("date_created", getDateAndTime());

    if (TANGO_GIT_COMMIT_ABBREV != "GITDIR-N") {
        nc.putAtt("git_commit", TANGO_GIT_COMMIT_ABBREV);
    }

    std::string command_line { argc > 0 ? argv[0] : "" };
    for (int i { 1 }; i < argc; ++i) {
        command_line += ' ' + std::string(argv[i]);
    }
    nc.putAtt("history", command_line);
    const std::string processing_version { YAML::Load(config)["processing_version"].as<std::string>() };
    nc.putAtt("processing_version", processing_version);
}

void L1Measurement::writeMetaData(netCDF::NcFile& nc) {
    // Edwards splitup functions copied from io.cpp with l1_products -> l1_measurement
    spdlog::info("Writing Image Attributes");
    const auto nc_images { nc.addDim("along_track", n_images) };

    auto nc_grp { nc.addGroup("image_attributes") };

    NcVar nc_img_time = addVariable(nc_grp, "image_time", "image time", "seconds since 2022-03-21", fill::d, 0.0, 92304.0, {nc_images});
    std::vector<double> buf {};
    buf.assign(n_images, 0.0);
    nc_img_time.putVar(buf.data());

    auto nc_var = nc_grp.addVar("binning_table", netCDF::ncByte, { nc_images });
    nc_var.putAtt("long_name", "binning table");
    buf.assign(n_images, l1_measurement.front().binning_table_id);
    nc_var.putVar(buf.data());

    nc_var = nc_grp.addVar("nr_coadditions", netCDF::ncUshort, { nc_images });
    nc_var.putAtt("long_name", "number of coadditions");
    buf.assign(n_images, l1_measurement.front().nr_coadditions);
    nc_var.putVar(buf.data());

    nc_var = nc_grp.addVar("exposure_time", netCDF::ncDouble, { nc_images });
    nc_var.putAtt("long_name", "exposure time");
    buf.assign(n_images, l1_measurement.front().exposure_time);
    nc_var.putVar(buf.data());

}

// May be obsolete?
void L1Measurement::writeObservationData(netCDF::NcFile& nc) {
    // Edwards splitup functions copied from io.cpp with l1_products -> l1_measurement
    spdlog::info("Writing Observation Data");
    const auto n_across_track { l1_measurement.front().observation_sig.size() };    
    const auto nc_across_track { nc.addDim("across_track", n_across_track) };
    const auto nc_images = nc.getDim("along_track");
    auto& l1 { l1_measurement[0] };

    const auto& wavelengths = *l1_measurement.front().wavelength; 
    const auto n_wavelength { wavelengths.front().size() };
    const auto nc_wavelength { nc.addDim("wavelength", n_wavelength) };
    
    auto nc_grp { nc.addGroup("observation_data") };

    // add wavelength
    NcVar nc_wl = addVariable(nc_grp, "wavelength", "radiance wavelengths", "nm", fill::f, 0.0f, 999.0f, {nc_across_track, nc_wavelength});
    std::vector<float> buf3(n_across_track * n_wavelength);
    for (size_t i {}; i < n_across_track; ++i) {
        for (size_t j {}; j < n_wavelength; ++j) {
            buf3[i * n_wavelength + j] = static_cast<float>(wavelengths[i][j]);
        }
    }
    nc_wl.putVar(buf3.data());

    std::string& units = l1_measurement.front().units;
    NcVar nc_rad = addVariable(nc_grp, "radiance", "radiance", units , fill::f, 0.0f, 1e20f, {nc_images, nc_across_track, nc_wavelength});
    NcVar nc_std = addVariable(nc_grp, "radiance_stdev", "standard deviation of radiance in bin",units, fill::f, 0.0f, 1e20f, {nc_images, nc_across_track, nc_wavelength});

    std::vector<double> buf(n_images * n_across_track * n_wavelength);
    std::vector<double> buf2(n_images * n_across_track * n_wavelength);

    for (size_t i {}; i < n_images; i++) {
        for (size_t j {}; j < n_across_track; j++) {
            for (size_t k {}; k < n_wavelength; k++) {
                int index = i * n_across_track * n_wavelength + j * n_wavelength + k;
                buf[index] = l1_measurement[i].observation_sig[j][k];
                buf2[index] = l1_measurement[i].observation_std[j][k];
            }
        }
    }
    spdlog::info("Writing radiance and radiance_std with units: {}", units);
    nc_rad.putVar(buf.data());
    nc_std.putVar(buf2.data());

}

void L1Measurement::writeScienceData(netCDF::NcFile& nc){
    // Edwards splitup functions copied from io.cpp with l1_products -> l1_measurement
    spdlog::info("Writing Science Data");
    // Note: n_images is member of L1Measurements
    //    const auto n_images { l1_products.size() };
    const auto n_bins { l1_measurement.front().image.size() };

    const auto nc_images = nc.getDim("along_track");
    const auto nc_detector_bin { nc.addDim("detector_bin", n_bins) };

    const auto n_rows {(*l1_measurement.front().wavelength).size()};
    const auto nc_rows { nc.addDim("rows", n_rows) };

    const auto& wavelengths = *l1_measurement.front().wavelength; 
    const auto n_cols { wavelengths.front().size() };
    const auto nc_cols { nc.addDim("cols", n_cols) };
    
    auto nc_grp { nc.addGroup("science_data") };
    std::string& units = l1_measurement.front().units;
    // add image and stdev
    NcVar nc_img = addVariable(nc_grp, "detector_image", "detector images", units , fill::d, -1e100, 1e100, {nc_images, nc_rows, nc_cols});
    NcVar nc_std = addVariable(nc_grp, "detector_stdev", "standard deviation of detector bin",units, fill::d, 0.0,1e100, {nc_images, nc_rows, nc_cols});
    std::vector<double> buf(n_images * n_rows * n_cols);
    std::vector<double> buf2(n_images * n_rows * n_cols);
    for (size_t i {}; i < n_images; i++) {
        for (size_t j {}; j < n_rows; j++) {
            for (size_t k {}; k < n_cols; k++) {
                buf[i * n_rows * n_cols + j * n_cols + k] = l1_measurement[i].image[j * n_cols + k];
                buf2[i * n_rows * n_cols + j * n_cols + k] = l1_measurement[i].stdev[j * n_cols + k];
            }
        }
    }
    nc_img.putVar(buf.data());
    nc_std.putVar(buf2.data());

    // Add wavelength
    NcVar nc_wl = addVariable(nc_grp, "wavelength", "radiance wavelengths", "nm", fill::f, 0.0f, 999.0f, {nc_rows, nc_cols});
    std::vector<float> buf3(n_rows * n_cols);
    for (size_t i {}; i < n_rows; i++) {
        for (size_t j {}; j < n_cols; j++) {
            buf3[i * n_cols + j] = static_cast<float>(wavelengths[i][j]);
        }
    }
    nc_wl.putVar(buf3.data());
}

void L1Measurement::writeGeolocationData(netCDF::NcFile& nc) {
    // Edwards splitup functions copied from io.cpp with l1_products -> l1_measurement
    spdlog::info("Writing Geolocation Data");

    const auto n_across_track { l1_measurement.front().spectra.size() };
    const auto nc_across_track { nc.addDim("across_track", n_across_track) };
    const auto nc_images = nc.getDim("along_track");
    const std::vector<netCDF::NcDim> geometry_shape { nc_images, nc_across_track };

    std::vector<float> lat(n_images * n_across_track);
    std::vector<float> lon(n_images * n_across_track);
    std::vector<float> height(n_images * n_across_track);
    std::vector<float> vza(n_images * n_across_track);
    std::vector<float> vaa(n_images * n_across_track);
    std::vector<float> sza(n_images * n_across_track);
    std::vector<float> saa(n_images * n_across_track);
    for (int i_alt {}; i_alt < n_images; ++i_alt) {
        const auto copy { [i_alt](const std::vector<float>& in,
                                    std::vector<float>& out) {
            std::copy(
                in.cbegin(), in.cend(), out.begin() + i_alt * in.size());
        } };
        copy(l1_measurement[i_alt].geo.lat, lat);
        copy(l1_measurement[i_alt].geo.lon, lon);
        copy(l1_measurement[i_alt].geo.height, height);
        copy(l1_measurement[i_alt].geo.vza, vza);
        copy(l1_measurement[i_alt].geo.vaa, vaa);
        copy(l1_measurement[i_alt].geo.sza, sza);
        copy(l1_measurement[i_alt].geo.saa, saa);
    }

    auto nc_grp = nc.addGroup("geolocation_data");

    NcVar nc_lat = addVariable(nc_grp, "latitude", "latitude at bin locations", "degrees_north", fill::f, -90.0f, 90.0f, geometry_shape);
    NcVar nc_lon = addVariable(nc_grp, "longitude", "longitude at bin locations", "degrees_east", fill::f, -180.0f, 180.0f, geometry_shape);
    NcVar nc_hei = addVariable(nc_grp, "height", "height at bin locations", "m", fill::f, -1000.0f, 1000.0f, geometry_shape);
    NcVar nc_vza = addVariable(nc_grp, "viewingzenithangle", "sensor zenith angle at bin locations", "degrees", fill::f, -90.0f, 90.0f, geometry_shape);
    NcVar nc_vaa = addVariable(nc_grp, "viewingazimuthangle", "sensor azimuth angle at bin locations", "degrees", fill::f, -180.0f, 180.0f, geometry_shape);
    NcVar nc_sza = addVariable(nc_grp, "solarzenithangle", "solar zenith angle at bin locations", "degrees", fill::f, -90.0f, 90.0f, geometry_shape);
    NcVar nc_saa = addVariable(nc_grp, "solarazimuthangle", "solar azimuth angle at bin locations", "degrees", fill::f, -180.0f, 180.0f, geometry_shape);
    
    nc_lon.putVar(lon.data());
    nc_lat.putVar(lat.data());
    nc_hei.putVar(height.data());
    nc_vza.putVar(vza.data());
    nc_vaa.putVar(vaa.data());
    nc_sza.putVar(sza.data());
    nc_saa.putVar(saa.data());

}

void L1Measurement::readMetaData(const netCDF::NcFile& nc, const std::string& config) {
    // Edwards Meta data part from readL1 fct from io_nitro.cpp?
    spdlog::info("Reading Metadata");
    if (l1_level == "SGM") {
        uint8_t detector_binning_table_id { YAML::Load(config)["detector"]["binning_table_id"].as<uint8_t>() };
        uint16_t nr_coadditions { YAML::Load(config)["detector"]["nr_coadditions"].as<uint16_t>() };
        float exposure_time { YAML::Load(config)["detector"]["exposure_time"].as<float>() };
        for (int i_alt {}; i_alt < static_cast<int>(n_images); ++i_alt) {
            l1_measurement[i_alt].binning_table_id =  detector_binning_table_id;
            l1_measurement[i_alt].nr_coadditions = nr_coadditions;
            l1_measurement[i_alt].exposure_time =  exposure_time;
        }
    } else {
        if (const auto grp { nc.getGroup("image_attributes") }; !grp.isNull()) {
            grp.getVar("image_time").getVar({ alt_beg }, { n_images }, image_time.data());
            grp.getVar("binning_table").getVar({ alt_beg }, { n_images }, binning_table.data());
            grp.getVar("nr_coadditions").getVar({ alt_beg }, { n_images }, nr_coadditions.data());
            grp.getVar("exposure_time").getVar({ alt_beg }, { n_images }, exposure_time.data());
            for (int i_alt {}; i_alt < static_cast<int>(n_images); ++i_alt) {
                l1_measurement[i_alt].image_time = image_time[i_alt];
                l1_measurement[i_alt].binning_table_id = binning_table[i_alt];
                l1_measurement[i_alt].nr_coadditions = nr_coadditions[i_alt];
                l1_measurement[i_alt].exposure_time = exposure_time[i_alt];
            }
        }
    }
    
}

void L1Measurement::readScienceData(const netCDF::NcFile& nc) {
    // Edwards Science data part from readL1 fct from io_nitrp.cpp
    spdlog::info("Reading Science Data");
    const auto grp { nc.getGroup("science_data") };
    const auto n_bins { nc.getDim("detector_bin").getSize() }; // == rows * cols
    const auto n_rows { nc.getDim("rows").getSize() };
    const auto n_cols { nc.getDim("cols").getSize() };

    std::vector<double> detector_images(n_images * n_bins);
    std::vector<double> detector_stdev(n_images * n_bins);
    std::vector<double> flattened_wl(n_bins);
    std::string units {};

    auto nc_img =  grp.getVar("detector_image");
    nc_img.getVar({ alt_beg, 0, 0 }, { n_images, n_rows, n_cols }, detector_images.data());
    nc_img.getAtt("units").getValues(units);

    grp.getVar("detector_stdev").getVar({ alt_beg, 0, 0 }, { n_images, n_rows, n_cols }, detector_stdev.data());
    grp.getVar("wavelength").getVar(flattened_wl.data());
   
    int pix {};
    for (size_t i_alt {}; i_alt < n_images; ++i_alt) {
        auto& l1 { l1_measurement[i_alt] };
        l1.image.resize(n_bins);
        l1.stdev.resize(n_bins);
        l1.units = units;
        for (size_t i {}; i < n_rows; i++) {
            for (size_t j {}; j < n_cols; j++) {
                pix = i * n_cols + j; // pixel index
                l1.image[pix] = detector_images[i_alt * n_rows * n_cols + pix];
                l1.stdev[pix] = detector_stdev[i_alt * n_rows * n_cols + pix];
            }
        }
    }
    
    // Import wavelength
    auto& l1 {l1_measurement.front()};
    std::vector<std::vector<double>> wavelength2D(n_rows, std::vector<double>(n_cols));
    for (int i {}; i < n_rows; ++i) {
        for (int j {}; j < n_cols; ++j) {
            wavelength2D[i][j] = flattened_wl[i * n_cols + j];
        }
    } 
    l1.wavelength = std::make_shared<std::vector<std::vector<double>>>(wavelength2D);
}

void L1Measurement::readSceneData(const netCDF::NcFile& nc){
    // Edwards read Scene data part from readL1 fct from io_nitro.cpp
    spdlog::info("Reading Scene Data");
    const auto n_act { nc.getDim("across_track").getSize() };
    const auto n_wavelength { nc.getDim("wavelength").getSize() };
    const auto n_alt { nc.getDim("along_track").getSize() };

    std::vector<double> spectra(n_images * n_act * n_wavelength);
    std::vector<double> wavelength(n_wavelength);
    std::string units {};

    auto nc_rad =  nc.getVar("radiance");
    nc_rad.getVar({ alt_beg, 0, 0 }, { n_images, n_act, n_wavelength }, spectra.data());
    nc_rad.getAtt("units").getValues(units);
//    nc.getVar("radiance").getVar({ alt_beg, 0, 0 }, { n_images, n_act, n_wavelength }, spectra.data());
//    nc_img.getAtt("units").getValues(units);
    nc.getVar("wavelength").getVar(wavelength.data());

    // import wavelength
    auto wavelength_lbl {std::make_shared<std::vector<std::vector<double>>>(n_act, std::vector<double>(n_wavelength))};
    for (size_t i_act {}; i_act < n_act; ++i_act) {
        for (size_t i {}; i < n_wavelength; ++i) {
            (*wavelength_lbl)[i_act][i] = wavelength[i];
        }
    }

    // set wavelengths
    for (L1& l1 : l1_measurement) {
        l1.wavelength = wavelength_lbl;
        l1.observation_wl = wavelength_lbl;
    }

    //import spectra
    for (size_t i_alt {}; i_alt < n_images; ++i_alt) {
        auto& l1 { l1_measurement[i_alt] };
        l1.observation_sig.resize(n_act);
        l1.observation_std.resize(n_act);
        l1.units = units;
        for (size_t i_act {}; i_act < n_act; ++i_act) {
            l1.observation_sig[i_act].resize(n_wavelength);
            l1.observation_std[i_act].resize(n_wavelength);
            for (size_t i {}; i < n_wavelength; ++i) {
                const size_t idx { (i_alt * n_act + i_act) * n_wavelength + i };
                l1.observation_sig[i_act][i] = spectra[idx];
            }
        }
    }
}

void L1Measurement::copyGeometry(const std::string& config)
{
    // The following line gives a problem. Not sure why.
    const std::string geo_filename { YAML::Load(config)["io"]["geometry"].as<std::string>() };

    const netCDF::NcFile nc_geo { geo_filename, netCDF::NcFile::read };
    const netCDF::NcGroup nc { nc_geo.getGroup("/") };
    const auto n_alt { nc_geo.getDim("along_track").getSize() };
    const auto n_act { nc_geo.getDim("across_track").getSize() };
    std::vector<double> lat(n_alt * n_act);
    std::vector<double> lon(n_alt * n_act);
    std::vector<double> vza(n_alt * n_act);
    std::vector<double> vaa(n_alt * n_act);
    std::vector<double> sza(n_alt * n_act);
    std::vector<double> saa(n_alt * n_act);
    nc_geo.getVar("latitude").getVar(lat.data());
    nc_geo.getVar("longitude").getVar(lon.data());
    nc_geo.getVar("viewingzenithangle").getVar(vza.data());
    nc_geo.getVar("viewingazimuthangle").getVar(vaa.data());
    nc_geo.getVar("solarzenithangle").getVar(sza.data());
    nc_geo.getVar("solarazimuthangle").getVar(saa.data());
    int i_alt_start = static_cast<int>(alt_beg);
    for (int i_alt {}; i_alt < static_cast<int>(l1_measurement.size()); ++i_alt) {
        const auto copy { [i_alt_start, i_alt, n_act](
                            const std::vector<double>& in,
                            std::vector<float>& out) {
            out = { in.begin() + (i_alt_start + i_alt) * n_act,
                    in.begin() + (i_alt_start + i_alt + 1) * n_act };
        } };
        copy(lat, l1_measurement[i_alt].geo.lat);
        copy(lon, l1_measurement[i_alt].geo.lon);
        copy(vza, l1_measurement[i_alt].geo.vza);
        copy(vaa, l1_measurement[i_alt].geo.vaa);
        copy(sza, l1_measurement[i_alt].geo.sza);
        copy(saa, l1_measurement[i_alt].geo.saa);
        l1_measurement[i_alt].geo.height = std::vector<float>(n_act, 0.0);
    }
}

} // namespace tango
