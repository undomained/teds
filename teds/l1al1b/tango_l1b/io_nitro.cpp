// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.


#include "io_nitro.h"
#include "io.h"

#include "l1.h"
#include "time.h"

#include <filesystem>
#include <fstream>
#include <omp.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>


namespace tango {

void writeL1product(
    const std::string& filename, 
    const std::string& level,
    const std::string& config, 
    const std::vector<L1>& l1_products, 
    const int argc, const char* const argv[]) 
    {
    
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    spdlog::info("Output data calibration level: {}", level);

    writeGlobalAttributes(nc, level, config,argc, argv);
    writeMetaData(nc, l1_products);

    if (level == "L1B"){ 
        writeObservationData(nc, l1_products);
    }
    writeScienceData(nc, l1_products);
    //writeGeolocationData(nc, l1_products);

    spdlog::info("Succesfully written to file: {}", filename);
}

void writeGlobalAttributes(netCDF::NcFile& nc, const std::string& level, const std::string& config, const int argc, const char* const argv[]) {
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

void writeMetaData(netCDF::NcFile& nc, const std::vector<L1>& l1_products) {
    spdlog::info("Writing Image Attributes");
    const auto n_images { l1_products.size() };
    const auto nc_images { nc.addDim("along_track", n_images) };

    auto nc_grp { nc.addGroup("image_attributes") };

    NcVar nc_img_time = addVariable(nc_grp, "image_time", "image time", "seconds since 2022-03-21", fill::d, 0.0, 92304.0, {nc_images});
    std::vector<double> buf {};
    buf.assign(n_images, 0.0);
    nc_img_time.putVar(buf.data());

    auto nc_var = nc_grp.addVar("binning_table", netCDF::ncByte, { nc_images });
    nc_var.putAtt("long_name", "binning table");
    buf.assign(n_images, l1_products.front().binning_table_id);
    nc_var.putVar(buf.data());

    nc_var = nc_grp.addVar("nr_coadditions", netCDF::ncUshort, { nc_images });
    nc_var.putAtt("long_name", "number of coadditions");
    buf.assign(n_images, l1_products.front().nr_coadditions);
    nc_var.putVar(buf.data());

    nc_var = nc_grp.addVar("exposure_time", netCDF::ncDouble, { nc_images });
    nc_var.putAtt("long_name", "exposure time");
    buf.assign(n_images, l1_products.front().exposure_time);
    nc_var.putVar(buf.data());
}

void writeScienceData(netCDF::NcFile& nc, const std::vector<L1>& l1_products) {
    spdlog::info("Writing Science Data");
    const auto n_images { l1_products.size() };
    const auto n_bins { l1_products.front().image.size()};

    const auto nc_images = nc.getDim("along_track");
    const auto nc_detector_bin { nc.addDim("detector_bin", n_bins) };

    const auto n_rows {(*l1_products.front().wavelength).size()};
    const auto nc_rows { nc.addDim("rows", n_rows) };

    const auto& wavelengths = *l1_products.front().wavelength; 
    const auto n_cols { wavelengths.front().size() };
    const auto nc_cols { nc.addDim("cols", n_cols) };

    const std::vector<netCDF::NcDim> nc_flatimg_shape { nc_images, nc_detector_bin };
    
    auto nc_grp { nc.addGroup("science_data") };
    
    // Add image
    NcVar nc_img = addVariable(nc_grp, "detector_image", "detector images", "counts" , fill::d, -1e100, 1e100, nc_flatimg_shape);
    std::vector<double> buf(n_images * n_bins);
    for (size_t i {}; i < n_images; ++i) {
        for (size_t j {}; j < n_bins; ++j) {
            buf[i * n_bins + j] = l1_products[i].image[j];
        }
    }
    nc_img.putVar(buf.data());
    
    // Add stdev
    NcVar nc_std = addVariable(nc_grp, "detector_stdev", "standard deviation of detector bin","counts", fill::d, 0.0,1e100, nc_flatimg_shape);
    for (size_t i {}; i < n_images; ++i) {
        for (size_t j {}; j < n_bins; ++j) {
            buf[i * n_bins + j] = l1_products[i].stdev[j];
        }
    }
    nc_std.putVar(buf.data());
    
    // Add wavelength
    NcVar nc_wl = addVariable(nc_grp, "wavelength", "radiance wavelengths", "nm", fill::f, 0.0f, 999.0f, {nc_rows, nc_cols});
    std::vector<float> buf2(n_rows * n_cols);
    for (size_t i {}; i < n_rows; ++i) {
        for (size_t j {}; j < n_cols; ++j) {
            buf2[i * n_cols + j] = static_cast<float>(wavelengths[i][j]);
        }
    }
    nc_wl.putVar(buf2.data());
}

void writeObservationData(netCDF::NcFile& nc, const std::vector<L1>& l1_products) {
    spdlog::info("Writing Observation Data");
    const auto n_across_track { l1_products.front().observation_sig.size() };    
    const auto nc_across_track { nc.addDim("across_track", n_across_track) };
    const auto n_images = l1_products.size();
    const auto nc_images = nc.getDim("along_track");
    auto& l1 { l1_products[0] };

    const auto& wavelengths = *l1_products.front().wavelength; 
    const auto n_wavelength { wavelengths.front().size() };
    const auto nc_wavelength { nc.addDim("wavelength", n_wavelength) };
    
    auto nc_grp { nc.addGroup("observation_data") };

    // add wavelength
    NcVar nc_wl = addVariable(nc_grp, "wavelength", "radiance wavelengths", "nm", fill::f, 0.0f, 999.0f, {nc_across_track, nc_wavelength});
    std::vector<float> buf(n_across_track * n_wavelength);
    for (size_t i {}; i < n_across_track; ++i) {
        for (size_t j {}; j < n_wavelength; ++j) {
            buf[i * n_wavelength + j] = static_cast<float>(wavelengths[i][j]);
        }
    }
    nc_wl.putVar(buf.data());

    // add radiance, flatten first
    const auto flatsignal { [&](auto& buf) {
        for (size_t i {}; i < n_images; ++i) {
            for (size_t j {}; j < n_across_track; ++j) {
                for (size_t k {}; k < n_wavelength; ++k) {
                    buf[(i * n_across_track + j) * n_wavelength + k] = static_cast<float>(l1_products[i].observation_sig[j][k]);
                }
            }
        }
    } };
    NcVar nc_rad = addVariable(nc_grp, "radiance", "radiance", "ph nm-1 s-1 sr-1 m-2", fill::f, 0.0f, 1e20f, {nc_images, nc_across_track, nc_wavelength});
    flatsignal(buf);
    nc_rad.putVar(buf.data());

    const auto flatstd { [&](auto& buf) {
        for (size_t i {}; i < n_images; ++i) {
            for (size_t j {}; j < n_across_track; ++j) {
                for (size_t k {}; k < n_wavelength; ++k) {
                    buf[(i * n_across_track + j) * n_wavelength + k] = static_cast<float>(l1_products[i].observation_std[j][k]);
                }
            }
        }
    }}; 
    
    NcVar nc_std = addVariable(nc_grp, "radiance_stdev", "standard deviation of radiance in bin", "ph nm-1 s-1 sr-1 m-2", fill::f, 0.0f, 1e20f, {nc_images, nc_across_track, nc_wavelength});
    flatstd(buf);
    nc_std.putVar(buf.data());
}

void writeGeolocationData(netCDF::NcFile& nc, const std::vector<L1>& l1_products) {
    spdlog::info("Writing Geolocation Data");

    const auto n_across_track { l1_products.front().spectra.size() };
    const auto nc_across_track { nc.getDim("across_track") };
    const auto n_images { l1_products.size() };
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
        copy(l1_products[i_alt].geo.lat, lat);
        copy(l1_products[i_alt].geo.lon, lon);
        copy(l1_products[i_alt].geo.height, height);
        copy(l1_products[i_alt].geo.vza, vza);
        copy(l1_products[i_alt].geo.vaa, vaa);
        copy(l1_products[i_alt].geo.sza, sza);
        copy(l1_products[i_alt].geo.saa, saa);
    }

    auto nc_grp = nc.addGroup("geolocation_data");

    NcVar nc_lat = addVariable(nc_grp, "latitude", "latitude at bin locations", "degrees_north", fill::f, -90.0f, 90.0f, geometry_shape);
    NcVar nc_lon = addVariable(nc_grp, "longitude", "longitude at bin locations", "degrees_east", fill::f, -180.0f, 180.0f, geometry_shape);
    NcVar nc_hei = addVariable(nc_grp, "height", "height at bin locations", "m", fill::f, -1000.0f, 1000.0f, geometry_shape);
    NcVar nc_vza = addVariable(nc_grp, "viewingzenithangle", "sensor zenith angle at bin locations", "degrees", fill::f, -90.0f, 90.0f, geometry_shape);
    NcVar nc_vaa = addVariable(nc_grp, "viewingazimuthangle", "sensor azimuth angle at bin locations", "degrees", fill::f, -180.0f, 180.0f, geometry_shape);
    NcVar nc_sza = addVariable(nc_grp, "solarzenithangle", "solar zenith angle at bin locations", "degrees", fill::f, -90.0f, 90.0f, geometry_shape);
    NcVar nc_saa = addVariable(nc_grp, "solarazimuthangle", "solar azimuth angle at bin locations", "degrees", fill::f, -180.0f, 180.0f, geometry_shape);
    
    nc_lon.putVar(lat.data());
    nc_lat.putVar(lat.data());
    nc_hei.putVar(height.data());
    nc_vza.putVar(vza.data());
    nc_vaa.putVar(vaa.data());
    nc_sza.putVar(sza.data());
    nc_saa.putVar(saa.data());

}


void readL1product(const std::string& filename, const int image_start, const int image_end, std::vector<L1>& l1_products) {
    // Read netCDF
    
    spdlog::info("Reading: {}", filename);
    const netCDF::NcFile nc {filename, netCDF::NcFile::read };

    // Determine data level
    std::string level {};
    try {
        nc.getAtt("product_type").getValues(level);
    } catch (const netCDF::exceptions::NcBadId&) {
        spdlog::error("Can't find product_type variable in input data, aborting");
        std::exit(EXIT_FAILURE);
    }
    spdlog::info("Input data calibration level: {}", level);

    // Set given range of images
    const size_t alt_beg { static_cast<size_t>(image_start) };
    size_t alt_end {};
    if (image_end == fill::i){
        alt_end = nc.getDim("along_track").getSize();
    } else {
        alt_end = static_cast<size_t>(image_end) + 1;
    }
    const auto n_images { alt_end - alt_beg };
    spdlog::info("Number of images: {}", n_images);
    l1_products.resize(n_images);

    // Read the data
    readMetaData(nc, alt_beg, n_images, l1_products);
    if (level == "SGM") {
        readSceneData(nc, alt_beg, n_images, l1_products);
    } else {
        readScienceData(nc, alt_beg, n_images, l1_products);
    }
}

void readMetaData(const netCDF::NcFile& nc, const size_t alt_beg, const size_t n_images, std::vector<L1>& l1_products){
    spdlog::info("Reading Metadata");
    if (const auto grp { nc.getGroup("image_attributes") }; !grp.isNull()) {
        std::vector<double> image_time(n_images);
        std::vector<uint8_t> binning_table(n_images);
        std::vector<uint16_t> nr_coadditions(n_images);
        std::vector<double> exposure_time(n_images);
        grp.getVar("image_time").getVar({ alt_beg }, { n_images }, image_time.data());
        grp.getVar("binning_table").getVar({ alt_beg }, { n_images }, binning_table.data());
        grp.getVar("nr_coadditions").getVar({ alt_beg }, { n_images }, nr_coadditions.data());
        grp.getVar("exposure_time").getVar({ alt_beg }, { n_images }, exposure_time.data());
        for (int i_alt {}; i_alt < static_cast<int>(n_images); ++i_alt) {
            l1_products[i_alt].image_time = image_time[i_alt];
            l1_products[i_alt].binning_table_id = binning_table[i_alt];
            l1_products[i_alt].nr_coadditions = nr_coadditions[i_alt];
            l1_products[i_alt].exposure_time = exposure_time[i_alt];
        }
    }
}

void readSceneData(const netCDF::NcFile& nc, const size_t alt_beg, const size_t n_images, std::vector<L1>& l1_products){
    spdlog::info("Reading Scene Data");
    const auto n_act { nc.getDim("across_track").getSize() };
    const auto n_wavelength { nc.getDim("wavelength").getSize() };
    std::vector<double> spectra(n_images * n_act * n_wavelength);
    std::vector<double> wavelength(n_wavelength);

    nc.getVar("radiance").getVar({ alt_beg, 0, 0 }, { n_images, n_act, n_wavelength }, spectra.data());
    nc.getVar("wavelength").getVar(wavelength.data());

    // import wavelength
    auto wavelength_lbl {std::make_shared<std::vector<std::vector<double>>>(n_act, std::vector<double>(n_wavelength))};
    for (size_t i_act {}; i_act < n_act; ++i_act) {
        for (size_t i {}; i < n_wavelength; ++i) {
            (*wavelength_lbl)[i_act][i] = wavelength[i];
        }
    }

    // set first l1 wavelength, is shared between other l1 objects
    l1_products.front().wavelength = wavelength_lbl;

    //import spectra
    for (size_t i_alt {}; i_alt < n_images; ++i_alt) {
        auto& l1 { l1_products[i_alt] };
        l1.observation_sig.resize(n_act);
        for (size_t i_act {}; i_act < n_act; ++i_act) {
            l1.observation_sig[i_act].resize(n_wavelength);
            for (size_t i {}; i < n_wavelength; ++i) {
                const size_t idx { (i_alt * n_act + i_act) * n_wavelength + i };
                l1.observation_sig[i_act][i] = spectra[idx];
            }
        }
    }
}

void readScienceData(const netCDF::NcFile& nc, const size_t alt_beg, const size_t n_images, std::vector<L1>& l1_products){
    spdlog::info("Reading Science Data");
    const auto grp { nc.getGroup("science_data") };
    const auto n_bins { nc.getDim("detector_bin").getSize() }; // == rows * cols
    const auto n_rows { nc.getDim("rows").getSize() };
    const auto n_cols { nc.getDim("cols").getSize() };

    std::vector<double> detector_images(n_images * n_bins);
    std::vector<double> detector_stdev(n_images * n_bins);
    std::vector<double> flattened_wl(n_bins);

    grp.getVar("detector_image").getVar({ alt_beg, 0 }, { n_images, n_bins }, detector_images.data());
    grp.getVar("detector_stdev").getVar({ alt_beg, 0 }, { n_images, n_bins }, detector_stdev.data());
    grp.getVar("wavelength").getVar(flattened_wl.data());
    for (size_t i_alt {}; i_alt < n_images; ++i_alt) {
        auto& l1 { l1_products[i_alt] };
        l1.image.resize(n_bins);
        l1.stdev.resize(n_bins);
        for (size_t i {}; i < n_bins; ++i) {
            l1.image[i] = detector_images[i_alt * n_bins + i];
        }
    }
    
    // Import wavelength
    auto& l1 {l1_products.front()};
    std::vector<std::vector<double>> wavelength2D(n_rows, std::vector<double>(n_cols));
    for (int i {}; i < n_rows; ++i) {
        for (int j {}; j < n_cols; ++j) {
            wavelength2D[i][j] = flattened_wl[i * n_cols + j];
        }
    } 
    l1.wavelength = std::make_shared<std::vector<std::vector<double>>>(wavelength2D);
}
} // namespace tango
