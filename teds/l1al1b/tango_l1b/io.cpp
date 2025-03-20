// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "io.h"

#include "cubic_spline.h"
#include "geometry.h"
#include "l1.h"
#include "quaternion.h"
#include "time.h"

#include <filesystem>
#include <fstream>
#include <netcdf>
#include <omp.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

namespace tango {

auto tango_formatter_flag::format(const spdlog::details::log_msg& log_msg,
                                  const std::tm& /* tm_time */,
                                  spdlog::memory_buf_t& dest) -> void
{
    std::string text {};
    switch (log_msg.level) {
    case spdlog::level::info:
        break;
    case spdlog::level::warn:
        text = " [warning]";
        break;
    case spdlog::level::err:
        text = " [error]";
        break;
    case spdlog::level::debug:
        text = " [debug]";
        break;
    case spdlog::level::off:
    case spdlog::level::trace:
    case spdlog::level::critical:
    case spdlog::level::n_levels:
    default:
        text = " [unknown]";
    }
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    dest.append(text.data(), text.data() + text.size());
}

auto tango_formatter_flag::clone() const
  -> std::unique_ptr<custom_flag_formatter>
{
    return spdlog::details::make_unique<tango_formatter_flag>();
}

auto initLogging(const bool set_debug_level) -> void
{
    // Only if the logger does not already exist
    if (!spdlog::get("plain")) {
        auto formatter { std::make_unique<spdlog::pattern_formatter>() };
        formatter->add_flag<tango_formatter_flag>('*').set_pattern(
          "[%H:%M:%S]%* %v");
        spdlog::set_formatter(std::move(formatter));
        spdlog::stdout_color_mt("plain");
        spdlog::get("plain")->set_pattern("%v");
    }
    if (set_debug_level) {
        spdlog::set_level(spdlog::level::debug);
    } else {
        spdlog::set_level(spdlog::level::info);
    }
}

auto printHeading(const std::string& heading,
                  const bool incl_empty_line) -> void
{
    if (incl_empty_line) {
        spdlog::get("plain")->info("");
    }
    std::string hash_line(heading.size() + 4, '#');
    spdlog::get("plain")->info(hash_line);
    spdlog::get("plain")->info("# " + heading + " #");
    spdlog::get("plain")->info(hash_line);
}

auto printSystemInfo(const std::string& project_version,
                     const std::string& git_commit,
                     const std::string& cmake_host_system,
                     const std::string& executable,
                     const std::string& compiler,
                     const std::string& compiler_flags,
                     const std::string& libraries,
                     const std::string& binning_table) -> void
{
    spdlog::get("plain")->info("Version                 : {}", project_version);
    if (git_commit != "GITDIR-N") {
        spdlog::get("plain")->info("Commit hash             : {}", git_commit);
    }
    spdlog::get("plain")->info("Date and timezone       : {}", getDate());
    spdlog::get("plain")->info(
      "Contacts                : raullaasner@gmail.com\n"
      "                          bitbucket.org/sron_earth/teds/issues");
    spdlog::get("plain")->info("Host system             : {}",
                               cmake_host_system);
    spdlog::get("plain")->info("Executable location     : {}", executable);
    spdlog::get("plain")->info("C++ compiler            : {}", compiler);
    spdlog::get("plain")->info("C++ compiler flags      : {}", compiler_flags);
    spdlog::get("plain")->info("Number of threads       : {}",
                               omp_get_max_threads());
    for (bool first_line { true };
         const auto& lib : splitString(libraries, ' ')) {
        if (first_line) {
            spdlog::get("plain")->info("Linking against         : {}", lib);
            first_line = false;
        } else {
            spdlog::get("plain")->info("                          {}", lib);
        }
    }
    spdlog::get("plain")->info("Binning table file      : {}", binning_table);
}

auto checkPresenceOfFile(const Setting<std::string>& setting,
                         const bool required) -> void
{
    if (!required && setting.empty()) {
        return;
    }
    try {
        std::ifstream file { setting };
        file.exceptions(std::ifstream::failbit);
    } catch (const std::ifstream::failure& e) {
        if (setting.empty()) {
            std::string yaml_keys_str {};
            for (const auto& yaml_key : setting.yaml_keys) {
                yaml_keys_str += '[' + yaml_key + ']';
            }
            const std::string msg { "missing " + yaml_keys_str };
            throw std::runtime_error { msg };
        }
        const std::string msg { "\nCould not open file: " + setting };
        throw std::ifstream::failure { e.what() + msg, e.code() };
    }
}

auto checkFileWritable(const std::string& filename) -> void
{
    if (!filename.empty()) {
        std::string dir { std::filesystem::path { filename }.parent_path() };
        if (!dir.empty() && static_cast<bool>(access(dir.c_str(), W_OK))) {
            throw std::runtime_error { filename + " is not writable" };
        }
    }
}

auto splitString(const std::string& list,
                 const char delimiter) -> std::vector<std::string>
{
    std::stringstream ss { list };
    std::string name {};
    std::vector<std::string> strings {};
    while (getline(ss, name, delimiter)) {
        strings.push_back(name);
    }
    return strings;
}

[[nodiscard]] auto procLevelToString(const ProcLevel proc_level) -> std::string
{
    switch (proc_level) {
    case ProcLevel::l1a:
        return "L1A";
    case ProcLevel::raw:
        return "raw";
    case ProcLevel::dark_offset:
        return "dark offset";
    case ProcLevel::noise:
        return "noise";
    case ProcLevel::dark_current:
        return "dark current";
    case ProcLevel::nonlin:
        return "nonlinearity";
    case ProcLevel::prnu:
        return "PRNU";
    case ProcLevel::stray:
        return "stray light";
    case ProcLevel::swath:
        return "swath";
    case ProcLevel::l1b:
        return "L1B";
    case ProcLevel::sgm:
        return "SGM";
    case ProcLevel::n_levels:
    default:
        throw std::invalid_argument { "invalid process identifier" };
    }
}

auto readL1(const std::string& filename,
            const size_t alt_beg,
            const std::optional<size_t> alt_end,
            L1& l1_prod,
            const bool in_memory) -> void
{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };

    // Determine data level
    std::string product_type_str {};
    try {
        nc.getAtt("product_type").getValues(product_type_str);
    } catch (const netCDF::exceptions::NcBadId&) {
        product_type_str = "L1B";
    }
    if (product_type_str == "L1A") {
        l1_prod.level = ProcLevel::l1a;
    } else if (product_type_str == "L1B") {
        l1_prod.level = ProcLevel::l1b;
    } else if (product_type_str == "SGM") {
        l1_prod.level = ProcLevel::sgm;
    } else if (product_type_str == "L1X") {
        nc.getAtt("l1x_level").getValues(&l1_prod.level);
    }
    spdlog::info("Input data calibration level: {}",
                 procLevelToString(l1_prod.level));

    const size_t n_alt { alt_end.value_or(
                           nc.getDim("along_track_sample").getSize() - 1)
                         + 1 - alt_beg };
    l1_prod.n_alt = static_cast<int>(n_alt);
    spdlog::info("Number of along-track bins: {}", n_alt);

    // Meta data
    if (const auto grp { nc.getGroup("image_attributes") }; !grp.isNull()) {
        l1_prod.time.resize(n_alt);
        l1_prod.tai_seconds.resize(n_alt);
        l1_prod.tai_subsec.resize(n_alt);
        grp.getVar("time").getVar({ alt_beg }, { n_alt }, l1_prod.time.data());
        grp.getVar("tai_seconds")
          .getVar({ alt_beg }, { n_alt }, l1_prod.tai_seconds.data());
        std::vector<uint16_t> tai_subsec(n_alt);
        grp.getVar("tai_subsec")
          .getVar({ alt_beg }, { n_alt }, tai_subsec.data());
        for (int i {}; i < l1_prod.n_alt; ++i) {
            l1_prod.tai_subsec[i] = tai_subsec[i] / 65535.0;
        }
        grp.getVar("binning_table").getVar(&l1_prod.binning_table_id);
        grp.getVar("nr_coadditions").getVar(&l1_prod.nr_coadditions);
        grp.getVar("exposure_time").getVar(&l1_prod.exposure_time);
    }

    // Science data
    if (!nc.getGroup("science_data").isNull()) {
        const auto grp { nc.getGroup("science_data") };
        const auto n_bins { nc.getDim("bin").getSize() };
        l1_prod.signal.resize(n_alt * n_bins);
        l1_prod.noise.resize(n_alt * n_bins, 1.0);
        grp.getVar("detector_image")
          .getVar({ alt_beg, 0 }, { n_alt, n_bins }, l1_prod.signal.data());
        if (const auto var { grp.getVar("detector_stdev") }; !var.isNull()) {
            var.getVar({ alt_beg, 0 }, { n_alt, n_bins }, l1_prod.noise.data());
        }
    } else {
        const auto n_act { nc.getDim("across_track_sample").getSize() };
        const auto n_wavelength { nc.getDim("wavelength").getSize() };
        l1_prod.wavelengths.resize(n_act, std::vector<double>(n_wavelength));
        if (l1_prod.level == ProcLevel::sgm
            || nc.getGroup("observation_data").isNull()) {
            std::vector<double> wavelengths(n_wavelength);
            nc.getVar("wavelength").getVar(wavelengths.data());
            for (size_t i_act {}; i_act < n_act; ++i_act) {
                std::copy(wavelengths.begin(),
                          wavelengths.begin() + n_wavelength,
                          l1_prod.wavelengths[i_act].begin());
            }
            if (in_memory) {
                l1_prod.spectra.resize(n_alt * n_act * n_wavelength);
                nc.getVar("radiance")
                  .getVar({ alt_beg, 0, 0 },
                          { n_alt, n_act, n_wavelength },
                          l1_prod.spectra.data());
            }
        } else {
            const auto grp { nc.getGroup("observation_data") };
            std::vector<double> wavelengths(n_act * n_wavelength);
            grp.getVar("wavelength").getVar(wavelengths.data());
            for (size_t i_act {}; i_act < n_act; ++i_act) {
                std::copy(wavelengths.begin() + i_act * n_wavelength,
                          wavelengths.begin() + (i_act + 1) * n_wavelength,
                          l1_prod.wavelengths[i_act].begin());
            }
            l1_prod.spectra.resize(n_alt * n_act * n_wavelength);
            grp.getVar("radiance")
              .getVar({ alt_beg, 0, 0 },
                      { n_alt, n_act, n_wavelength },
                      l1_prod.spectra.data());
            l1_prod.spectra_noise.assign(n_alt * n_act * n_wavelength, 1.0);
            if (const auto var { grp.getVar("radiance_stdev") };
                !var.isNull()) {
                var.getVar({ alt_beg, 0, 0 },
                           { n_alt, n_act, n_wavelength },
                           l1_prod.spectra_noise.data());
            }
        }
    }

    // Navigation data
    if (const auto grp { nc.getGroup("navigation_data") }; !grp.isNull()) {
        const auto n_time { grp.getDim("time").getSize() };
        std::vector<double> times(n_time);
        std::vector<double> orb_pos(n_time * dims::vec);
        std::vector<double> att_quat_raw(n_time * dims::quat);
        grp.getVar("time").getVar(times.data());
        grp.getVar("orb_pos").getVar(orb_pos.data());
        grp.getVar("att_quat").getVar(att_quat_raw.data());
        // Interpolate orbit positions to detector image timestamps
        l1_prod.orb_pos.resize(l1_prod.n_alt * dims::vec);
        for (int i_dir {}; i_dir < dims::vec; ++i_dir) {
            // Orbit position component (x/y/z)
            std::vector<double> pos_comp(n_time);
            for (int i {}; i < static_cast<int>(n_time); ++i) {
                pos_comp[i] = orb_pos[i * dims::vec + i_dir];
            }
            const CubicSpline spline { times, pos_comp };
            for (int i_alt {}; i_alt < l1_prod.n_alt; ++i_alt) {
                l1_prod.orb_pos[i_alt * dims::vec + i_dir] =
                  spline.eval(l1_prod.time[i_alt]);
            }
        }
        // Interpolate attitude quaternions to detector image timestamps
        std::vector<Quaternion> att_quat(n_time);
        for (int i {}; i < static_cast<int>(n_time); ++i) {
            att_quat[i] = { att_quat_raw[i * dims::quat + 0],
                            att_quat_raw[i * dims::quat + 1],
                            att_quat_raw[i * dims::quat + 2],
                            att_quat_raw[i * dims::quat + 3] };
            att_quat[i].normalize();
        }
        l1_prod.att_quat.resize(l1_prod.n_alt);
        for (int i_alt {}; i_alt < l1_prod.n_alt; ++i_alt) {
            interpolateQuaternion(
              times, att_quat, l1_prod.time[i_alt], l1_prod.att_quat[i_alt]);
        }
    }
}

auto writeL1(const std::string& filename,
             const std::string& config,
             const L1& l1_prod,
             const int argc,
             const char* const argv[]) -> void
{
    // L1B or lower level product
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    spdlog::info("Writing output data (calibration level: {})",
                 procLevelToString(l1_prod.level));

    // Global attributes
    nc.putAtt("Conventions", "CF-1.11");
    if (l1_prod.level == ProcLevel::l1b) {
        nc.putAtt("title", "Tango Carbon level 1B data");
    } else {
        nc.putAtt("title", "Tango Carbon level 1A data");
    }
    if (l1_prod.level == ProcLevel::l1a || l1_prod.level == ProcLevel::l1b) {
        nc.putAtt("product_type", procLevelToString(l1_prod.level));
    } else {
        nc.putAtt("product_type", "L1X");
        nc.putAtt("l1x_level", netCDF::ncInt, static_cast<int>(l1_prod.level));
    }
    nc.putAtt("project", "TANGO");
    nc.putAtt("instrument", "Carbon");
    nc.putAtt("product_name", filename);
    nc.putAtt("creator_name", "SRON/Earth Science");
    nc.putAtt("creator_url", "https://earth.sron.nl/project/tango");
    nc.putAtt("date_created", getDateAndTime());
    const std::string no_git_commit { "GITDIR-N" };
    if (TANGO_GIT_COMMIT_ABBREV != no_git_commit) {
        nc.putAtt("git_commit", TANGO_GIT_COMMIT_ABBREV);
    }
    std::string command_line { argc > 0 ? argv[0] : "" };
    for (int i { 1 }; i < argc; ++i) {
        command_line = command_line + ' ' + argv[i];
    }
    nc.putAtt("history", command_line);
    // Extract processing version from the input file
    const std::string processing_version {
        YAML::Load(config)["processing_version"].as<std::string>()
    };
    nc.putAtt("processing_version", processing_version);

    // Compression will be enabled only for official products (L1A, L1B)
    constexpr int compression_level { 5 };

    const size_t n_alt { static_cast<size_t>(l1_prod.n_alt) };
    const auto nc_alt { nc.addDim("along_track_sample", n_alt) };

    auto nc_var { nc.addVar("configuration", netCDF::ncString) };
    nc_var.putAtt("comment",
                  "configuration parameters used for producing this file");
    const char* conf_char { config.c_str() };
    nc_var.putVar(&conf_char);

    // Detector image attributes
    if (l1_prod.level < ProcLevel::l1b) {
        // Image timestamps come from the geometry file
        const netCDF::NcFile nc_geo {
            YAML::Load(config)["io_files"]["geometry"].as<std::string>(),
            netCDF::NcFile::read
        };
        std::string time_units {};
        nc_geo.getVar("time").getAtt("units").getValues(time_units);
        const size_t alt_beg { YAML::Load(config)["alt_beg"].as<size_t>() };

        auto nc_grp { nc.addGroup("image_attributes") };

        nc_var = nc_grp.addVar("time", netCDF::ncDouble, { nc_alt });
        nc_var.putAtt("long_name", "detector image time");
        nc_var.putAtt("description",
                      "integration start time in seconds of day");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("units", time_units);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 172800.0); // 2 x day
        std::vector<double> time(n_alt);
        nc_geo.getVar("time").getVar({ alt_beg }, { n_alt }, time.data());
        nc_var.putVar(time.data());

        nc_var = nc_grp.addVar("tai_seconds", netCDF::ncUint, { nc_alt });
        nc_var.putAtt("long_name", "detector image TAI time (seconds)");
        nc_var.putAtt("units", "seconds since 1958-01-01 00:00:00 TAI");
        nc_var.putAtt("valid_min", netCDF::ncUint, 1956528000u);
        nc_var.putAtt("valid_max", netCDF::ncUint, 2493072000u); // 2 x day
        std::vector<uint32_t> tai_seconds(n_alt);
        nc_geo.getVar("tai_seconds")
          .getVar({ alt_beg }, { n_alt }, tai_seconds.data());
        nc_var.putVar(tai_seconds.data());

        nc_var = nc_grp.addVar("tai_subsec", netCDF::ncUshort, { nc_alt });
        nc_var.putAtt("long_name", "detector image TAI time (subseconds)");
        nc_var.putAtt("units", "1/65536 s");
        nc_var.putAtt("valid_min", netCDF::ncUint, 0);
        nc_var.putAtt("valid_max", netCDF::ncUint, 65535);
        std::vector<uint16_t> tai_subsec(n_alt);
        nc_geo.getVar("tai_subsec")
          .getVar({ alt_beg }, { n_alt }, tai_subsec.data());
        nc_var.putVar(tai_subsec.data());

        nc_var = nc_grp.addVar("binning_table", netCDF::ncByte);
        nc_var.putAtt("long_name", "binning table ID");
        nc_var.putVar(&l1_prod.binning_table_id);

        nc_var = nc_grp.addVar("nr_coadditions", netCDF::ncUshort);
        nc_var.putAtt("long_name", "coaddition factor");
        nc_var.putAtt("comment", "number of detector read-outs summed");
        nc_var.putVar(&l1_prod.nr_coadditions);

        nc_var = nc_grp.addVar("exposure_time", netCDF::ncDouble);
        nc_var.putAtt("long_name", "exposure time");
        nc_var.putAtt("comment", "exposure time per detector read-out");
        nc_var.putAtt("units", "s");
        nc_var.putVar(&l1_prod.exposure_time);
    }

    // Science or observation data
    if (l1_prod.level <= ProcLevel::stray) {
        auto nc_grp { nc.addGroup("science_data") };
        const auto n_bins { l1_prod.signal.size() / n_alt };
        const auto nc_bin { nc.addDim("bin", n_bins) };
        const std::vector<netCDF::NcDim> nc_detector_shape { nc_alt, nc_bin };
        if (l1_prod.level == ProcLevel::l1a) {
            nc_var =
              nc_grp.addVar("detector_image", netCDF::ncInt, nc_detector_shape);
            nc_var.putAtt("long_name", "detector image signal");
            nc_var.putAtt("_FillValue", netCDF::ncInt, fill::i);
            nc_var.putAtt("units", "counts");
            nc_var.putAtt("valid_min", netCDF::ncInt, 0);
            nc_var.putAtt("valid_max", netCDF::ncInt, 60000);
            nc_var.setCompression(true, true, compression_level);
            std::vector<size_t> chunksize { 1, n_bins };
            nc_var.setChunking(netCDF::NcVar::nc_CHUNKED, chunksize);
            nc_var.putVar(l1_prod.signal.data());
        } else {
            nc_var = nc_grp.addVar(
              "detector_image", netCDF::ncDouble, nc_detector_shape);
            nc_var.putAtt("long_name", "detector image signal");
            nc_var.putAtt("units", "counts");
            nc_var.putAtt("valid_min", netCDF::ncDouble, -1e10);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e30);
            nc_var.putVar(l1_prod.signal.data());
        }
        if (!l1_prod.noise.empty()) {
            nc_var = nc_grp.addVar(
              "detector_stdev", netCDF::ncDouble, nc_detector_shape);
            nc_var.putAtt("long_name", "standard deviation of detector bin");
            nc_var.putAtt("units", "counts");
            nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e100);
            nc_var.putVar(l1_prod.noise.data());
        }
    } else {
        auto nc_grp { nc.addGroup("observation_data") };

        const auto n_act { l1_prod.wavelengths.size() };
        const auto n_wavelength { l1_prod.wavelengths.front().size() };
        const auto nc_act { nc.addDim("across_track_sample", n_act) };
        const auto nc_wavelength { nc.addDim("wavelength", n_wavelength) };

        nc_var = nc_grp.addVar(
          "wavelength", netCDF::ncDouble, { nc_act, nc_wavelength });
        nc_var.putAtt("long_name", "wavelength");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 2000.0);
        nc_var.putAtt("units", "nm");
        std::vector<double> buf(n_act * n_wavelength);
        for (size_t i {}; i < n_act; ++i) {
            for (size_t j {}; j < n_wavelength; ++j) {
                buf[i * n_wavelength + j] = l1_prod.wavelengths[i][j];
            }
        }
        nc_var.putVar(buf.data());

        nc_var = nc_grp.addVar(
          "radiance", netCDF::ncDouble, { nc_alt, nc_act, nc_wavelength });
        nc_var.putAtt("long_name", "spectral photon radiance");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("units", "nm-1 s-1 sr-1 m-2");
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1e20);
        nc_var.putVar(l1_prod.spectra.data());

        if (!l1_prod.spectra_noise.empty()) {
            nc_var = nc_grp.addVar("radiance_stdev",
                                   netCDF::ncDouble,
                                   { nc_alt, nc_act, nc_wavelength });
            nc_var.putAtt("long_name",
                          "standard deviation of radiance in bin ");
            nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
            nc_var.putAtt("units", "nm-1 s-1 sr-1 m-2");
            nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e20);
            nc_var.putVar(l1_prod.spectra_noise.data());
        }
    }
    // Geolocation data
    if (l1_prod.geo.lat.empty()) {
        return;
    }
    // Need a copy of geometry because l1_prod is read-only but we
    // need to convert the angles to degrees.
    Geometry geo { l1_prod.geo };
    const auto n_act { l1_prod.geo.lat.size() / n_alt };
    const auto nc_act { l1_prod.level >= ProcLevel::swath
                          ? nc.getDim("across_track_sample")
                          : nc.addDim("across_track_sample", n_act) };
    auto nc_grp { nc.addGroup("geolocation_data") };
    const std::vector<netCDF::NcDim> geometry_shape { nc_alt, nc_act };

    nc_var = nc_grp.addVar("latitude", netCDF::ncDouble, geometry_shape);
    nc_var.putAtt("long_name", "latitude at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -90.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 90.0);
    nc_var.putAtt("units", "degrees_north");
    moduloAndConvertAngles(geo.lat);
    nc_var.putVar(geo.lat.data());

    nc_var = nc_grp.addVar("longitude", netCDF::ncDouble, geometry_shape);
    nc_var.putAtt("long_name", "longitude at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -180.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 180.0);
    nc_var.putAtt("units", "degrees_east");
    moduloAndConvertAngles(geo.lon);
    nc_var.putVar(geo.lon.data());

    nc_var = nc_grp.addVar("height", netCDF::ncDouble, geometry_shape);
    nc_var.putAtt("long_name", "height from sea level at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -1000.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 10000.0);
    nc_var.putAtt("units", "m");
    nc_var.putVar(geo.height.data());

    nc_var = nc_grp.addVar("sensor_zenith", netCDF::ncDouble, geometry_shape);
    nc_var.putAtt("long_name", "sensor zenith angle at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -90.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 90.0);
    nc_var.putAtt("units", "degrees");
    moduloAndConvertAngles(geo.vza);
    nc_var.putVar(geo.vza.data());

    nc_var = nc_grp.addVar("sensor_azimuth", netCDF::ncDouble, geometry_shape);
    nc_var.putAtt("long_name", "sensor azimuth angle at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -180.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 180.0);
    nc_var.putAtt("units", "degrees");
    moduloAndConvertAngles(geo.vaa);
    nc_var.putVar(geo.vaa.data());

    nc_var = nc_grp.addVar("solar_zenith", netCDF::ncDouble, geometry_shape);
    nc_var.putAtt("long_name", "solar zenith angle at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -90.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 90.0);
    nc_var.putAtt("units", "degrees");
    moduloAndConvertAngles(geo.sza);
    nc_var.putVar(geo.sza.data());

    nc_var = nc_grp.addVar("solar_azimuth", netCDF::ncDouble, geometry_shape);
    nc_var.putAtt("long_name", "solar azimuth angle at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -180.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 180.0);
    nc_var.putAtt("units", "degrees");
    moduloAndConvertAngles(geo.saa);
    nc_var.putVar(geo.saa.data());
}

// Fetch the image_start value that was used in the instrument model
static auto getIMImageStart(const std::string& l1a_filename) -> size_t
{
    char* buf[10'000];
    netCDF::NcFile { l1a_filename, netCDF::NcFile::read }
      .getVar("configuration")
      .getVar(buf);
    YAML::Node config = YAML::Load(*buf);
    return config["alt_beg"].as<size_t>();
}

auto copyGeometry(const std::string& l1a_filename,
                  const std::string& geo_filename,
                  const size_t alt_beg,
                  L1& l1_prod) -> void
{
    const netCDF::NcFile nc_geo { geo_filename, netCDF::NcFile::read };
    const auto n_act { nc_geo.getDim("across_track_sample").getSize() };
    l1_prod.geo.lat.resize(l1_prod.n_alt * n_act);
    l1_prod.geo.lon.resize(l1_prod.n_alt * n_act);
    l1_prod.geo.vza.resize(l1_prod.n_alt * n_act);
    l1_prod.geo.vaa.resize(l1_prod.n_alt * n_act);
    l1_prod.geo.sza.resize(l1_prod.n_alt * n_act);
    l1_prod.geo.saa.resize(l1_prod.n_alt * n_act);
    // Starts and counts
    std::vector<size_t> s { alt_beg + getIMImageStart(l1a_filename), 0 };
    std::vector<size_t> c { static_cast<size_t>(l1_prod.n_alt), n_act };
    nc_geo.getVar("latitude").getVar(s, c, l1_prod.geo.lat.data());
    nc_geo.getVar("longitude").getVar(s, c, l1_prod.geo.lon.data());
    nc_geo.getVar("sensor_zenith").getVar(s, c, l1_prod.geo.vza.data());
    nc_geo.getVar("sensor_azimuth").getVar(s, c, l1_prod.geo.vaa.data());
    nc_geo.getVar("solar_zenith").getVar(s, c, l1_prod.geo.sza.data());
    nc_geo.getVar("solar_azimuth").getVar(s, c, l1_prod.geo.saa.data());
    for (int i {}; i < static_cast<int>(l1_prod.geo.lat.size()); ++i) {
        l1_prod.geo.lat[i] *= math::deg_to_rad;
        l1_prod.geo.lon[i] *= math::deg_to_rad;
        l1_prod.geo.vza[i] *= math::deg_to_rad;
        l1_prod.geo.vaa[i] *= math::deg_to_rad;
        l1_prod.geo.sza[i] *= math::deg_to_rad;
        l1_prod.geo.saa[i] *= math::deg_to_rad;
    }
    l1_prod.geo.height.assign(l1_prod.n_alt * n_act, 0.0);
}

auto copyNavigationData(const std::string& navigation_filename,
                        const std::string& l1a_filename) -> void
{
    const netCDF::NcFile nc_nav { navigation_filename, netCDF::NcFile::read };
    netCDF::NcFile nc_l1a { l1a_filename, netCDF::NcFile::write };
    auto grp { nc_l1a.addGroup("navigation_data") };
    const auto n_time { nc_nav.getDim("time").getSize() };
    const auto nc_time_dim { grp.addDim("time", n_time) };
    const auto nc_vec_dim { grp.addDim("vector_elements", 3) };
    const auto nc_quat_dim { grp.addDim("quaternion_elements", 4) };

    const auto copyAtts { [](const netCDF::NcVar& in, netCDF::NcVar& out) {
        void* value_buf[32768] {};
        for (const auto& [name, att] : in.getAtts()) {
            att.getValues(value_buf);
            out.putAtt(name, att.getType(), att.getAttLength(), value_buf);
        }
    } };
    netCDF::NcVar var {};
    std::vector<double> buf(n_time * dims::quat);

    var = grp.addVar("time", netCDF::ncDouble, nc_time_dim);
    copyAtts(nc_nav.getVar("time"), var);
    nc_nav.getVar("time").getVar(buf.data());
    var.putVar(buf.data());

    var = grp.addVar("orb_pos", netCDF::ncDouble, { nc_time_dim, nc_vec_dim });
    copyAtts(nc_nav.getVar("orb_pos"), var);
    nc_nav.getVar("orb_pos").getVar(buf.data());
    var.putVar(buf.data());

    var =
      grp.addVar("att_quat", netCDF::ncDouble, { nc_time_dim, nc_quat_dim });
    copyAtts(nc_nav.getVar("att_quat"), var);
    nc_nav.getVar("att_quat").getVar(buf.data());
    var.putVar(buf.data());
}

} // namespace tango
