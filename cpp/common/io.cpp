// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "io.h"

#include "cubic_spline.h"
#include "geometry.h"
#include "l1.h"
#include "l2.h"
#include "time.h"

#include <filesystem>
#include <fstream>
#include <netcdf>
#include <omp.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>

namespace tango {

// Compression will be enabled only for official products (L1A, L1B, L2)
constexpr int compression_level { 5 };

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

auto initLogging() -> void
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
    spdlog::set_level(spdlog::level::info);
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
                     const std::string& libraries) -> void
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
}

auto printPercentage(const int iteration,
                     const size_t work_size,
                     const std::string_view text) -> void
{
    if (omp_get_thread_num() != 0) {
        return;
    }
    // Number of times this function has been called by the first
    // process. The actual iteration number should always be larger
    // than the adjusted iteration number. If not then the adjusted
    // iteration number is from a previous parallel for loop and needs
    // to be reset.
    static int adjusted_iteration { fill::i };
    static size_t prev_work_size {};
    if (iteration < adjusted_iteration || iteration < omp_get_num_threads()
        || work_size != prev_work_size) {
        adjusted_iteration = fill::i;
        prev_work_size = work_size;
    }
    ++adjusted_iteration;
    spdlog::info(
      "{} {:6.2f}%",
      text,
      std::min(100.0, 1e2 * iteration / static_cast<double>(work_size)));
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

auto lower(const std::string& str) -> std::string
{
    std::string str_l { str };
    for (size_t i {}; i < str_l.length(); ++i) {
        str_l[i] = std::tolower(str_l[i]);
    }
    return str_l;
}

auto upper(const std::string& str) -> std::string
{
    std::string str_u { str };
    for (size_t i {}; i < str_u.length(); ++i) {
        str_u[i] = std::toupper(str_u[i]);
    }
    return str_u;
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
    case ProcLevel::l2:
        return "L2";
    case ProcLevel::n_levels:
    default:
        throw std::invalid_argument { "invalid process identifier" };
    }
}

auto readGeometry(const netCDF::NcGroup& grp,
                  const size_t n_alt,
                  const size_t n_act,
                  const size_t alt_beg,
                  Geometry& geo) -> void
{
    geo.lat.resize(n_alt, n_act);
    geo.lon.resize(n_alt, n_act);
    geo.height.resize(n_alt, n_act);
    geo.vza.resize(n_alt, n_act);
    geo.vaa.resize(n_alt, n_act);
    geo.sza.resize(n_alt, n_act);
    geo.saa.resize(n_alt, n_act);
    // Starts and counts
    std::vector<size_t> starts { alt_beg, 0 };
    std::vector<size_t> counts { n_alt, n_act };
    grp.getVar("latitude").getVar(starts, counts, geo.lat.data());
    grp.getVar("longitude").getVar(starts, counts, geo.lon.data());
    grp.getVar("sensor_zenith").getVar(starts, counts, geo.vza.data());
    grp.getVar("sensor_azimuth").getVar(starts, counts, geo.vaa.data());
    grp.getVar("solar_zenith").getVar(starts, counts, geo.sza.data());
    grp.getVar("solar_azimuth").getVar(starts, counts, geo.saa.data());
    geo.lat *= math::deg_to_rad;
    geo.lon *= math::deg_to_rad;
    geo.vza *= math::deg_to_rad;
    geo.vaa *= math::deg_to_rad;
    geo.sza *= math::deg_to_rad;
    geo.saa *= math::deg_to_rad;
    if (const auto var { grp.getVar("height") }; var.isNull()) {
        geo.height = 0.0;
    } else {
        var.getVar(starts, counts, geo.height.data());
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

    const size_t nc_n_alt { nc.getDim("along_track_sample").getSize() };
    const size_t n_alt { std::min(alt_end.value_or(nc_n_alt - 1) + 1 - alt_beg,
                                  nc_n_alt) };
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
        l1_prod.signal.resize(n_alt, n_bins);
        l1_prod.noise.resize(n_alt, n_bins);
        grp.getVar("detector_image")
          .getVar({ alt_beg, 0 }, { n_alt, n_bins }, l1_prod.signal.data());
        if (const auto var { grp.getVar("detector_stdev") }; var.isNull()) {
            l1_prod.noise = 1.0;
        } else {
            var.getVar({ alt_beg, 0 }, { n_alt, n_bins }, l1_prod.noise.data());
        }
    } else {
        const size_t n_act { nc.getDim("across_track_sample").getSize() };
        const auto n_wavelength { nc.getDim("wavelength").getSize() };
        l1_prod.wavelengths.resize(n_wavelength);
        if (const auto grp { nc.getGroup("observation_data") }; grp.isNull()) {
            nc.getVar("wavelength").getVar(l1_prod.wavelengths.data());
            if (in_memory || l1_prod.level < ProcLevel::sgm) {
                l1_prod.spectra.resize(n_alt * n_act, n_wavelength);
                nc.getVar("radiance")
                  .getVar({ alt_beg, 0, 0 },
                          { n_alt, n_act, n_wavelength },
                          l1_prod.spectra.data());
            }
        } else {
            grp.getVar("wavelength").getVar(l1_prod.wavelengths.data());
            l1_prod.spectra.resize(n_alt * n_act, n_wavelength);
            grp.getVar("radiance")
              .getVar({ alt_beg, 0, 0 },
                      { n_alt, n_act, n_wavelength },
                      l1_prod.spectra.data());
            l1_prod.spectra_noise.resize(n_alt * n_act, n_wavelength);
            if (const auto var { grp.getVar("radiance_stdev") }; var.isNull()) {
                l1_prod.spectra_noise = 1.0;
            } else {
                var.getVar({ alt_beg, 0, 0 },
                           { n_alt, n_act, n_wavelength },
                           l1_prod.spectra_noise.data());
            }
        }
    }

    // Geometry
    if (const auto grp { nc.getGroup("geolocation_data") }; !grp.isNull()) {
        const size_t n_act { nc.getDim("across_track_sample").getSize() };
        readGeometry(grp, n_alt, n_act, alt_beg, l1_prod.geo);
    }

    // Navigation data
    if (const auto grp { nc.getGroup("navigation_data") }; !grp.isNull()) {
        const auto n_time { grp.getDim("time").getSize() };
        Eigen::ArrayXd times(n_time);
        ArrayXXd orb_pos(n_time, dims::vec);
        ArrayXXd att_quat_raw(n_time, dims::quat);
        grp.getVar("time").getVar(times.data());
        grp.getVar("orb_pos").getVar(orb_pos.data());
        grp.getVar("att_quat").getVar(att_quat_raw.data());
        // Interpolate orbit positions to detector image timestamps
        l1_prod.orb_pos0.resize(l1_prod.n_alt * dims::vec);
        l1_prod.orb_pos.resize(l1_prod.n_alt, dims::vec);
        for (int i_dir {}; i_dir < dims::vec; ++i_dir) {
            const CubicSpline spline { times, orb_pos.col(i_dir) };
            l1_prod.orb_pos.col(i_dir) = spline.eval(l1_prod.time);
            for (int i_alt {}; i_alt < l1_prod.n_alt; ++i_alt) {
                l1_prod.orb_pos0[i_alt * dims::vec + i_dir] =
                  l1_prod.orb_pos(i_alt, i_dir);
            }
        }
        // Interpolate attitude quaternions to detector image timestamps
        std::vector<Eigen::Quaterniond> att_quat(n_time);
        for (int i {}; i < static_cast<int>(n_time); ++i) {
            att_quat[i] = Eigen::Quaterniond(att_quat_raw(i, 3),
                                             att_quat_raw(i, 0),
                                             att_quat_raw(i, 1),
                                             att_quat_raw(i, 2))
                            .normalized();
        }
        l1_prod.att_quat.resize(l1_prod.n_alt);
        for (int i_alt {}; i_alt < l1_prod.n_alt; ++i_alt) {
            l1_prod.att_quat[i_alt] =
              interpolateQuaternion(times, att_quat, l1_prod.time[i_alt]);
        }
    }
}

// Apply compression to a NetCDF variable
static auto setCompression(const bool compress,
                           netCDF::NcVar& nc_var,
                           std::vector<size_t> chunksizes) -> void
{
    if (compress) {
        nc_var.setCompression(true, true, compression_level);
        nc_var.setChunking(netCDF::NcVar::nc_CHUNKED, chunksizes);
    }
}

// Write global attributes to a NetCDF file
static auto writeHeader(netCDF::NcFile& nc,
                        const ProcLevel level,
                        const std::string& filename,
                        const std::string& config,
                        const int argc,
                        const char* const argv[]) -> void
{
    // Global attributes
    nc.putAtt("Conventions", "CF-1.11");
    if (level == ProcLevel::l2) {
        nc.putAtt("title", "Tango Carbon level 2 data");
    } else if (level == ProcLevel::l1b) {
        nc.putAtt("title", "Tango Carbon level 1B data");
    } else {
        nc.putAtt("title", "Tango Carbon level 1A data");
    }
    if (level == ProcLevel::l1a || level == ProcLevel::l1b
        || level == ProcLevel::l2) {
        nc.putAtt("product_type", procLevelToString(level));
    } else {
        nc.putAtt("product_type", "L1X");
        nc.putAtt("l1x_level", netCDF::ncInt, static_cast<int>(level));
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

    auto nc_var { nc.addVar("configuration", netCDF::ncString) };
    nc_var.putAtt("comment",
                  "configuration parameters used for producing this file");
    const char* conf_char { config.c_str() };
    nc_var.putVar(&conf_char);
}

// Write latitudes to a NetCDF group. Values are converted to degrees
// first and brought to the range -180..180.
static auto writeLat(netCDF::NcGroup& nc_grp,
                     const std::vector<netCDF::NcDim>& geometry_shape,
                     ArrayXXd& lat) -> void
{
    auto nc_var { nc_grp.addVar("latitude", netCDF::ncDouble, geometry_shape) };
    nc_var.putAtt("long_name", "latitude at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -90.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 90.0);
    nc_var.putAtt("units", "degrees_north");
    moduloAndConvertAngles(lat);
    nc_var.putVar(lat.data());
}

// Write longitudes to a NetCDF group
static auto writeLon(netCDF::NcGroup& nc_grp,
                     const std::vector<netCDF::NcDim>& geometry_shape,
                     ArrayXXd& lon) -> void
{
    auto nc_var { nc_grp.addVar(
      "longitude", netCDF::ncDouble, geometry_shape) };
    nc_var.putAtt("long_name", "longitude at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -180.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 180.0);
    nc_var.putAtt("units", "degrees_east");
    moduloAndConvertAngles(lon);
    nc_var.putVar(lon.data());
}

// Write heights (elevation from sea level) to a NetCDF group
static auto writeHeight(netCDF::NcGroup& nc_grp,
                        const std::vector<netCDF::NcDim>& geometry_shape,
                        const ArrayXXd& height) -> void
{
    auto nc_var { nc_grp.addVar("height", netCDF::ncDouble, geometry_shape) };
    nc_var.putAtt("long_name", "height from sea level at bin locations");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, -1000.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 10000.0);
    nc_var.putAtt("units", "m");
    nc_var.putVar(height.data());
}

auto writeL1(const std::string& filename,
             const std::string& config,
             const L1& l1_prod,
             const bool compress,
             const int argc,
             const char* const argv[]) -> void
{
    // L1B or lower level product
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    spdlog::info("Writing output data (calibration level: {})",
                 procLevelToString(l1_prod.level));

    writeHeader(nc, l1_prod.level, filename, config, argc, argv);

    const size_t n_alt { static_cast<size_t>(l1_prod.n_alt) };
    const auto nc_alt { nc.addDim("along_track_sample", n_alt) };
    netCDF::NcVar nc_var {};

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
        const auto n_bins { static_cast<size_t>(l1_prod.signal.cols()) };
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
            setCompression(compress, nc_var, { 1, n_bins });
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
        if (l1_prod.noise.size() > 0) {
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

        const size_t n_act { static_cast<size_t>(l1_prod.spectra.rows()
                                                 / l1_prod.n_alt) };
        const size_t n_wavelength { static_cast<size_t>(
          l1_prod.wavelengths.size()) };
        const auto nc_act { nc.addDim("across_track_sample", n_act) };
        const auto nc_wavelength { nc.addDim("wavelength", n_wavelength) };

        nc_var = nc_grp.addVar("wavelength", netCDF::ncDouble, nc_wavelength);
        nc_var.putAtt("long_name", "wavelength");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 1550.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1700.0);
        nc_var.putAtt("units", "nm");
        nc_var.putVar(l1_prod.wavelengths.data());

        nc_var = nc_grp.addVar(
          "radiance", netCDF::ncDouble, { nc_alt, nc_act, nc_wavelength });
        nc_var.putAtt("long_name", "spectral photon radiance");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("units", "nm-1 s-1 sr-1 m-2");
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1e20);
        std::vector<size_t> chunksizes { 1, n_act, n_wavelength };
        setCompression(compress, nc_var, chunksizes);
        nc_var.putVar(l1_prod.spectra.data());

        if (l1_prod.spectra_noise.size() > 0) {
            nc_var = nc_grp.addVar("radiance_stdev",
                                   netCDF::ncDouble,
                                   { nc_alt, nc_act, nc_wavelength });
            nc_var.putAtt("long_name",
                          "standard deviation of radiance in bin ");
            nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
            nc_var.putAtt("units", "nm-1 s-1 sr-1 m-2");
            nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e20);
            setCompression(compress, nc_var, chunksizes);
            nc_var.putVar(l1_prod.spectra_noise.data());
        }
    }
    // Geolocation data
    if (l1_prod.geo.lat.size() == 0) {
        return;
    }
    // Need a copy of geometry because l1_prod is read-only but we
    // need to convert the angles to degrees.
    Geometry geo { l1_prod.geo };
    const auto n_act { l1_prod.geo.lat.cols() };
    const auto nc_act { l1_prod.level >= ProcLevel::swath
                          ? nc.getDim("across_track_sample")
                          : nc.addDim("across_track_sample", n_act) };
    auto nc_grp { nc.addGroup("geolocation_data") };
    const std::vector<netCDF::NcDim> geometry_shape { nc_alt, nc_act };

    writeLat(nc_grp, geometry_shape, geo.lat);

    writeLon(nc_grp, geometry_shape, geo.lon);

    writeHeight(nc_grp, geometry_shape, geo.height);

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

auto writeL2(const std::string& filename,
             const std::string& config,
             const L2& l2,
             Geometry& geometry,
             const Eigen::ArrayXd& z_lay,
             const std::map<std::string, Eigen::ArrayXd>& ref_profiles,
             const bool compress,
             const int argc,
             const char* const argv[]) -> void
{
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };

    writeHeader(nc, ProcLevel::l2, filename, config, argc, argv);

    const auto nc_alt { nc.addDim("along_track_sample", geometry.lat.rows()) };
    const auto nc_act { nc.addDim("across_track_sample", geometry.lat.cols()) };
    const auto nc_lay { nc.addDim("layers",
                                  l2.col_avg_kernels.at("CO2").cols()) };

    // Work variables
    netCDF::NcGroup nc_grp {};
    netCDF::NcVar nc_var {};
    std::vector<double> buf {};
    std::vector<size_t> chunksizes { nc_alt.getSize(), nc_act.getSize() };

    const std::map<std::string, double> scales { { "CO2", 1e6 },
                                                 { "CH4", 1e9 },
                                                 { "H2O", 1e6 } };

    nc_var = nc.addVar("zlay", netCDF::ncDouble, nc_lay);
    nc_var.putAtt("long_name", "central layer height");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 1e+5);
    nc_var.putAtt("units", "m");
    nc_var.putVar(z_lay.data());

    for (const std::string gas : { "CO2", "CH4" }) {
        const std::string xgas { 'x' + lower(gas) };
        const std::string xgas_proxy { 'X' + gas + " proxy" };

        nc_var =
          nc.addVar(xgas + "_proxy", netCDF::ncDouble, { nc_alt, nc_act });
        setCompression(compress, nc_var, chunksizes);
        nc_var.putAtt("long_name", xgas_proxy + " dry air column mixing ratio");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 2e20);
        nc_var.putAtt("units", gas == "CH4" ? "ppbv" : "ppmv");
        nc_var.putVar((l2.proxys.at(gas) * scales.at(gas)).eval().data());

        nc_var = nc.addVar(
          "precision_" + xgas + "_proxy", netCDF::ncDouble, { nc_alt, nc_act });
        setCompression(compress, nc_var, chunksizes);
        nc_var.putAtt("long_name",
                      xgas_proxy + " dry air column mixing ratio precision");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 2e20);
        nc_var.putAtt("units", gas == "CH4" ? "ppbv" : "ppmv");
        nc_var.putVar(
          (l2.proxy_precisions.at(gas) * scales.at(gas)).eval().data());

        nc_var = nc.addVar(
          "accuracy_" + xgas + "_proxy", netCDF::ncDouble, { nc_alt, nc_act });
        setCompression(compress, nc_var, chunksizes);
        nc_var.putAtt("long_name",
                      xgas_proxy + " dry air column mixing ratio accuracy");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 300.0);
        nc_var.putAtt("units", gas == "CH4" ? "ppbv" : "ppmv");
        buf.assign(nc_alt.getSize() * nc_act.getSize(), 99999.0);
        nc_var.putVar(buf.data());

        nc_var = nc.addVar(
          "qa_value_" + xgas + "_proxy", netCDF::ncDouble, { nc_alt, nc_act });
        setCompression(compress, nc_var, chunksizes);
        nc_var.putAtt("long_name",
                      xgas_proxy
                        + " dry air column mixing ratio quality value");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 300.0);
        buf.assign(nc_alt.getSize() * nc_act.getSize(), 100.0);
        nc_var.putVar(buf.data());
    }

    for (const std::string gas : { "CO2", "CH4", "H2O" }) {
        const std::string xgas { 'x' + lower(gas) };
        nc_var = nc.addVar("col_avg_kernel_" + xgas,
                           netCDF::ncDouble,
                           { nc_alt, nc_act, nc_lay });
        setCompression(
          compress, nc_var, { 1, nc_act.getSize(), nc_lay.getSize() });
        nc_var.putAtt("long_name",
                      lower(gas) + " column averaging kernel of " + xgas);
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1e28);
        nc_var.putAtt("units", "molecules / cm2");
        nc_var.putVar(l2.col_avg_kernels.at(gas).data());
    }

    nc_var = nc.addVar("aquisition_time", netCDF::ncDouble, nc_alt);
    nc_var.putAtt("long_name", "aquisition time of sensing in UTC");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 1e+25);
    nc_var.putAtt("units", "s");
    buf.assign(nc_alt.getSize(), 99999.0);
    nc_var.putVar(buf.data());

    nc_var = nc.addVar("albedo", netCDF::ncDouble, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "wavelength-independent component of albedo");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 1.0);
    nc_var.putAtt("units", "s");
    nc_var.putVar(l2.albedo0.data());

    nc_var = nc.addVar("spectral_shift", netCDF::ncDouble, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "constant shift of spectral calibration");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 10.0);
    nc_var.putAtt("units", "nm");
    nc_var.putVar(l2.spec_shift.data());

    nc_var =
      nc.addVar("spectral_squeeze", netCDF::ncDouble, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "squeeze of spectral calibration");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 10.0);
    nc_var.putVar(l2.spec_squeeze.data());

    nc_grp = nc.addGroup("non_scattering_retrieval");
    for (const std::string gas : { "CO2", "CH4", "H2O" }) {
        const std::string xgas { 'x' + lower(gas) };
        nc_var = nc_grp.addVar(xgas, netCDF::ncDouble, { nc_alt, nc_act });
        setCompression(compress, nc_var, chunksizes);
        nc_var.putAtt("long_name",
                      lower(gas) + " dry air column mixing ratio " + xgas);
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1e4);
        nc_var.putAtt("units", gas == "CH4" ? "ppbv" : "ppmv");
        nc_var.putVar(
          (l2.mixing_ratios.at(gas) * scales.at(gas)).eval().data());

        nc_var = nc_grp.addVar(
          "precision_" + xgas, netCDF::ncDouble, { nc_alt, nc_act });
        setCompression(compress, nc_var, chunksizes);
        nc_var.putAtt("long_name",
                      lower(gas) + " precision of dry air column mixing ratio "
                        + xgas);
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1e4);
        nc_var.putAtt("units", gas == "CH4" ? "ppbv" : "ppmv");
        nc_var.putVar((l2.precisions.at(gas) * scales.at(gas)).eval().data());
    }

    nc_grp = nc.addGroup("geolocation_data");
    writeLat(nc_grp, { nc_alt, nc_act }, geometry.lat);
    writeLon(nc_grp, { nc_alt, nc_act }, geometry.lon);
    writeHeight(nc_grp, { nc_alt, nc_act }, geometry.height);

    nc_grp = nc.addGroup("prior");
    for (const std::string gas : { "CO2", "CH4", "H2O" }) {
        nc_var =
          nc_grp.addVar("apri_prof_" + lower(gas), netCDF::ncDouble, nc_lay);
        nc_var.putAtt("long_name",
                      "a priori profile " + lower(gas)
                        + " in layer column density");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1e28);
        nc_var.putAtt("units", "molecules / cm2");
        nc_var.putVar(ref_profiles.at(gas).data());
    }
    nc_var =
      nc_grp.addVar("surface_pressure", netCDF::ncDouble, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "surface pressure");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 1400.0);
    buf.assign(nc_alt.getSize() * nc_act.getSize(), 1013.0);
    nc_var.putVar(buf.data());

    nc_grp = nc.addGroup("diagnostics");

    nc_var =
      nc_grp.addVar("processing_flag", netCDF::ncShort, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name",
                  "processing flag to indicate processor anomalies");
    nc_var.putAtt("_FillValue", netCDF::ncShort, fill::i);
    nc_var.putAtt("valid_min", netCDF::ncShort, 0);
    nc_var.putAtt("valid_max", netCDF::ncShort, 101);
    std::vector<int16_t> buf_i16(nc_alt.getSize() * nc_act.getSize(), 100);
    nc_var.putVar(buf_i16.data());

    nc_var = nc_grp.addVar("convergence", netCDF::ncByte, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "whether pixel converged");
    nc_var.putAtt("_FillValue", netCDF::ncByte, 7);
    nc_var.putAtt("valid_min", netCDF::ncByte, 0);
    nc_var.putAtt("valid_max", netCDF::ncByte, 1);
    std::vector<uint8_t> buf_u8(nc_alt.getSize() * nc_act.getSize());
    for (int i {}; i < static_cast<int>(buf_u8.size()); ++i) {
        buf_u8[i] =
          static_cast<bool>(l2.converged.reshaped<Eigen::RowMajor>()(i));
    }
    nc_var.putVar(buf_u8.data());

    nc_var = nc_grp.addVar("iterations", netCDF::ncInt, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "number of iterations");
    nc_var.putAtt("_FillValue", netCDF::ncInt, fill::i);
    nc_var.putAtt("valid_min", netCDF::ncInt, 0);
    nc_var.putAtt("valid_max", netCDF::ncInt, 100);
    nc_var.putVar(l2.iterations.data());

    nc_var = nc_grp.addVar("chi2", netCDF::ncDouble, { nc_alt, nc_act });
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "spectral chi square of the fit");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::i);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 100.0);
    nc_var.putVar(l2.chi2.data());
}

auto writeL2Diagnostics(const std::string& filename,
                        const std::string& config,
                        const L1& l1b,
                        const L2& l2,
                        const bool compress,
                        const int argc,
                        const char* const argv[]) -> void
{
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };

    writeHeader(nc, ProcLevel::l2, filename, config, argc, argv);
    nc.putAtt("title", "Tango Carbon level 2 diagnostics data");

    const auto nc_alt { nc.addDim("along_track_sample", l1b.n_alt) };
    const auto nc_act { nc.addDim("across_track_sample",
                                  l1b.spectra.rows() / l1b.n_alt) };
    const auto nc_wav { nc.addDim("wavelength", l1b.wavelengths.size()) };

    // Work variables
    netCDF::NcVar nc_var {};
    const std::vector<netCDF::NcDim> nc_dims { nc_alt, nc_act, nc_wav };
    std::vector<size_t> chunksizes { 1, nc_act.getSize(), nc_wav.getSize() };

    nc_var = nc.addVar("wavelength", netCDF::ncDouble, nc_wav);
    nc_var.putAtt("long_name", "wavelengths");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 8000.0);
    nc_var.putAtt("units", "nm");
    nc_var.putVar(l1b.wavelengths.data());

    nc_var = nc.addVar("measurement", netCDF::ncDouble, nc_dims);
    setCompression(compress, nc_var, chunksizes);
    nc_var.putAtt("long_name", "spectral photon radiance");
    nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
    nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
    nc_var.putAtt("valid_max", netCDF::ncDouble, 1e20);
    nc_var.putAtt("units", "nm-1 s-1 sr-1 m-2");
    nc_var.putVar(l1b.spectra.data());

    for (const std::string gas : { "CO2", "CH4", "H2O" }) {
        nc_var = nc.addVar("gain_" + lower(gas), netCDF::ncDouble, nc_dims);
        setCompression(compress, nc_var, chunksizes);
        nc_var.putAtt("long_name", gas + " spectral gain vector");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, 1e37);
        nc_var.putAtt("valid_min", netCDF::ncDouble, -1e30);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 1e30);
        nc_var.putAtt("units", "ppm/(photons/(nm m2 s sr))");
        nc_var.putVar(l2.gains.at(gas).data());
    }
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
    readGeometry(nc_geo,
                 static_cast<size_t>(l1_prod.n_alt),
                 n_act,
                 alt_beg + getIMImageStart(l1a_filename),
                 l1_prod.geo);
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
