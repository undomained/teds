// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "io.h"

#include "l1.h"
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

auto printHeading(const std::string& heading, const bool incl_empty_line)
  -> void
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
      "                          bitbucket.org/sron_earth/teds/issues "
      "(request permission)");
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
      "{} {:5.1f}%",
      text,
      std::min(100.0, 1e2 * iteration / static_cast<double>(work_size)));
}

auto checkPresenceOfFile(const std::string& filename, const bool required)
  -> void
{
    if (!required && filename.empty()) {
        return;
    }
    try {
        std::ifstream file { filename };
        file.exceptions(std::ifstream::failbit);
    } catch (const std::ifstream::failure& e) {
        const std::string msg { "\nCould not open file: " + filename };
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

auto splitString(const std::string& list, const char delimiter)
  -> std::vector<std::string>
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
    case ProcLevel::dark:
        return "dark";
    case ProcLevel::noise:
        return "noise";
    case ProcLevel::nonlin:
        return "nonlinearity";
    case ProcLevel::prnu:
        return "PRNU";
    case ProcLevel::stray:
        return "stray light";
    case ProcLevel::swath:
        return "swath";
    case ProcLevel::wave:
        return "spectral";
    case ProcLevel::rad:
        return "radiometric";
    case ProcLevel::l1b:
        return "L1B";
    case ProcLevel::n_levels:
    default:
        throw std::invalid_argument { "invalid process identifier" };
    }
}

auto readL1(const std::string& filename,
            const int image_start,
            const int image_end,
            std::vector<L1>& l1_products) -> void
{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };

    // Determine data level
    std::string processing_level_str {};
    try {
        nc.getAtt("processing_level").getValues(processing_level_str);
    } catch (const netCDF::exceptions::NcBadId&) {
        processing_level_str = "L1B";
    }
    ProcLevel level {};
    if (processing_level_str == "L1A") {
        level = ProcLevel::l1a;
    } else if (processing_level_str == "L1B") {
        level = ProcLevel::l1b;
    } else if (processing_level_str == "L1X") {
        ProcLevel l1x_level {};
        nc.getAtt("l1x_level").getValues(&l1x_level);
        level = l1x_level;
    }
    spdlog::info("Input data calibration level: {}", procLevelToString(level));

    // The along-track dimension is called along_track for a L1B
    // products and detector_image otherwise.
    const std::string alt_dim_name { level <= ProcLevel::stray
                                       ? "detector_image"
                                       : "along_track" };
    const size_t alt_beg { static_cast<size_t>(image_start) };
    const size_t alt_end { image_end == fill::i
                             ? nc.getDim(alt_dim_name).getSize()
                             : static_cast<size_t>(image_end + 1) };

    const auto n_images { alt_end - alt_beg };
    spdlog::info("Number of images: {}", n_images);
    l1_products.resize(n_images);
    for (L1& l1 : l1_products) {
        l1.level = level;
    }

    // Meta data
    if (const auto grp { nc.getGroup("image_attributes") }; !grp.isNull()) {
        std::vector<double> image_time(n_images);
        std::vector<uint8_t> binning_table(n_images);
        std::vector<uint16_t> nr_coadditions(n_images);
        std::vector<double> exposure_time(n_images);
        grp.getVar("image_time")
          .getVar({ alt_beg }, { n_images }, image_time.data());
        grp.getVar("binning_table")
          .getVar({ alt_beg }, { n_images }, binning_table.data());
        grp.getVar("nr_coadditions")
          .getVar({ alt_beg }, { n_images }, nr_coadditions.data());
        grp.getVar("exposure_time")
          .getVar({ alt_beg }, { n_images }, exposure_time.data());
        for (int i_alt {}; i_alt < static_cast<int>(n_images); ++i_alt) {
            l1_products[i_alt].image_time = image_time[i_alt];
            l1_products[i_alt].binning_table_id = binning_table[i_alt];
            l1_products[i_alt].nr_coadditions = nr_coadditions[i_alt];
            l1_products[i_alt].exposure_time = exposure_time[i_alt];
        }
    }

    // Science data
    if (level <= ProcLevel::stray) {
        const auto n_bins { nc.getDim("detector_bin").getSize() };
        const auto grp { nc.getGroup("science_data") };
        if (level == ProcLevel::l1a) {
            std::vector<int> detector_images(n_images * n_bins);
            grp.getVar("detector_image")
              .getVar(
                { alt_beg, 0 }, { n_images, n_bins }, detector_images.data());
            for (size_t i_alt {}; i_alt < n_images; ++i_alt) {
                auto& l1 { l1_products[i_alt] };
                l1.image.resize(n_bins);
                l1.stdev.resize(n_bins);
                for (size_t i {}; i < n_bins; ++i) {
                    l1.image[i] =
                      static_cast<double>(detector_images[i_alt * n_bins + i])
                      / l1.nr_coadditions;
                }
            }
        } else if (level <= ProcLevel::stray) {
            std::vector<double> detector_images(n_images * n_bins);
            std::vector<double> detector_stdev(n_images * n_bins);
            grp.getVar("detector_image")
              .getVar(
                { alt_beg, 0 }, { n_images, n_bins }, detector_images.data());
            grp.getVar("detector_stdev")
              .getVar(
                { alt_beg, 0 }, { n_images, n_bins }, detector_stdev.data());
            for (size_t i_alt {}; i_alt < n_images; ++i_alt) {
                auto& l1 { l1_products[i_alt] };
                l1.image.resize(n_bins);
                l1.stdev.resize(n_bins);
                for (size_t i {}; i < n_bins; ++i) {
                    l1.image[i] = detector_images[i_alt * n_bins + i];
                }
            }
        }
    }
    if (level > ProcLevel::stray) {
        const auto n_act { nc.getDim("across_track").getSize() };
        const auto n_wavelength { nc.getDim("wavelength").getSize() };
        std::vector<double> spectra(n_images * n_act * n_wavelength);
        std::vector<double> spectra_stdev(n_images * n_act * n_wavelength);
        std::vector<double> wavelength(n_act * n_wavelength);
        const auto nc_grp { nc.getGroup("observation_data") };
        nc_grp.getVar("i").getVar(
          { alt_beg, 0, 0 }, { n_images, n_act, n_wavelength }, spectra.data());
        nc_grp.getVar("i_stdev").getVar({ alt_beg, 0, 0 },
                                        { n_images, n_act, n_wavelength },
                                        spectra_stdev.data());
        nc_grp.getVar("wavelength").getVar(wavelength.data());
        auto wavelength_lbl {
            std::make_shared<std::vector<std::vector<double>>>(
              n_act, std::vector<double>(n_wavelength))
        };
        for (size_t i_act {}; i_act < n_act; ++i_act) {
            for (size_t i {}; i < n_wavelength; ++i) {
                (*wavelength_lbl)[i_act][i] =
                  wavelength[i_act * n_wavelength + i];
            }
        }
        for (size_t i_alt {}; i_alt < n_images; ++i_alt) {
            auto& l1 { l1_products[i_alt] };
            l1.spectra.resize(n_act);
            l1.wavelength = wavelength_lbl;
            for (size_t i_act {}; i_act < n_act; ++i_act) {
                l1.spectra[i_act].signal.resize(n_wavelength);
                l1.spectra[i_act].stdev.resize(n_wavelength);
                for (size_t i {}; i < n_wavelength; ++i) {
                    const size_t idx { (i_alt * n_act + i_act) * n_wavelength
                                       + i };
                    l1.spectra[i_act].signal[i] = spectra[idx];
                    l1.spectra[i_act].stdev[i] = spectra_stdev[idx];
                }
            }
        }
    }
}

auto writeL1(const std::string& filename,
             const std::string& config,
             const std::vector<L1>& l1_products,
             const int argc,
             const char* const argv[]) -> void
{
    // L1B or lower level product
    netCDF::NcFile nc { filename, netCDF::NcFile::replace };
    const auto& level { l1_products.front().level };
    spdlog::info("Output data calibration level: {}", procLevelToString(level));

    // Global attributes
    nc.putAtt("Conventions", "CF-1.11");
    if (level == ProcLevel::l1b) {
        nc.putAtt("title", "Tango Carbon level 1B data");
    } else {
        nc.putAtt("title", "Tango Carbon level 1A data");
    }
    if (level == ProcLevel::l1a || level == ProcLevel::l1b) {
        nc.putAtt("processing_level", procLevelToString(level));
    } else {
        nc.putAtt("processing_level", "L1X");
        nc.putAtt("l1x_level", netCDF::ncInt, static_cast<int>(level));
    }
    nc.putAtt("project", "TANGO");
    nc.putAtt("instrument", "TANGO");
    nc.putAtt("product_name", filename);
    nc.putAtt("creator_name", "SRON/Earth Science");
    nc.putAtt("creator_url", "https://www.sron.nl/missions-earth");
    nc.putAtt("date_created", getDateAndTime());
    if (TANGO_GIT_COMMIT_ABBREV != "GITDIR-N") {
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

    // Along track dimension
    const auto n_images { l1_products.size() };
    const auto nc_images { nc.addDim(
      (level <= ProcLevel::stray ? "detector_image" : "along_track"),
      n_images) };

    auto nc_var { nc.addVar("configuration", netCDF::ncString) };
    nc_var.putAtt("comment",
                  "configuration parameters used for producing this file");
    const char* conf_char { config.c_str() };
    nc_var.putVar(&conf_char);

    // Image attributes
    if (level < ProcLevel::l1b) {
        auto nc_grp { nc.addGroup("image_attributes") };
        nc_var = nc_grp.addVar("image_time", netCDF::ncDouble, { nc_images });
        nc_var.putAtt("name", "image time");
        nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::d);
        nc_var.putAtt("units", "seconds since 2022-03-21");
        nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
        nc_var.putAtt("valid_max", netCDF::ncDouble, 92304.0);
        std::vector<double> buf {};
        buf.assign(n_images, 0.0);
        nc_var.putVar(buf.data());

        nc_var = nc_grp.addVar("binning_table", netCDF::ncByte, { nc_images });
        nc_var.putAtt("name", "binning table");
        buf.assign(n_images, l1_products.front().binning_table_id);
        nc_var.putVar(buf.data());

        nc_var =
          nc_grp.addVar("nr_coadditions", netCDF::ncUshort, { nc_images });
        nc_var.putAtt("name", "number of coadditions");
        buf.assign(n_images, l1_products.front().nr_coadditions);
        nc_var.putVar(buf.data());

        nc_var =
          nc_grp.addVar("exposure_time", netCDF::ncDouble, { nc_images });
        nc_var.putAtt("name", "exposure time");
        buf.assign(n_images, l1_products.front().exposure_time);
        nc_var.putVar(buf.data());
    }

    // Science or observation data
    if (level <= ProcLevel::stray) {
        auto nc_grp { nc.addGroup("science_data") };
        const auto n_bins { level == ProcLevel::l1a
                              ? l1_products.front().image_i32.size()
                              : l1_products.front().image.size() };
        const auto nc_detector_bin { nc.addDim("detector_bin", n_bins) };
        const std::vector<netCDF::NcDim> nc_detector_shape { nc_images,
                                                             nc_detector_bin };
        if (level == ProcLevel::l1a) {
            nc_var =
              nc_grp.addVar("detector_image", netCDF::ncInt, nc_detector_shape);
            nc_var.putAtt("name", "detector images");
            nc_var.putAtt("_FillValue", netCDF::ncInt, fill::i);
            nc_var.putAtt("units", "counts");
            nc_var.putAtt("valid_min", netCDF::ncInt, 0);
            nc_var.putAtt("valid_max", netCDF::ncInt, 60000);
            nc_var.setCompression(true, true, compression_level);
            std::vector<size_t> chunksize { 1, n_bins };
            nc_var.setChunking(netCDF::NcVar::nc_CHUNKED, chunksize);
            std::vector<int> buf_i32(n_images * n_bins);
            for (size_t i {}; i < n_images; ++i) {
                for (size_t j {}; j < n_bins; ++j) {
                    buf_i32[i * n_bins + j] = l1_products[i].image_i32[j];
                }
            }
            nc_var.putVar(buf_i32.data());
        } else {
            nc_var = nc_grp.addVar(
              "detector_image", netCDF::ncDouble, nc_detector_shape);
            nc_var.putAtt("name", "detector images");
            nc_var.putAtt("units", "counts");
            nc_var.putAtt("valid_min", netCDF::ncDouble, -1e100);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e100);
            std::vector<double> buf(n_images * n_bins);
            for (size_t i {}; i < n_images; ++i) {
                for (size_t j {}; j < n_bins; ++j) {
                    buf[i * n_bins + j] = l1_products[i].image[j];
                }
            }
            nc_var.putVar(buf.data());
            nc_var = nc_grp.addVar(
              "detector_stdev", netCDF::ncDouble, nc_detector_shape);
            nc_var.putAtt("name", "standard deviation of detector bin");
            nc_var.putAtt("units", "counts");
            nc_var.putAtt("valid_min", netCDF::ncDouble, 0.0);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e100);
            std::ranges::fill(buf, 1.0);
            nc_var.putVar(buf.data());
        }
    } else {
        const auto n_across_track { l1_products.front().spectra.size() };
        const auto wavelengths { *l1_products.front().wavelength };
        const auto n_wavelength { wavelengths.front().size() };
        const auto nc_across_track { nc.addDim("across_track",
                                               n_across_track) };
        const auto nc_wavelength { nc.addDim("wavelength", n_wavelength) };

        // Observation data. The difference between L1B and lower
        // levels is mostly in the data format - float vs double.
        auto nc_grp { nc.addGroup("observation_data") };

        nc_var = nc_grp.addVar(
          "wavelength", netCDF::ncFloat, { nc_across_track, nc_wavelength });
        nc_var.putAtt("name", "radiance wavelengths");
        nc_var.putAtt("_FillValue", netCDF::ncFloat, fill::f);
        nc_var.putAtt("valid_min", netCDF::ncFloat, 0.0f);
        nc_var.putAtt("valid_max", netCDF::ncFloat, 999.0f);
        nc_var.putAtt("units", "nm");
        std::vector<float> buf(n_across_track * n_wavelength);
        for (size_t i {}; i < n_across_track; ++i) {
            for (size_t j {}; j < n_wavelength; ++j) {
                buf[i * n_wavelength + j] =
                  static_cast<float>(wavelengths[i][j]);
            }
        }
        nc_var.putVar(buf.data());

        const auto flattenSignal { [&](auto& buf) {
            for (size_t i {}; i < n_images; ++i) {
                for (size_t j {}; j < n_across_track; ++j) {
                    for (size_t k {}; k < n_wavelength; ++k) {
                        buf[(i * n_across_track + j) * n_wavelength + k] =
                          static_cast<float>(
                            l1_products[i].spectra[j].signal[k]);
                    }
                }
            }
        } };
        const auto flattenStdev { [&](auto& buf) {
            for (size_t i {}; i < n_images; ++i) {
                for (size_t j {}; j < n_across_track; ++j) {
                    for (size_t k {}; k < n_wavelength; ++k) {
                        buf[(i * n_across_track + j) * n_wavelength + k] =
                          static_cast<float>(
                            l1_products[i].spectra[j].stdev[k]);
                    }
                }
            }
        } };
        if (level == ProcLevel::l1b) {
            nc_var =
              nc_grp.addVar("i",
                            netCDF::ncFloat,
                            { nc_images, nc_across_track, nc_wavelength });
            nc_var.putAtt("name", "radiance");
            nc_var.putAtt("_FillValue", netCDF::ncFloat, fill::f);
            nc_var.putAtt("units", "ph nm-1 s-1 sr-1 m-2");
            nc_var.putAtt("valid_min", netCDF::ncFloat, 0.0f);
            nc_var.putAtt("valid_max", netCDF::ncFloat, 1e20f);
            std::vector<float> buf(n_images * n_across_track * n_wavelength);
            flattenSignal(buf);
            nc_var.putVar(buf.data());
            nc_var =
              nc_grp.addVar("i_stdev",
                            netCDF::ncFloat,
                            { nc_images, nc_across_track, nc_wavelength });
            nc_var.putAtt("name", "standard deviation of radiance in bin ");
            nc_var.putAtt("_FillValue", netCDF::ncFloat, fill::f);
            nc_var.putAtt("units", "ph nm-1 s-1 sr-1 m-2");
            nc_var.putAtt("valid_min", netCDF::ncFloat, 0.0f);
            nc_var.putAtt("valid_max", netCDF::ncFloat, 1e20f);
            flattenStdev(buf);
            nc_var.putVar(buf.data());
        } else {
            nc_var =
              nc_grp.addVar("i",
                            netCDF::ncDouble,
                            { nc_images, nc_across_track, nc_wavelength });
            nc_var.putAtt("name", "radiance");
            nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::f);
            nc_var.putAtt("units", "ph nm-1 s-1 sr-1 m-2");
            nc_var.putAtt("valid_min", netCDF::ncDouble, -1e20);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e20);
            std::vector<double> buf(n_images * n_across_track * n_wavelength);
            flattenSignal(buf);
            nc_var.putVar(buf.data());
            nc_var =
              nc_grp.addVar("i_stdev",
                            netCDF::ncDouble,
                            { nc_images, nc_across_track, nc_wavelength });
            nc_var.putAtt("name", "standard deviation of radiance in bin ");
            nc_var.putAtt("_FillValue", netCDF::ncDouble, fill::f);
            nc_var.putAtt("units", "ph nm-1 s-1 sr-1 m-2");
            nc_var.putAtt("valid_min", netCDF::ncDouble, -1e20);
            nc_var.putAtt("valid_max", netCDF::ncDouble, 1e20);
            flattenStdev(buf);
            nc_var.putVar(buf.data());
        }
    }
}

} // namespace tango
