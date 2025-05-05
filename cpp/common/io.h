// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Functions for formatting the output to stdout, for manipulating
// input/output files, and tests for things like whether a file is
// readable/writable. Also, provide the main functions for reading and
// writing L1A-L1B products.

#pragma once

#include "constants.h"
#include "setting.h"

#include <Eigen/Dense>
#include <map>
#include <spdlog/pattern_formatter.h>

namespace netCDF {

class NcGroup;

} // namespace netCDF

namespace tango {

class Geometry;
class L1;
class L2;

// Define a new spdlog formatter flag. The primary purpose is to show
// labels such as [warning] for warnings but no label for regular
// (info) messages.
class tango_formatter_flag : public spdlog::custom_flag_formatter
{
public:
    auto format(const spdlog::details::log_msg& log_msg,
                const std::tm&,
                spdlog::memory_buf_t& dest) -> void override;
    auto clone() const -> std::unique_ptr<custom_flag_formatter> override;
};

// Set up two loggers: one with a verbose pattern (mostly used
// throughout the code) and a plain one (clean pattern, i.e. prints
// just the message text).
auto initLogging() -> void;

// Print the name of a processing section. For example, the CKD
// initialization would start with
//
// ######################
// # CKD initialization #
// ######################
auto printHeading(const std::string& heading,
                  const bool incl_empty_line = true) -> void;

// Print information about the host system and how the executable or
// library was built and some runtime options.
auto printSystemInfo(const std::string& project_version,
                     const std::string& git_commit,
                     const std::string& cmake_host_system,
                     const std::string& executable,
                     const std::string& compiler,
                     const std::string& compiler_flags,
                     const std::string& libraries) -> void;

// Print the percentage of work done (iteration / work_size). Call
// this in OpenMP parallel for loops with dynamic scheduling.
auto printPercentage(const int iteration,
                     const size_t work_size,
                     const std::string_view text) -> void;

// Check if filename exists. If required == true and filename is an
// empty string raise an error, otherwise return without checking.
auto checkPresenceOfFile(const Setting<std::string>& setting,
                         const bool required) -> void;

// Check if destination is writable
auto checkFileWritable(const std::string& filename) -> void;

auto splitString(const std::string& list,
                 const char delimiter) -> std::vector<std::string>;

// Convert word to lower or upper case
auto lower(const std::string& str) -> std::string;
auto upper(const std::string& str) -> std::string;

// Convert a process level enum to string suitable for displaying in
// output.
[[nodiscard]] auto procLevelToString(const ProcLevel proc_level) -> std::string;

auto readGeometry(const netCDF::NcGroup& grp,
                  const size_t n_alt,
                  const size_t n_act,
                  const size_t alt_beg,
                  Geometry& geo) -> void;

// Read a list of L1 products from a single NetCDF file. The input
// data level may be L1A, L1B, or anything in between. alt_beg/end
// specify a subrange of along-track positions to process. To process
// all use values 0 and fill::i. If the input is SGM spectra, memory
// usage can be reduced with in_meomry=false. The spectra are then not
// read into the L1 product and are only read later one by one for the
// ISRF convolution.
auto readL1(const std::string& filename,
            const size_t alt_beg,
            const std::optional<size_t> alt_end,
            L1& l1_prod,
            const bool in_memory = false) -> void;

// Write a L1 product to file. The NetCDF structure of the product
// depends on the data product level (L1A-L1B). argc and argv are used
// to record how the program was invoked.
auto writeL1(const std::string& filename,
             const std::string& config,
             const L1& l1_prod,
             const bool compress,
             const int argc = 0,
             const char* const argv[] = nullptr) -> void;

// Write L2 product to file
auto writeL2(const std::string& filename,
             const std::string& config,
             const L2& l2,
             Geometry& geometry,
             const Eigen::ArrayXd& z_lay,
             const std::map<std::string, Eigen::ArrayXd>& ref_profiles,
             const bool compress,
             const int argc,
             const char* const argv[] = nullptr) -> void;

// Write additional L2 diagnostics information to a separate file
auto writeL2Diagnostics(const std::string& filename,
                        const std::string& config,
                        const L1& l1b,
                        const L2& l2,
                        const bool compress,
                        const int argc,
                        const char* const argv[]) -> void;

// Copy geolocation data from the geometry file directly to the L1B
// product. This is a placeholder function until geolocation is
// properly implemented.
auto copyGeometry(const std::string& l1a_filename,
                  const std::string& geo_filename,
                  const size_t alt_beg,
                  L1& l1_prod) -> void;

auto copyNavigationData(const std::string& navigation_filename,
                        const std::string& l1a_filename) -> void;

} // namespace tango
