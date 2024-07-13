// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Functions for formatting the output to stdout, for manipulating
// input/output files, and tests for things like whether a file is
// readable/writable. Also, provide the main functions for reading and
// writing L1A-L1B products.

#pragma once

#include "constants.h"
#include "setting.h"

#include <spdlog/pattern_formatter.h>

namespace tango {

class L1;

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
// just the message text). The argument sets the logger level to
// either debug or info.
auto initLogging(const bool set_debug_level) -> void;

// Print the name of a major section at its beginning. For example,
// the CKD initialization would start with
//
// ######################
// # CKD initialization #
// ######################
auto printHeading(const std::string& heading, const bool incl_empty_line = true)
  -> void;

// Print information about the host system and how the executable or
// library was built and some runtime options.
auto printSystemInfo(const std::string& project_version,
                     const std::string& git_commit,
                     const std::string& cmake_host_system,
                     const std::string& executable,
                     const std::string& compiler,
                     const std::string& compiler_flags,
                     const std::string& libraries,
                     const std::string& binning_table) -> void;

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

auto splitString(const std::string& list, const char delimiter)
  -> std::vector<std::string>;

// Convert a process level enum to string suitable for displaying in
// output.
[[nodiscard]] auto procLevelToString(const ProcLevel proc_level) -> std::string;

// Read a list of L1 products from a single NetCDF file. The input
// data level may be L1A, L1B, or anything in between. image_start/end
// specify a subrange to process. To process all images, use values 0
// and fill::i.
auto readL1(const std::string& filename,
            const int image_start,
            const int image_end,
            std::vector<L1>& l1_products) -> void;

// Write a L1 product to file. The NetCDF structure of the product
// depends on the data product level (L1A-L1B). argc and argv are used
// to record how the program was invoked.
auto writeL1(const std::string& filename,
             const std::string& config,
             const std::vector<L1>& l1_products,
             const int argc = 0,
             const char* const argv[] = nullptr) -> void;

// Copy geolocation data from the geometry file directly to the L1B
// product. This is a placeholder function until geolocation is
// properly implemented.
auto copyGeometry(const std::string& filename,
                  const int i_alt_start,
                  std::vector<L1>& l1_products) -> void;

} // namespace tango
