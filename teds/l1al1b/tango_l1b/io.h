// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Functions for formatting the output to stdout, for manipulating
// input/output files, and tests for things like whether a file is
// readable/writable. Also, provide the main functions for reading and
// writing L1A-L1B products.

#pragma once

#include "constants.h"
#include "setting.h"
#include <netcdf>
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

auto splitString(const std::string& list,
                 const char delimiter) -> std::vector<std::string>;

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
auto copyGeometry(const std::string& l1a_filename,
                  const std::string& geo_filename,
                  int i_alt_start,
                  std::vector<L1>& l1_products) -> void;


void writeL1product(const std::string& filename, const std::string& level, const std::string& config, const std::vector<L1>& l1_products, const int argc, const char* const argv[]);
void writeGlobalAttributes(netCDF::NcFile& nc, const std::string& level, const std::string& config, const int argc, const char* const argv[]);
void writeImageAttributes(netCDF::NcFile& nc, const std::vector<L1>& l1_products);
void writeScienceData(netCDF::NcFile& nc, const std::vector<L1>& l1_products);
void writeObservationData(netCDF::NcFile& nc, const std::vector<L1>& l1_products);
void writeGeolocationData(netCDF::NcFile& nc, const std::vector<L1>& l1_products);



// Template function to get netCDF type from standard dtype
template<typename T>
netCDF::NcType getNetCDFType() {
    if constexpr (std::is_same<T, double>::value) {
        return netCDF::ncDouble;
    } else if constexpr (std::is_same<T, float>::value) {
        return netCDF::ncFloat;
    } else if constexpr (std::is_same<T, int>::value) {
        return netCDF::ncInt;
    } else if constexpr (std::is_same<T, unsigned short>::value) {
        return netCDF::ncUshort;
    } else if constexpr (std::is_same<T, char>::value) {
        return netCDF::ncByte;
    }
}

typedef netCDF::NcVar NcVar;

// Template function for adding variables to netcdf
template<typename T>
netCDF::NcVar addVariable(    
    netCDF::NcGroup& nc_grp, 
    const std::string& varname, 
    const std::string& long_name,
    const std::string& units,  
    T fill_value, 
    T min_val,
    T max_val,
    const std::vector<netCDF::NcDim>& dims) 
{
    // Use the helper function to get the appropriate netCDF type
    netCDF::NcType nc_type = getNetCDFType<T>();

    // Create the variable
    netCDF::NcVar nc_var = nc_grp.addVar(varname, nc_type, dims);

    // Set attributes
    nc_var.putAtt("_FillValue", nc_type, fill_value);
    nc_var.putAtt("long_name", long_name);
    nc_var.putAtt("units", units);
    nc_var.putAtt("valid_min", nc_type, min_val);
    nc_var.putAtt("valid_max", nc_type, max_val);

    return nc_var;
}

template<typename T>
std::vector<T> flatten3DVector(const std::vector<std::vector<std::vector<T>>>& vec3D) {
    std::vector<T> flatVec;

    // Iterate over the first dimension (outer vector)
    for (const auto& vec2D : vec3D) {
        // Iterate over the second dimension (2D vector)
        for (const auto& vec1D : vec2D) {
            // Append elements from the third dimension (1D vector)
            flatVec.insert(flatVec.end(), vec1D.begin(), vec1D.end());
        }
    }

    return flatVec;
}

template<typename T>
std::vector<T> flatten2DVector(const std::vector<std::vector<T>>& vec2D) {
    std::vector<T> flatVec;

    // Iterate over the 2D vector (each row)
    for (const auto& row : vec2D) {
        // Append elements from the current row to the flattened vector
        flatVec.insert(flatVec.end(), row.begin(), row.end());
    }

    return flatVec;
}



} // namespace tango
