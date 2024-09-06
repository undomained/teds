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
#include "settings_l1b.h"

namespace tango {

class L1;

// Write functions
void writeL1product(const std::string& filename, const std::string& level, const std::string& config, const std::vector<L1>& l1_products, const int argc, const char* const argv[]);
void writeGlobalAttributes(netCDF::NcFile& nc, const std::string& level, const std::string& config, const int argc, const char* const argv[]);
void writeMetaData(netCDF::NcFile& nc, const std::vector<L1>& l1_products);
void writeScienceData(netCDF::NcFile& nc, const std::vector<L1>& l1_products);
void writeObservationData(netCDF::NcFile& nc, const std::vector<L1>& l1_products);
void writeGeolocationData(netCDF::NcFile& nc, const std::vector<L1>& l1_products);

// Read functions
void readL1product(const std::string& filename, const int image_start, const int image_end, std::vector<L1>& l1_products);
void readMetaData(const netCDF::NcFile& nc, const size_t alt_beg, const size_t n_images, std::vector<L1>& l1_products);
void readScienceData(const netCDF::NcFile& nc, const size_t alt_beg, const size_t n_images, std::vector<L1>& l1_products);
void readSceneData(const netCDF::NcFile& nc, const size_t alt_beg, const size_t n_images, std::vector<L1>& l1_products);

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











} // namespace tango
