// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "constants.h"
#include "spectrum.h"
#include "setting.h"
#include <netcdf>
#include <memory>
#include "ckd.h"

namespace tango {

class L1;
class CKD;

class L1Measurement {
private:
    std::string l1_filename;
    std::string l1_level;
    std::vector<L1> l1_measurement;
    std::vector<double> image_time;
    std::vector<uint8_t> binning_table;
    std::vector<uint16_t> nr_coadditions;
    std::vector<double> exposure_time;
    size_t alt_beg {};
    size_t alt_end {};
    size_t n_images {};

    void setLevel(const netCDF::NcFile& nc);
    void readMetaData(const netCDF::NcFile& nc, const std::string& config);
    void readScienceData(const netCDF::NcFile& nc);
    void readSceneData(const netCDF::NcFile& nc);
    void writeGlobalAttributes(netCDF::NcFile& nc, const std::string& level, const std::string& config, const int argc, const char* const argv[]);
    void writeMetaData(netCDF::NcFile& nc);
    void writeScienceData(netCDF::NcFile& nc);
    void writeObservationData(netCDF::NcFile& nc);
    void writeGeolocationData(netCDF::NcFile& nc);

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

public:
    /// Constructor.
    L1Measurement(
            const std::string filename,
            const int image_start,
            const int image_end,
            const std::string config);

    /// Destructor.
    ~L1Measurement() = default;

    void read(const std::string& filename,
        const std::string& config);

    void write(const std::string& filename,
        const std::string& level,
        const std::string& config, 
        const int argc, const char* const argv[]);

    // Some functions defined public for accessing the measurement points/images/L1 objects
    // Limited acces. Do not want the outside world to be able to delete, resize, remove L1 objects
    // add L1 objects etc.
    L1& operator[](const int img);
    int size();
    L1& front();
    L1& back();

    // Copy geolocation data from the geometry file directly to the L1B
    // product. This is a placeholder function until geolocation is
    // properly implemented.
    void copyGeometry(const std::string& config);
    //void copyGeometry(const std::string& l1a_filename,
    //                  const std::string& geo_filename);
    //                  int i_alt_start,
    //                  std::vector<L1>& l1_products) -> void;

    void extractSpectraWavelength(L1& l1, const CKD& ckd);
};


} // namespace tango
