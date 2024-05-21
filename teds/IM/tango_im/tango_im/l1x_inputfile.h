#pragma once

#include "header.h"
#include "logger.h"

namespace tango {

// Forward declaration.
class Netcdf_object;
class Settings_main;
class CKD;

// Class for one image and everything that is being done with it.
class L1X_inputfile : public Logger { // {{{

    public:
    L1X_inputfile(
        Logger *creator
    );
    ~L1X_inputfile();

    int init(
        Settings_main *set,
        CKD *ckd
    );

    shared_ptr<NetCDF_object> nc;

    // General information.
    size_t nframe;
    vector<size_t> subset; // Eventual subset of L1X file.
    l1x_t il1x_start;

    // Truth L1X dimensions that are independent of CKD.
    size_t dim_spat_truth; // Across-track sampling.
    size_t dim_spec_truth; // Spectral sampling.
    // Frame-independent science data.
    vector<double> wavelength; // Wavelength grid, it does not change per observations. Dimensions: dim_vp,dim_spec_truth (dim_vp comes from CKD).

    bool navi_exists; // Flag for if navigation data exists.

    // File GSE data.
    bool gse_exists; // Flag whether there is GSE data.
    double illumination_level = NC_FILL_DOUBLE; // Scalar illumination level for detector measurements.
    uint8_t viewportinteger = NC_FILL_UBYTE; // Number representing illumated viewports.
    size_t dim_refspec; // Size of the reference spectrum.
    vector<double> refspec_wavelength; // Wavelengths of reference spectrum.
    vector<double> refspec_radiance; // Radiance values or reference spectrum.
    double refspec_dolp = NC_FILL_DOUBLE; // Constant over the entire spectrum.
    double refspec_aolp = NC_FILL_DOUBLE; // Constant over the entire spectrum.
    double act_angle = NC_FILL_DOUBLE; // Across-track rotation angle of rotation stage.
    double alt_angle = NC_FILL_DOUBLE; // Along-track rotation angle of rotation stage.

}; // }}}

} // namespace tango
