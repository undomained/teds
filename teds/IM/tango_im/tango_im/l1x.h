#pragma once

#include "header.h"
#include "logger.h"

namespace tango {

// Forward declaration.
class NetCDF_object;
class Settings_noise;
class CKD;
class Frame;
class L1X_inputfile;

// L1X file administration.
class L1X: public Logger { // {{{

    public:
    // Constructor.
    L1X(
        Logger *creator
    );
    ~L1X(); // Destructor.

    // Only for L1A.
    int binning_table_id;
    size_t npix_binned;

    // Initialization.
    int init(
        l1x_t il1x, // L1X level identifier.
        size_t nframe, // Number of frames.
        string &filename, // L1X file name.
        size_t a_npix_binned, // Number of pixels in binned format. Equal to entire detector for non-L1A.
        int a_binning_table_id, // Binning table ID. Zero for non-L1A.
        CKD *ckd, // Calibration key data (for detector and spectrum dimensions).
        L1X_inputfile *l1x_inputfile // L1X input file structure.
    );

    int write_metadata(
        Frame *frm // Frame to write L1X from.
    );

    // Writes observation in L1X, contents depend on step.
    int write(
        Frame *frm, // Frame to write L1X from.
        Settings_noise *set, // Detector settings.
        CKD *ckd
    );

    // Member variables.
    private:
    size_t npix; // So that the CKD need not be dragged all over the place.
    l1x_t step; // Save step number, so it knows what it can do and what not.
    bool execute = false; // Flag for this L1X step to be executed or not. Will be set to true if it reaches its init function.
    unique_ptr<NetCDF_object> nc; // NetCDF object of the L1X, can be shared for e.g. flight-style L1X.
    // NetCDF variables. {{{
    NcVar var_image_time;
    NcVar var_binning_table;
    NcVar var_nr_coadditions;
    NcVar var_exposure_time;
    NcVar var_img; // Used for both binned and unbinned. TODO: This could be improved.
    NcVar var_img_mask; // Used for both binned and unbinned. TODO: This could be improved.
    NcVar var_img_noise; // Used for both binned and unbinned. TODO: This could be improved.
    NcVar var_img_with_current;
    NcVar var_spectra_sizes;
    NcVar var_spectra;
    NcVar var_spectra_noise;
    NcVar var_spectra_mask;
    NcVar var_spectra_target;
    NcVar var_spectra_target_noise;
    NcVar var_spectra_target_mask;
    NcVar var_att_time;
    NcVar var_att_qaut;
    NcVar var_orb_time;
    NcVar var_orb_pos;
    NcVar var_orb_vel;
    /*
     //The follows has been added for the simplified GM
    vector<double> nav_latitude;
    vector<double> nav_longitude;
    vector<double> nav_velocity_lat;
    vector<double> nav_velocity_lon;
    vector<double> nav_altitude;
    vector<double> nav_roll_ang;
    size_t nav_dim_navi;
    vector<double> nav_timestamps;
    */
    double nav_latitude;
    double nav_longitude;
    double nav_velocity_lat;
    double nav_velocity_lon;
    double nav_altitude;
    double nav_roll_ang;
    size_t nav_dim_navi;
    vector<double> nav_timestamps;
    //The end of the added for the simple GM
    
    // }}}

}; // }}}

} // namespace tango
