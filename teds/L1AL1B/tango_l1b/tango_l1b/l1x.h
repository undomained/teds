// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "header.h"
#include "logger.h"

namespace tango {

// Forward declaration.
class NetCDF_object;
class CKD;
class L1A;

// L1X file administration.
class L1X: public Logger { // {{{

    public:
    // Constructor.
    L1X(
        Logger *creator
    );
    ~L1X(); // Destructor.

    // Initialization.
    int init(
        l1x_t il1x, // L1X level identifier.
        size_t nframe, // Number of frames (assessed by L1A manager).
        string &filename, // L1X file name.
        L1A_file_metadata *file_meta, // L1A file metadata.
        CKD *ckd, // Calibration key data (for detector and spectrum dimensions).
        size_t largest_spectrum_size // NetCDF dimension size for extracted spectra (only for FOV and later L1X steps).
    );

    int write_metadata(
        L1A *l1a // L1A instance, which contains all the metadata.
    );

    // Writes observation in L1X, contents depend on step.
    int write(
        L1A *l1a, // L1A instance that contains the information.
        size_t ifov=NC_FILL_UINT64 // Field-of-view coordinate (only for spectra, L1X_FOV and L1X_RAD).
    );

    // Member variables.
    private:
    l1x_t step; // Save step number, so it knows what it can do and what not.
    bool execute = false; // Flag for this L1X step to be executed or not. Will be set to true if it reaches its init function.
    unique_ptr<NetCDF_object> nc; // NetCDF object of the L1X, can be shared for e.g. flight-style L1X.
    // NetCDF variables. {{{
    NcVar var_image_time;
    NcVar var_image_CCSDS_sec;
    NcVar var_image_CCSDS_usec;
    NcVar var_binning_table;
    NcVar var_digital_offset;
    NcVar var_nr_coadditions;
    NcVar var_exposure_time;
    NcVar var_hk_tlm_time;
    NcVar var_temp_detector;
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
    // }}}

}; // }}}

} // namespace tango
