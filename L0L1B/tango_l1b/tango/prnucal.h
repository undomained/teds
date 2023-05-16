// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef PRNUCAL_H
#define PRNUCAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;

class Settings_prnucal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_prnucal(
        Logger *creator
    );
    ~Settings_prnucal();
    // Specific settings.
    size_t order_spat = 0; // Low-order polynomial DOFs for spatial dimension (one higher than order).
    size_t order_spec = 0; // Low-order polynomial DOFs for spectral dimension (one higher than order).
    double outlier_cutoff = NC_FILL_DOUBLE; // Outlier-protection threshold for low-order polynomial fit.
    bool iterative = false; // Flag for letting the pixel mask and the PRNU iteratively converge.
    // Pixel mask criteria.
    double mask_prnu_min = NC_FILL_DOUBLE; // Minimum pixel response non-uniformity.
    double mask_prnu_max = NC_FILL_DOUBLE; // Maximum pixel response non-uniformity.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

class Prnucal : public Processor { // {{{

    // A public constructor.
    public:
    Prnucal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Prnucal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_prnucal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(
        size_t ibatch
    ) override;

    private:
    // Working variables that are saved if detailed output is selected.
    vector<double> photons; // Light exposure (illumination strength multiplied by exposure time).
    vector<double> sensitivity; // Pixel sensitivity, still including low-order terms.
    vector<double> sensitivity_fit; // Low-order part of the pixel sensitivity.
    // Dimension size of number of measurement in s less clumsy way.
    size_t det_nl1a;
    int write_detailed_output(
        NetCDF_object *nc,
        NcGroup &grp
    ) override;

}; // }}}

#endif
