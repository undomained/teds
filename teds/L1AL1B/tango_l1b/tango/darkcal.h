// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef DARKCAL_H
#define DARKCAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;

class Settings_darkcal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_darkcal(
        Logger *creator
    );
    ~Settings_darkcal();
    // Specific settings.
    int order = 0; // B-spline order for temperature fit (one higher than polynomial order).
    double nominal_temperature = NC_FILL_DOUBLE; // Standard temperature.
    double saturation_fraction_start = NC_FILL_DOUBLE; // Fraction of maximum signal to include in initial fit.
    double chi2_increment_threshold = NC_FILL_DOUBLE; // Threshold of chi square increment where to stop
    // Pixel mask criteria.
    double mask_chi2_max = NC_FILL_DOUBLE; // Maximum fit residual chi squared for dark fit.
    double mask_offset_min = NC_FILL_DOUBLE; // Minimum fixed dark signal (independent of integration time).
    double mask_offset_max = NC_FILL_DOUBLE; // Maximum fixed dark signal (independent of integration time).
    double mask_current_min = NC_FILL_DOUBLE; // Minimum dark signal added per second of integration time.
    double mask_current_max = NC_FILL_DOUBLE; // Maximum dark signal added per second of integration time.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;
}; // }}}

class Darkcal : public Processor { // {{{

    // A public constructor.
    public:
    Darkcal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Darkcal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_darkcal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(size_t ibatch, const Calibration_options& opt) override;

    private:
    // Diagnostic CKD.
    vector<double> dark_chi2; // Fit residual chi squared for dark fit.

    // Detailed output.
    size_t det_nl1a; // Copy of nl1a in batch or nl1a_total, but resillient against possible structure change in the algorithms.
    vector<double> medians; // Also used as working variable, but set as member because of detailed output option.
    vector<int> det_inclusion; // 1: In initial set. >1: Iteratively included and survived filter (in order 2,n). -1: Triggered filter, thus excluded. 0: Higher than measurement that triggered filter, so not even considered.
    vector<double> det_medianfit_chi2;
    int write_detailed_output(
        NetCDF_object *nc,
        NcGroup &grp
    ) override;

}; // }}}

#endif
