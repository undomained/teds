// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef NONLINCAL_H
#define NONLINCAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class Bspline;
class CKD;

class Settings_nonlincal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_nonlincal(
        Logger *creator
    );
    ~Settings_nonlincal();
    // Specific settings.
    size_t nfrac = 1; // Number of fraction to break up L1A images.
    double lmin = NC_FILL_DOUBLE; // Lower border (in measured signal) of linear domain.
    double lmax = NC_FILL_DOUBLE; // Upper border (in measured signal) of linear domain.
    vector<double> exposure_time; // Exposure times for which the non-linearity is evaluated.
    double exptime_tol = 0.0; // Relative tolerance for exposure time grouping.
    bool scale_median = true; // Flag to scale the spline knots with the median of the saturated image. Set false to use individual-pixel saturation levels.
    vector<double> knots; // B-spline knots in signal domain.
    int order = 0; // B-spline order.
    // Pixel mask criteria.
    double mask_lin_chi2_max = NC_FILL_DOUBLE; // Maximum residual chi squared on the linear fit.
    double mask_chi2_max = NC_FILL_DOUBLE; // Maximum residual chi squared on non-linearity correction.

    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

// Detailed output of one series.
struct DetailedOutputNonlin { // {{{
    size_t nl1a; // Number of measurements in the exposure time series.
    vector<double> signal_meas; // Measured signal, for all pixels.
    vector<double> signal_corr; // Desired corrected signal, for all pixels.
}; // }}}

enum phase_t {
    phase_median = 0,
    phase_lin = 1,
    phase_nonlin = 2
};

class Nonlincal : public Processor { // {{{

    // A public constructor.
    public:
    Nonlincal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Nonlincal();

    // Overwritten virtual functions.
    private:
    unique_ptr<Settings_nonlincal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(size_t ibatch, const Calibration_options& opt) override;
    int process_finalize() override;

    int process_median(
        size_t iexptime
    );
    int process_lin(
        size_t ipix_start,
        size_t ipix_end
    );
    int process_nonlin(
        size_t ipix_start,
        size_t ipix_end,
        size_t iexptime
    );

    private:
    // NetCDF contents that is written, but not read.
    vector<double> nonlin_lin_slope; // Slope fit during linear fit.
    vector<double> nonlin_lin_chi2; // Residual chi squared on the linear fit.
    vector<double> nonlin_chi2; // Residual chi squared on non-linearity correction.

    // Batch administration.
    vector<phase_t> phase_batches; // Phase (0=get median, 1=find linear relationship, 2=fit non-linearity)
    vector<size_t> iexptime_batches; // Exposure time index that belongs to the batch.
    vector<size_t> ifrac_batches; // Fraction indices of batches.

    // Abscissa (exposure time times illumination level, for all L1A instances).
    vector<double> abscissa;

    // Fitting spline.
    unique_ptr<Bspline> b;

    // Optional output (level 1).
    vector<DetailedOutputNonlin> detailed_output_series;
    int write_detailed_output(
        NetCDF_object *nc,
        NcGroup &grp
    ) override;

}; // }}}

#endif
