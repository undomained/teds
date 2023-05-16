// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef FOVCAL_H
#define FOVCAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class Bspline;
class Splitgaussfit;
class CKD;

class Settings_fovcal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_fovcal(
        Logger *creator
    );
    ~Settings_fovcal();
    // Specific settings.
    int order = 0; // Polynomial DOFs to fit the curves on the detector.
    int insignificant_island = 0; // Setting for peak detection: Size of peaks or non-peaks to be discarded.
    int maxreso = 0; // Highest denominator in setting threshold to detect peaks.
    size_t backgroundorder = 1; // Degrees of freedom of background polynomial for peak detection.
    bool gaussfit = true; // Flag for performing gauss fit.
    double domain_range = NC_FILL_DOUBLE; // Range in a-priori FWHM that a Gauss contains in split splitgaussfit.
    size_t maxiter = 0; // Number of iterations for fitting Gauss peak shape.
    double magtol = NC_FILL_DOUBLE; // Tolerance of peak magnitude (in units of a-priori magnitude).
    double postol = NC_FILL_DOUBLE; // Tolerance of peak positions (in indices).
    double widthtol = NC_FILL_DOUBLE; // Tolerance of peak FWHM (in indices).
    double backgroundtol = NC_FILL_DOUBLE; // Tolerance of background B-spline coefficients (in median signals).
    double steptolerance = NC_FILL_DOUBLE; // Inversion step acceptance criterion for fitting Gauss peak shape.
    double initial_step_reducer = NC_FILL_DOUBLE; // Starting value of the Levenberg-Marquardt step control parameter.
    double reducerfactor_fail = NC_FILL_DOUBLE; // Multiplication of step reducer for a failed step.
    double reducerfactor_success = NC_FILL_DOUBLE; // Multiplication of step reducer for a successful step.
    double reducer_limit = NC_FILL_DOUBLE; // Maximum step control reducer for convergence.
    size_t ispec_start = 0; // First spectral pixel index to use for fitting polynomial.
    size_t ispec_end = NC_FILL_UINT64; // First spectral pixel index not to use for fitting polynomial.
    double outlier_cutoff = NC_FILL_DOUBLE; // As long as the difference between polynomial and data is this many times sigma, we remove the worst point.
    double hitthreshold = NC_FILL_DOUBLE; // Threshold for spectrum median intensity to be considered detected (fraction of highest).
    vector<double> act_sampling; // ACT angle sampling per viewport.
    vector<double> force_actangle_min; // Force minimum ACT angle as L1B FOV. Warning if not hit. This is one element per viewport.
    vector<double> force_actangle_max; // Force maximum ACT angle as L1B FOV. Warning if not hit. This is one element per viewport.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

// Detailed output per viewport.
struct DetailedOutputFov {
    size_t nl1a; // Number of measurements.
    // All arrays are stored in the sorted way, thus sorted in ascending ACT angles.
    vector<double> act_angles; // ACT angles of the files.
    vector<double> ispat_raw; // Raw fit results. Dimension (L1A,pol,detector_spec).
    vector<double> ispat_smooth; // Smoothed fit results, still per L1A file. Dimension (L1A,pol,detector_spec).
    vector<double> spectrummedians; // Median spectra for visibility assessment. Dimension (L1A).
    vector<uint8_t> visible; // Visibility flag (L1A).
};

class Fovcal : public Processor { // {{{

    // A public constructor.
    public:
    Fovcal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Fovcal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_fovcal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(
        size_t ibatch
    ) override;
    int process_finalize() override;

    private:
    // Batch-to-viewport coupling. Non-trivial if viewports are skipped.
    vector<size_t> ivp_batch;
    // Saved content from process_init to process_batch.
    unique_ptr<Splitgaussfit> spl;
    unique_ptr<Bspline> b;
    vector<double> mat_full;

    // Optional output (level 1).
    vector<DetailedOutputFov> detailed_output_viewports;
    int write_detailed_output(
        NetCDF_object *nc,
        NcGroup &grp
    ) override;

}; // }}}

#endif
