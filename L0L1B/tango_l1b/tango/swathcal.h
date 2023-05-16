// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef SWATHCAL_H
#define SWATHCAL_H

#include "header.h"
#include "vector.h" // Because of member of type Quaternion.
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;
class Theodolite;

class Settings_swathcal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_swathcal(
        Logger *creator
    );
    ~Settings_swathcal();
    // Specific settings.
    string theodolite_reference_filename; // NetCDF file name of the theodolite reference file (a lot of measurement on the wrong day).
    string theodolite_spotcheck_filename; // NetCDF file name of the theodolite spot-check file (one measurement on correct day).
    vector<double> viewport_alt_angles; // Expected ALT angles of each viewport (also used to put the viewports in the correct order).
    double theodolite_alt_angle_tol = -1.0; // Maximum allowed different between ALT angle of theodolite measurement and expected ALT angle of a viewport.
    double actscan_sampling = NC_FILL_DOUBLE; // Sampling for initial scan over swath to find the correct FOVcal ACT angle (default = just use swathcal stage ACT angle).
    double act_linesearch_tol = -1.0; // Convergence tolerance in line search over ACT angles.
    double act_angle_tol = -1.0; // Difference between across-track angles that are considered to be the same.
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
    vector<size_t> order; // Degrees of freedom for polynomial fit of ALT against ACT angles per viewport, each with a minimum number of successful associations. Start with highest.
    vector<size_t> min_fit; // Minimum required associations to perform the fit with the respective order. Start with highest.
    vector<double> rotmat; // Rotation matrix from spacecraft to instrument coordinates. Nine elements. Source direction is quick. Destination dimension is slow.
    double rot_warningthreshold = 0.0; // Allowed error in check if rotmat is a proper rotation matrix.
    size_t theo_altmodel_maxiter = 0; // Maximum number of iterations in non-linear inversion for ALT-model in the theodolite.
    double theo_altmodel_tol = NC_FILL_DOUBLE; // Convergence tolerance on vector elements in non-linear inversion for ALT-model in the theodolite.
    double theo_altmodel_steptolerance = NC_FILL_DOUBLE; // Inversion step acceptance criterion for non-linear inversion for ALT-model in the theodolite.
    double theo_altmodel_initial_step_reducer = NC_FILL_DOUBLE; // Starting value of the Levenberg-Marquardt step control parameter for non-linear inversion for ALT-model in the theodolite.
    double theo_altmodel_reducerfactor_fail = NC_FILL_DOUBLE; // Multiplication of step reducer for a failed step for non-linear inversion for ALT-model in the theodolite.
    double theo_altmodel_reducerfactor_success = NC_FILL_DOUBLE; // Multiplication of step reducer for a successful step for non-linear inversion for ALT-model in the theodolite.
    double theo_altmodel_reducer_limit = NC_FILL_DOUBLE; // Maximum step control reducer for convergence for non-linear inversion for ALT-model in the theodolite.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

class Swathcal : public Processor { // {{{

    // A public constructor.
    public:
    Swathcal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Swathcal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_swathcal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(
        size_t ibatch
    ) override;
    int process_finalize() override;

    private:
    // Administration.
    vector<size_t> ibatch_start;
    vector<size_t> nbatch_vp;
    vector<size_t> nfail_vp;
    // Raw fit results.
    vector<double> act_fovcal;
    vector<double> act_swathcal;
    vector<double> alt_swathcal;
    vector<bool> mask;
    vector<double> mag;
    // Intermediate output to write as detailed output.
    vector<double> l1b_act_angles;
    vector<double> l1b_alt_angles;
    // B-spline order attempts.
    size_t nattempt;
    size_t min_fit_poorest; // For easy access to determine certain failure.
    // Rotation quaternion.
    Quaternion q_inst; // Instrument rotation quaternion.
    unique_ptr<Theodolite> theo;

    // Detailed output.
    int write_detailed_output(
        NetCDF_object *nc,
        NcGroup &grp
    ) override;

    // Private function for spectrum extraction, an often-repeated action
    // during process_batch. This extraction does not work for binned images.
    int act_extract(
        size_t ivp, // Viewport index.
        double *image, // Detector image.
        vector<bool> &pixelmask, // Detector pixel mask.
        double act, // ACT-angle in FOVcal perspective.
        double &strength // Output spectrum strength (median of two spectra).
    );

}; // }}}

#endif
