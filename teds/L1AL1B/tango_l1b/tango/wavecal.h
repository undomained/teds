// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef WAVECAL_H
#define WAVECAL_H

#include "header.h"
#include "settings_proc.h"
#include "processor.h"

// Forward declaration.
class CKD;

class Settings_wavecal : public Settings_proc { // {{{

    public:
    // Constructor.
    Settings_wavecal(
        Logger *creator
    );
    ~Settings_wavecal();
    // Specific settings.
    int peak_insignificant_island = 0; // Setting for peak detection: Size of peaks or non-peaks to be discarded.
    int peak_maxreso = 0; // Highest denominator in setting threshold to detect peaks.
    size_t peak_backgroundorder = 1; // Degrees of freedom of background polynomial for peak detection.
    bool peak_gaussfit = true; // Flag for performing gauss fit.
    double peak_domain_range = NC_FILL_DOUBLE; // Range in a-priori FWHM that a Gauss contains in split splitgaussfit.
    size_t peak_maxiter = 0; // Number of iterations for fitting Gauss peak shape.
    double peak_magtol = NC_FILL_DOUBLE; // Tolerance of peak magnitude (in units of a-priori magnitude).
    double peak_postol = NC_FILL_DOUBLE; // Tolerance of peak positions (in indices).
    double peak_widthtol = NC_FILL_DOUBLE; // Tolerance of peak FWHM (in indices).
    double peak_backgroundtol = NC_FILL_DOUBLE; // Tolerance of background B-spline coefficients (in median signals).
    double peak_steptolerance = NC_FILL_DOUBLE; // Inversion step acceptance criterion for fitting Gauss peak shape.
    double peak_initial_step_reducer = NC_FILL_DOUBLE; // Starting value of the Levenberg-Marquardt step control parameter.
    double peak_reducerfactor_fail = NC_FILL_DOUBLE; // Multiplication of step reducer for a failed step.
    double peak_reducerfactor_success = NC_FILL_DOUBLE; // Multiplication of step reducer for a successful step.
    double peak_reducer_limit = NC_FILL_DOUBLE; // Maximum step control reducer for convergence.
    int refspec_insignificant_island = 0; // In reference spectrum: Setting for peak detection: Size of peaks or non-peaks to be discarded.
    int refspec_maxreso = 0; // In reference spectrum: Highest denominator in setting threshold to detect peaks.
    size_t refspec_backgroundorder = 1; // In reference spectrum: Degrees of freedom of background polynomial for peak detection.
    bool refspec_gaussfit = true; // In reference spectrum: Flag for performing gauss fit.
    double refspec_domain_range = NC_FILL_DOUBLE; // In reference spectrum: Range in a-priori FWHM that a Gauss contains in split splitgaussfit.
    size_t refspec_maxiter = 0; // In reference spectrum: Number of iterations for fitting Gauss peak shape.
    double refspec_magtol = NC_FILL_DOUBLE; // In reference spectrum: Tolerance of peak magnitude (in units of a-priori magnitude).
    double refspec_postol = NC_FILL_DOUBLE; // In reference spectrum: Tolerance of peak positions (in indices).
    double refspec_widthtol = NC_FILL_DOUBLE; // In reference spectrum: Tolerance of peak FWHM (in indices).
    double refspec_backgroundtol = NC_FILL_DOUBLE; // In reference spectrum: Tolerance of background B-spline coefficients (in median signals).
    double refspec_steptolerance = NC_FILL_DOUBLE; // In reference spectrum: Inversion step acceptance criterion for fitting Gauss peak shape.
    double refspec_initial_step_reducer = NC_FILL_DOUBLE; // In reference spectrum: Starting value of the Levenberg-Marquardt step control parameter.
    double refspec_reducerfactor_fail = NC_FILL_DOUBLE; // In reference spectrum: Multiplication of step reducer for a failed step.
    double refspec_reducerfactor_success = NC_FILL_DOUBLE; // In reference spectrum: Multiplication of step reducer for a successful step.
    double refspec_reducer_limit = NC_FILL_DOUBLE; // In reference spectrum: Maximum step control reducer for convergence.
    int order = 0; // Degrees of freedom in fit of wavelength against pixel index.
    size_t npeak = 0; // Number of peak wavelengths.

    // Overwritten virtual function.
    protected:
    int init_step(
        stringstream &stream, // A string stream to use (just initialize one).
        string &key, // Name of the setting.
        string &value, // Where the value will be stored.
        bool &recognized // Return flag whether the setting is successfully recognized.
    ) override;

}; // }}}

class Wavecal : public Processor { // {{{

    // A public constructor.
    public:
    Wavecal(
        Logger *creator,
        CKD *ckd_arg
    );
    ~Wavecal();

    // Overwritten virtual functions.
    protected:
    unique_ptr<Settings_wavecal> set; // To ensure that everyone knows that set in this instance is of this derived type.
    int process_init() override;
    int process_batch(size_t ibatch, const Calibration_options& opt) override;

    private:
    // Batch-to-viewport coupling. Non-trivial if viewports are skipped.
    vector<size_t> ivp_batch;

    // Working objects.
    unique_ptr<Splitgaussfit> spl; // Peak fitter for measured spectra.
    unique_ptr<Bspline> b; // For any fit of wavelength as function of spectral index.
    vector<double> evalmat; // Matrix to evaluate the wavelength as function of spectral index.

}; // }}}

#endif
