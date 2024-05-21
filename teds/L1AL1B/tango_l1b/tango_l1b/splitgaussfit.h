// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef SPLITGAUSSFIT_H
#define SPLITGAUSSFIT_H

#include "header.h"
#include "logger.h"

// This is an inversion routine that attempts to fit a set of Gaussians.
// This is used to recognize a number of peaks in a line.
class Splitgaussfit : public Logger { // {{{

    // Constructor.
    public:
    Splitgaussfit(
        Logger *creator, // Object that creates this one.
        size_t a_ngauss, // Number of peaks (gausses).
        size_t a_npoly, // Number of polynomial parameters for the background.
        size_t a_nmeas // Number of elements in the measurement vector.
    );
    ~Splitgaussfit(); // Destructor.
    private:
    size_t ngauss; // Total number of peaks.
    size_t npoly; // Number of polynomial terms, one higher than the order.
    size_t nmeas;
    vector<double> own_x;
    double *x;
    bool new_x;
    vector<double> own_meas;
    double *meas;
    vector<double> s_y;
    vector<bool> own_mask; // Pixel mask, created with copy or by default.
    vector<bool> *mask; // Pixel mask, eventually a pointer to its own mask, or to an external one.
    vector<double> mat;

    // Execution of the inversion.
    bool gaussfit = true;

    // Inversion properties.
    size_t inv_maxiter; // Maximum number of iterations to reach convergence.
    double inv_steptolerance; // Tolerance for a step to be recognized as successful.
    double inv_initial_step_reducer = NC_FILL_DOUBLE; // Starting value of the Levenberg-Marquardt step control parameter.
    double inv_reducerfactor_fail = NC_FILL_DOUBLE; // Multiplication of step reducer for a failed step.
    double inv_reducerfactor_success = NC_FILL_DOUBLE; // Multiplication of step reducer for a successful step.
    double inv_reducer_limit = NC_FILL_DOUBLE; // Maximum step control reducer for convergence.

    // Inversion convergence tolerances.
    double inv_magtol; // Tolerance in peak magnitude in units of a-priori magnitude.
    double inv_postol; // Tolerance in peak position in X-coordinate units.
    double inv_widthtol; // Tolerance in peak FWHM in X-coordinate units.
    double inv_backgroundtol; // Tolerance in background in line medians.

    // Setter.
    public:
    void setX(
        double *a_x, // Independent parameter of the line.
        bool copy=false // Flag for copying the input data. Leave false to pass pointer.
    );

    void setMeas(
        double *y,
        bool copy=false
    );

    void setNoise(
        double *y_noise,
        bool s=false // Flag for providing S-elements (no need to square anymore).
    );
    void removeNoise();
    void setMask(
        vector<bool> *a_mask,
        bool copy=false
    );
    void removeMask();

    // Set flag for executing iterative fit of a Gaussian.
    void setGaussfit(bool a_gaussfit) {gaussfit=a_gaussfit;}

    // Simple setters.
    void setMaxiter(size_t a_maxiter) {inv_maxiter=a_maxiter;}
    void setSteptolerance(double a_steptolerance) {inv_steptolerance=a_steptolerance;}
    void setInitialStepReducer(double a_initial_step_reducer) {inv_initial_step_reducer = a_initial_step_reducer;}
    void setReducerfactorFail(double a_reducerfactor_fail) {inv_reducerfactor_fail = a_reducerfactor_fail;}
    void setReducerfactorSuccess(double a_reducerfactor_success) {inv_reducerfactor_success = a_reducerfactor_success;}
    void setReducerLimit(double a_reducer_limit) {inv_reducer_limit = a_reducer_limit;}

    // Tolerance setters.
    void setMagtol(double a_magtol) {inv_magtol=a_magtol;}
    void setPostol(double a_postol) {inv_postol=a_postol;}
    void setWidthtol(double a_widthtol) {inv_widthtol=a_widthtol;}
    void setBackgroundtol(double a_backgroundtol) {inv_backgroundtol=a_backgroundtol;}

    int solve(
        int insignificant_island, // Islands of this many pixels are not recognized.
        double domain_range, // Size of a gauss tail to form domains.
        int maxreso, // Highest denominator in setting thresholds.
        double *res, // Retrieval result.
        bool debug=false // Logs debug information.
    );

}; // }}}

#endif
