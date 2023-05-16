// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef INVERSION_H
#define INVERSION_H

#include "header.h"
#include "logger.h"

// A class for a iterative non-linear inversion.
// TODO: Expand this with more features like in Sicor or RemoTeC.
// This class is to be inherited from. Each inversion system has to define
// her own forward model and the desired additional parameters.
// So, they inherit from this system.
class Inversion : public Logger {

    protected:
    size_t nstate = 0; // Number of state parameters.
    size_t nmeas = 0; // Number of elements in the measurement vector.
    vector<double> own_apriori; // A-priori state when transferring with copy.
    double *apriori; // Pointer to apriori state.
    vector<double> own_meas; // Measurement vector when transferring with copy.
    double *meas; // Measurement vector.
    vector<double> s_y; // Error covariance matrix (if used).
    vector<bool> own_mask; // Pixel mask, created with copy or by default.
    vector<bool> *mask; // Pixel mask, eventually a pointer to its own mask, or to an external one.
    vector<double> own_tol; // Tolerance for convergence in state vector when transfoerring with copy.
    double *tol; // Tolerance for convergence in state vector.
    size_t maxiter; // Option: Maximum number of iterations to reach convergence.
    double steptolerance; // Maximum increment factor of chi square to accept a step.
    double initial_step_reducer; // Starting value of the Levenberg-Marquardt step control parameter.
    double reducerfactor_fail; // Multiplication of step reducer for a failed step.
    double reducerfactor_success; // Multiplication of step reducer for a successful step.
    double reducer_limit = NC_FILL_DOUBLE; // Maximum step control reducer required for convergence.

    public:
    Inversion(
        Logger *creator, // Object that creates this instance.
        size_t a_nstate, // Number of state parameters.
        size_t a_nmeas // Number of measurement.
    );
    virtual ~Inversion() {} // Do-nothing virtual destructor.

    // Simple setters.
    void setMaxiter(size_t a_maxiter) {maxiter = a_maxiter;}
    void setSteptolerance(double a_steptolerance) {steptolerance = a_steptolerance;}
    void setInitialStepReducer(double a_initial_step_reducer) {initial_step_reducer = a_initial_step_reducer;}
    void setReducerfactorFail(double a_reducerfactor_fail) {reducerfactor_fail = a_reducerfactor_fail;}
    void setReducerfactorSuccess(double a_reducerfactor_success) {reducerfactor_success = a_reducerfactor_success;}
    void setReducerLimit(double a_reducer_limit) {reducer_limit = a_reducer_limit;}

    // Less simple setters (memcpy).
    void setApriori(
        double *x,
        bool copy=false
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

    void setTol(
        double *t,
        bool copy=false
    );

    // Virtual forward model, to be overwritten by any inheriting class.
    protected:
    virtual int fwd(
        double *state, // State vector.
        double *signal, // Modelled signal (outout).
        double *jacobian // Modelled jacobian (outout).
    ) = 0;

    // The main function. Perform the inverison.
    public:
    int inversion_calculate(
        double *res, // Retrieved state (output).
        bool debug=false // Logs debugging information.
    );

    // Execution of linear inverison.
    virtual int execute_linear_invert(
        double *jac,
        double *residual,
        vector<bool> *excluded,
        double *state_des
    );

    // Simple getter.
    int getNstate() {return nstate;}

};

#endif
