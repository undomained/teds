// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "header.h"
#include "matrix.h" // Required because of member variable of type squareoption_t.

namespace tango {

class Lininv {
    public:
    // Input.
    size_t nstate;
    size_t nmeas;
    double *jac;
    squareoption_t opt_s_y;
    double *s_y;
    vector<bool> excluded;
    double *residual;
    // Ouptut.
    bool opt_res = true;
    bool opt_sx = false;
    bool opt_gain = false;
    vector<double> res;
    vector<double> sx;
    vector<double> gain;

    int execute();
    ~Lininv();

};

// Out-of-class inversion routines, using the class internally.
int linear_invert(
    size_t nstate, // Number of state parameters.
    size_t nmeas, // Number of measurements, should be at least as high as the number of state parameters.
    double *jac, // Jacobian matrix.
    squareoption_t opt_s_y, // Type of error covariance matrix.
    double *s_y, // Error covariance matrix.
    vector<bool> *mask, // Pixel mask.
    double *residual, // Measurement residual.
    double *res, // Retrieved state.
    double *sx=NULL // Error covariance on retrieved state.
);

int linear_invert_gain(
    size_t nstate, // Number of state parameters.
    size_t nmeas, // Number of measurements, should be at least as high as the number of state parameters.
    double *jac, // Jacobian matrix.
    squareoption_t opt_s_y, // Type of error covariance matrix.
    double *s_y, // Error covariance matrix.
    vector<bool> *mask, // Pixel mask.
    double *gain // Gain matrix.
);

// Linear inversion with hard constraints. Works the same as a regular linear
// inversion, but as a boundary condition, an ill-posed equation jac*res = residual
// is solved and the regular equation performs a least square fit using the
// remaining degrees of freedom in res.
int constrained_linear_invert(
    size_t nstate, // Number of state paramters.
    size_t nmeas_soft, // Number of measurements for least-square fit.
    size_t nmeas_hard, // Number of hard constraints.
    double *jac_soft, // Jacobian matrix for least-square fit.
    double *jac_hard, // Jacobian matrix for hard constraints.
    squareoption_t opt_s_y_soft, // Definition option soft error covariance matrix.
    double *s_y_soft, // Soft error covariance matrix.
    vector<bool> *mask_soft, // Soft pixel mask.
    double *residual_soft, // Residual measurements for least-square fit.
    double *residual_hard, // Residual measurements for hard constraints.
    double *res // The result vector.
);

// Linear inversion with outlier protection. Outlier measurements are discarded based on a threshold.
int linear_invert_outlierprotect(
    size_t nstate, // Number of state parameters.
    size_t nmeas, // Number of measurements, should be at least as high as the number of state parameters.
    double *jac, // Jacobian matrix.
    squareoption_t opt_s_y, // Type of error covariance matrix.
    double *s_y, // Error covariance matrix.
    vector<bool> *mask, // Output array of which measurement points were excluded.
    double *residual, // Measurement residual.
    double outlier_cutoff, // Deviation criterion on which measurements are discarded (in SD).
    double *res // Retreived state.
);

} // namespace tango
