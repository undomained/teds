// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "matrix.h"
#include "lininv.h"

// Sorry, this interface is terrible, all arguments in a different form.
// Unfortunately, this is the only one that works.
int invert_s_y( // {{{
    size_t nmeas,
    squareoption_t opt_s_y,
    double *s_y,
    vector<bool> excluded,
    squareoption_t &opt_s_y_use,
    vector<double> &inv_s_y
)
{
    if (opt_s_y == OPT_DIAG) {
        inv_s_y.resize(nmeas);
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
            inv_s_y[imeas] = 1.0/s_y[imeas];
        }
        if (excluded.size() != 0) {
            // Kill off excluded points from inverse of s_y.
            for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                if (excluded[imeas]) inv_s_y[imeas] = 0.0;
            }
        }
    }
    if (opt_s_y == OPT_FULL) {
        inv_s_y.resize(nmeas*(nmeas+1)/2);
        memcpy(inv_s_y.data(),s_y,nmeas*(nmeas+1)/2*sizeof(double));
        if (excluded.size() != 0) {
            // Remove correlation terms with excluded points.
            double *row = inv_s_y.data(); // Currnetly inv_s_y is still a copy of s_y, which we can demolish.
            for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                row += imeas;
                if (excluded[imeas]) {
                    for (size_t icov=0 ; icov<imeas ; icov++) row[icov] = 0.0;
                }
                row[imeas] = 1.0;
            }
        }
        handle_nonlethal(chol_solve(nmeas,inv_s_y.data()));
        if (excluded.size() != 0) {
            // Set diagonal terms of excluded points to zero.
            double *row = inv_s_y.data();
            for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                row += imeas;
                if (excluded[imeas]) row[imeas] = 0.0;
            }
        }
    }
    if (opt_s_y == OPT_NONE && excluded.size() != 0) {
        opt_s_y_use = OPT_DIAG;
        inv_s_y.resize(nmeas);
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
            inv_s_y[imeas] = excluded[imeas]?0.0:1.0;
        }
    } else {
        opt_s_y_use = opt_s_y;
    }

    return 0;

} // }}}

Lininv::~Lininv() {}
int Lininv::execute( // {{{
)
{

    // Invert the error covariance matrix with all the options.
    squareoption_t opt_s_y_use; // Output variable.
    vector<double> inv_s_y; // Output variable.
    handle_nonlethal(invert_s_y(nmeas,opt_s_y,s_y,excluded,opt_s_y_use,inv_s_y));

    // If opt_s_y is OPT_NONE and no mask is given, the inv_s_y remains
    // uninitialized and inv_s_y.data() should become a null pointer.
    sx.resize(nstate*(nstate+1)/2);
    Matrix::matmul_tsk(nstate,nmeas,jac,opt_s_y_use,inv_s_y.data(),sx.data());
    handle_nonlethal(chol_solve(nstate,sx.data()));
    if (opt_res) {
        res.resize(nstate);
        Matrix::matmul_stsy(nstate,nmeas,jac,OPT_FULL,sx.data(),opt_s_y_use,inv_s_y.data(),residual,res.data());
    }
    if (opt_gain) {
        gain.resize(nstate*nmeas);
        Matrix::matmul_sts(nstate,nmeas,jac,OPT_FULL,sx.data(),opt_s_y_use,inv_s_y.data(),gain.data());
    }
    // Remove by-product sx if the user is not interested.
    if (!opt_sx) sx.resize(0);

    // Check for finite numbers.
    if (opt_res) {
        for (size_t istate=0 ; istate<nstate ; istate++) {
            if (!isfinite(res[istate])) return 1;
        }
    }
    if (opt_gain) {
        for (size_t iel=0 ; iel<nstate*nmeas ; iel++) {
            if (!isfinite(gain[iel])) return 1;
        }
    }
    if (opt_sx) {
        for (size_t iel=0 ; iel<nstate*(nstate+1)/2 ; iel++) {
            if (!isfinite(sx[iel])) return 1;
        }
    }
    return 0;
} // }}}

// Out-of-class linear inversion routines, using the class internally.
int linear_invert( // {{{
    size_t nstate, // Number of state parameters.
    size_t nmeas, // Number of measurements, should be at least as high as the number of state parameters.
    double *jac, // Jacobian matrix.
    squareoption_t opt_s_y, // Type of error covariance matrix.
    double *s_y, // Error covariance matrix.
    vector<bool> *mask, // Pixel mask.
    double *residual, // Measurement residual.
    double *res, // Retreived state.
    double *sx // Error covariance on retrieved state.
)
{
    Lininv l;
    l.nstate = nstate;
    l.nmeas = nmeas;
    l.jac = jac;
    l.opt_s_y = opt_s_y;
    l.s_y = s_y;
    l.residual = residual;
    if (mask != NULL) l.excluded = *mask; // Copies the vector, so that things can be changed wihtout changing the original. These changes happen during outlierportoect.
    if (sx != NULL) l.opt_sx = true;
    handle_nonlethal(l.execute());
    memcpy(res,l.res.data(),nstate*sizeof(double));
    if (sx != NULL) memcpy(sx,l.sx.data(),nstate*(nstate+1)/2*sizeof(double));
    return 0;
} // }}}
int linear_invert_gain( // {{{
    size_t nstate, // Number of state parameters.
    size_t nmeas, // Number of measurements, should be at least as high as the number of state parameters.
    double *jac, // Jacobian matrix.
    squareoption_t opt_s_y, // Type of error covariance matrix.
    double *s_y, // Error covariance matrix.
    vector<bool> *mask, // Pixel mask.
    double *gain // Gain matrix.
)
{
    Lininv l;
    l.nstate = nstate;
    l.nmeas = nmeas;
    l.jac = jac;
    l.opt_s_y = opt_s_y;
    l.s_y = s_y;
    if (mask != NULL) l.excluded = *mask; // Copies the vector, so that things can be changed wihtout changing the original. These changes happen during outlierportoect.
    l.opt_res = false;
    l.opt_gain = true;
    handle_nonlethal(l.execute());
    memcpy(gain,l.gain.data(),nstate*nmeas*sizeof(double));
    return 0;
} // }}}

// Linear inversion with hard constraints. Works the same as a regular linear
// inversion, but as a boundary condition, an ill-posed equation jac*res = residual
// is solved and the regular equation performs a least square fit using the
// remaining degrees of freedom in res.
int constrained_linear_invert( // {{{
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
)
{

    // Regularize with KTK of the hard.
    vector<double> reg_soft(nstate*(nstate+1)/2);
    Matrix::matmul_tsk(nstate,nmeas_hard,jac_hard,OPT_NONE,NULL,reg_soft.data());
    // Get an inverse of the soft S_y matrix.
    vector<double> mat_soft(nstate*(nstate+1)/2);
    squareoption_t opt_s_y_soft_use;
    vector<double> inv_s_y_soft;
    vector<bool> excluded_soft;
    if (mask_soft != NULL) excluded_soft = *mask_soft;
    handle_nonlethal(invert_s_y(nmeas_soft,opt_s_y_soft,s_y_soft,excluded_soft,opt_s_y_soft_use,inv_s_y_soft));
    Matrix::matmul_tsk(nstate,nmeas_soft,jac_soft,opt_s_y_soft_use,inv_s_y_soft.data(),mat_soft.data());
    for (size_t imat=0 ; imat<nstate*(nstate+1)/2 ; imat++) mat_soft[imat] += reg_soft[imat];
    handle_nonlethal(chol_solve(nstate,mat_soft.data()));
    vector<double> res_soft(nstate);
    Matrix::matmul_stsy(nstate,nmeas_soft,jac_soft,OPT_FULL,mat_soft.data(),opt_s_y_soft_use,inv_s_y_soft.data(),residual_soft,res_soft.data());
    vector<double> augmented_residual_hard(nmeas_hard);
    Matrix::matmul_rl_fold_slow(nstate,nmeas_hard,jac_hard,res_soft.data(),augmented_residual_hard.data());
    for (size_t imeas=0 ; imeas<nmeas_hard ; imeas++) augmented_residual_hard[imeas] = residual_hard[imeas] - augmented_residual_hard[imeas];
    // Now, we solve the hard (ill-posed) problem with this S_x.
    vector<double> mat_hard(nmeas_hard*(nmeas_hard+1)/2);
    Matrix::matmul_kst(nstate,nmeas_hard,jac_hard,OPT_FULL,mat_soft.data(),mat_hard.data());
    handle_nonlethal(chol_solve(nmeas_hard,mat_hard.data()));
    Matrix::matmul_stsy(nstate,nmeas_hard,jac_hard,OPT_FULL,mat_soft.data(),OPT_FULL,mat_hard.data(),augmented_residual_hard.data(),res);
    for (size_t istate=0 ; istate<nstate ; istate++) res[istate] += res_soft[istate];

    // Check for finite numbers.
    for (size_t istate=0 ; istate<nstate ; istate++) {
        if (!isfinite(res[istate])) return 1;
    }
    return 0;

} // }}}

int linear_invert_outlierprotect( // {{{
    size_t nstate, // Number of state parameters.
    size_t nmeas, // Number of measurements, should be at least as high as the number of state parameters.
    double *jac, // Jacobian matrix.
    squareoption_t opt_s_y, // Type of error covariance matrix.
    double *s_y, // Error covariance matrix.
    vector<bool> *mask, // Output array of which measurement points were excluded.
    double *residual, // Measurement residual.
    double outlier_cutoff, // Deviation criterion on which measurements are discarded (in SD).
    double *res // Retreived state.
)
{

    vector<bool> excluded_own;
    if (mask == NULL) {
        excluded_own.resize(nmeas,false);
        mask = &excluded_own;
    }
    vector<bool> &excluded = *mask; // To avoid clumsy constructions like (*mask)[imeas].
    // Thus, mask is the pointer, excluded is the array itself. Use mask to alter

    size_t nfit = nmeas;
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
        if (excluded[imeas]) nfit--;
    }

    bool ctd = true;
    while (ctd) {

        ctd = false;

        handle_nonlethal(linear_invert(nstate,nmeas,jac,opt_s_y,s_y,mask,residual,res));

        // There may be numeric problems if the number of living pixels approaches nstate.
        // In such a case, it is not the worst thing in the world if the fit is flagged
        // as unsuccessful. An outlier protect fit is generally done with a lot of measurements.

        double agg = 0.0;
        double aggsq = 0.0;
        vector<double> diff_squared(nmeas,0.0);
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
            if (excluded[imeas]) continue;
            double mod = 0.0;
            for (size_t istate=0 ; istate<nstate ; istate++) {
                mod += jac[istate*nmeas+imeas]*res[istate];
            }
            double sig = 1.0; // If opt_s_y == OPT_NONE
            if (opt_s_y == OPT_DIAG) sig = s_y[imeas];
            if (opt_s_y == OPT_FULL) sig = s_y[(imeas+1)*(imeas+2)/2-1];
            double diff = (residual[imeas] - mod) / sig;
            diff_squared[imeas] = pow(diff,2.0);
            agg += diff;
            aggsq += diff_squared[imeas];
        }

        // Evaluate the variance of this mismatch.
        double avg = agg/nfit;
        double avgsq = aggsq/nfit;
        double variance = avgsq - pow(avg,2.0);
        if (variance <= 0.0) break; // Zero sigma or numeric error if a real sigma is calculated with a square root.

        for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
            if (diff_squared[imeas] / variance > pow(outlier_cutoff,2.0)) {
                ctd = true;
                excluded[imeas] = true;
                nfit--;
            }
        }
    }

    return 0;

} // }}}

