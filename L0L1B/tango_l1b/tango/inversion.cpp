// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "matrix.h"
#include "lininv.h"
#include "inversion.h"

Inversion::Inversion( // {{{
    Logger *creator, // Object that creates this instance.
    size_t a_nstate, // Number of state parameters.
    size_t a_nmeas // Number of elements in measurement vector.
) : Logger(creator)
{

    // Copy dimension size.
    nstate = a_nstate;
    nmeas = a_nmeas;

    // The mask is the only one that has a default.
    own_mask.resize(nmeas,false);
    mask = &own_mask; // Do not use data() for vector<bool>.

} // }}}

// Setters.
void Inversion::setApriori( // {{{
    double *x,
    bool copy
)
{
    if (copy) {
        own_apriori.resize(nstate);
        apriori = own_apriori.data();
        // Copy array with the size that is already known.
        memcpy(apriori,x,nstate*sizeof(double));
    } else apriori = x;

} // }}}
void Inversion::setMeas( // {{{
    double *y,
    bool copy
)
{

    if (copy) {
        own_meas.resize(nmeas);
        meas = own_meas.data();
        // Copy array with the size that is already known.
        memcpy(meas,y,nmeas*sizeof(double));
    } else meas = y;

} // }}}
void Inversion::setNoise( // {{{
    double *y_noise,
    bool s
)
{

    s_y.resize(nmeas);
    if (s) {
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) s_y[imeas] = y_noise[imeas];
    } else {
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) s_y[imeas] = pow(y_noise[imeas],2.0);
    }

} // }}}
void Inversion::removeNoise( // {{{
)
{
    s_y.resize(0);
} // }}}
void Inversion::setMask( // {{{
    vector<bool> *a_mask,
    bool copy
)
{

    if (copy) {
        own_mask = *a_mask; // Copy of the vector<bool>.
        mask = &own_mask;
    } else mask = a_mask;

} // }}}
void Inversion::removeMask( // {{{
)
{
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) own_mask[imeas] = false;
    mask = &own_mask;
} // }}}
void Inversion::setTol( // {{{
    double *t,
    bool copy
)
{

    if (copy) {
        own_tol.resize(nstate);
        tol = own_tol.data();
        // Copy array with the size that is already known.
        memcpy(tol,t,nstate*sizeof(double));
    } else tol = t;

} // }}}

// Perform the inversion.
int Inversion::inversion_calculate( // {{{
    double *res, // Retrieved state.
    bool debug // Logs debugging information.
)
{

    // Verify that the inversion is not ill-posed.
    check_error(nmeas < nstate,"Cannot fit %zu parameters with only %zu measurement.",nstate,nmeas);

    // Construct local variables.
    vector<double> mod(nmeas); // Modelled signal.
    vector<double> residual(nmeas); // Measured minus model.
    vector<double> jac(nmeas*nstate); // Jacobian.
    vector<double> state(nstate); // Current state.
    vector<double> state_des(nstate); // Desired state change.

    // There may be an outlier protection mechanism in the linear inverison.
    // In that case, we keep track of all the measurement points that are
    // excluded on the way. It starts with the input pixel mask.
    vector<bool> excluded(*mask);

    // Go to a-priori state.
    memcpy(state.data(),apriori,nstate*sizeof(double));

    if (debug) {
        writelog(log_debug,"A-priori");
        string line = "";
        for (size_t istate=0 ; istate<nstate ; istate++) line += format("%21.12f",state[istate]);
        writelog(log_debug,"%s",line.c_str());
    }

    // Run forward model.
    fwd(state.data(),mod.data(),jac.data());
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) residual[imeas] = meas[imeas]-mod[imeas];

    int conv = 1; // Not converged.
    double step_reducer = initial_step_reducer;
    for (size_t iter=0 ; iter < maxiter ; iter++) {

        // Perform linear inversion algorithm.
        int ierr;
        ierr = execute_linear_invert(jac.data(),residual.data(),&excluded,state_des.data());
        if (ierr != 0) return 2; // Numeric error occurred. Not a lethal program error.
        // The state_des is now a movement.

        // Calculate initial chi squared using the mask provided by the inversion..
        double chi2_before = 0.0;
        if (s_y.size() == 0) {
            for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                if (!excluded[imeas]) chi2_before += pow(residual[imeas],2.0);
            }
        } else {
            for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                if (!excluded[imeas]) chi2_before += pow(residual[imeas],2.0)/s_y[imeas]; // In this inversion, the s_y can only be diagonal.
            }
        }
        // No need to reduce.
        if (debug && iter == 0) writelog(log_debug,"Initial unreduced chi square: %.12f",chi2_before);

        // Apply retrieval to state.
        for (size_t istate=0 ; istate<nstate ; istate++) state[istate] += (1.0/(1.0+step_reducer)) * state_des[istate];
        if (debug) {
            writelog(log_debug,"Attempted state");
            string line = "";
            for (size_t istate=0 ; istate<nstate ; istate++) line += format("%21.12f",state[istate]);
            writelog(log_debug,"%s",line.c_str());
        }

        bool acc = true;
        while (true) { // Step control loop.

            // Rerun forward model.
            fwd(state.data(),mod.data(),jac.data());
            for (size_t imeas=0 ; imeas<nmeas ; imeas++) residual[imeas] = meas[imeas]-mod[imeas];

            // Calculate new chi squared using the same mask as during the inversion.
            double chi2_after = 0.0;
            if (s_y.size() == 0) {
                for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                    if (!excluded[imeas]) chi2_after += pow(residual[imeas],2.0);
                }
            } else {
                for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                    if (!excluded[imeas]) chi2_after += pow(residual[imeas],2.0)/s_y[imeas];
                }
            }
            if (debug) writelog(log_debug,"Unreduced chi square: %.12f",chi2_after);

            if (isnan(chi2_after)) return 2; // Numeric error occurred. Not a lethal program error.
            if (steptolerance == NC_FILL_DOUBLE || chi2_after <= chi2_before*steptolerance) {
                step_reducer *= reducerfactor_success;
                break;
            }
            if (debug) writelog(log_debug,"Reject step.");

            acc = false; // First step was not accepted. No convergence can be triggered at this iteration.
            for (size_t istate=0 ; istate<nstate ; istate++) {
                // Undo step.
                state[istate] -= (1.0/(1.0+step_reducer)) * state_des[istate];
            }
            // Update reducer.
            step_reducer *= reducerfactor_fail;
            for (size_t istate=0 ; istate<nstate ; istate++) {
                state[istate] += (1.0/(1.0+step_reducer)) * state_des[istate];
            }
            if (debug) {
                writelog(log_debug,"New attempt");
                string line = "";
                for (size_t istate=0 ; istate<nstate ; istate++) line += format("%21.12f",state[istate]);
                writelog(log_debug,"%s",line.c_str());
            }

        }
        if (debug) writelog(log_debug,"Accept.");
        // Check for convergence.
        if (acc && (reducer_limit == NC_FILL_DOUBLE || step_reducer < reducer_limit)) { // Direct accept.
            double mag_diff_sq = 0.0;
            // The change of the state is (still) in state_des. It can be reduced by
            // a step reduction factor, but then the boolean acc had been false and
            // we would not be here.
            for (size_t istate=0 ; istate<nstate ; istate++) mag_diff_sq += pow(state_des[istate]/tol[istate],2.0);
            if (mag_diff_sq < 1.0) {
                conv = 0;
                break;
            }
        }

    }

    // Copy state to result.
    for (size_t istate=0 ; istate<nstate ; istate++) res[istate] = state[istate];

    // This is zero if the loop is terminated by convergence, and one if the loop
    // was terminated by iterations.
    return conv;

} // }}}

// Default linear inversion. Can be overwritten by child classes.
int Inversion::execute_linear_invert( // {{{
    double *jac,
    double *residual,
    vector<bool> *excluded,
    double *state_des
)
{
    squareoption_t opt_s_y = s_y.size() == 0?OPT_NONE:OPT_DIAG;
    handle_nonlethal(linear_invert(nstate,nmeas,jac,opt_s_y,s_y.data(),excluded,residual,state_des));

    return 0;

} // }}}

