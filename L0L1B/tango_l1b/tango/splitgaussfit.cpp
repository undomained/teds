// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "lininv.h"
#include "bspline.h"
#include "multigaussfit.h"
#include "splitgaussfit.h"

Splitgaussfit::Splitgaussfit( // {{{
    Logger *creator, // Object that creates this one.
    size_t a_ngauss, // Number of peaks (gausses).
    size_t a_npoly, // Number of polynomial parameters for the background.
    size_t a_nmeas // Number of elements in the measurement vector.
) : Logger(creator)
{

    // Copy dimensions.
    ngauss = a_ngauss;
    npoly = a_npoly;
    nmeas = a_nmeas;
    // Shape arrays.
    mat.resize(npoly*nmeas);
    x = 0; // Null-pointer for X. If the calculation is done with a null-pointer,
    // the default values are put there.

    // Set default mask. It is easier if it always exists.
    own_mask.resize(nmeas,false);
    mask = &own_mask; // Do not use data() for vector<bool>.

} // }}}
Splitgaussfit::~Splitgaussfit() {}

// Sets independent axis.
void Splitgaussfit::setX( // {{{
    double *a_x, // Independent axis array.
    bool copy // Flag for copying the input data. Leave false to pass pointer.
)
{

    if (copy) {
        own_x.resize(nmeas);
        x = own_x.data();
        memcpy(x,a_x,nmeas*sizeof(double));
    } else x = a_x;
    new_x = true;

} // }}}

void Splitgaussfit::setMeas( // {{{
    double *y,
    bool copy
)
{

    if (copy) {
        own_meas.resize(nmeas);
        meas = own_meas.data();
        // Copy array with the size that is already known..
        memcpy(meas,y,nmeas*sizeof(double));
    } else meas = y;

} // }}}

void Splitgaussfit::setNoise( // {{{
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
void Splitgaussfit::removeNoise( // {{{
)
{
    s_y.resize(0);
} // }}}

void Splitgaussfit::setMask( // {{{
    vector<bool> *a_mask,
    bool copy
)
{

    if (copy) {
        own_mask = *a_mask; // Copy of the vector<bool>.
        mask = &own_mask;
    } else mask = a_mask;

} // }}}
void Splitgaussfit::removeMask( // {{{
)
{
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) own_mask[imeas] = false;
    mask = &own_mask;
} // }}}

int Splitgaussfit::solve( // {{{
    int insignificant_island, // Islands of this many pixels are not recognized.
    double domain_range, // Size of a gauss tail to form domains.
    int maxreso, // Highest denominator in setting thresholds.
    double *res, // Retrieval result.
    bool debug // Logs debug information.
)
{

    // Null-pointer for x? Make new default value for x.
    if (x == 0) {
        own_x.resize(nmeas);
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) own_x[imeas] = (double)imeas;
        x = own_x.data();
        new_x = true;
    }

    // First, fit a polynomial. It will overestimate the background, because peaks
    // lift up the polynomial. But it is just an a-priori. It is actually a thing
    // that should be used as threshold to recognize peaks. Peaks become even worse,
    // so the inversion system still has work to do.
    if (new_x) {
        vector<double> knots = {x[0],x[nmeas-1]};
        Bspline b(this,npoly,2,knots.data()); // b.nspline = npoly.
        handle(b.jaccalc(nmeas,x,mat.data()));
        new_x = false;
    }
    vector<double> polypart(npoly);
    squareoption_t opt = s_y.size() == 0?OPT_NONE:OPT_DIAG;
    check_error(linear_invert(npoly,nmeas,mat.data(),opt,s_y.data(),mask,meas,polypart.data()) != 0,"Error: Initial background fit failed.");

    // Subtract this polynomial to get a residual that can be used to recognize peaks.
    vector<double> residual(nmeas);
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
        residual[imeas] = meas[imeas];
        for (size_t ipoly=0 ; ipoly<npoly ; ipoly++) residual[imeas] -= mat[ipoly*nmeas+imeas]*polypart[ipoly];
    }
    // Get signal median to base the convergence criteria.
    // The expectation is that the background goes to the median. The a-priori is
    // the mean, not the median, because that also works for higher-order polynomial
    // backgrounds. So the change during inversion is expected to be in the order of
    // the difference between mean and median, which is, if it is zeroth order, the
    // median of the residual.
    double med;
    median(nmeas,residual.data(),med);
    if (med >= 0.0) {
        if (debug) writelog(log_debug,"This seems to be a valley fitter.");
        return 3; // Median residual should be negative. Otherwise, it is not a peak fitter, but a valley fitter.
    }

    // Find the maximum to scale the peak threshold.
    double mx = 0.0;
    double mn = 0.0;
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
        if (residual[imeas] > mx) mx = residual[imeas];
        if (residual[imeas] < mn) mn = residual[imeas];
    }

    // We will consecutively try thresholds between minimum and maximum
    // at 1/2, 1/4, 3/4, 1/8, 1/3, 5/8, 7/8 etc until a maximum is
    // reached and then an error is raised.

    vector<bool> peak(nmeas); // This one should be available outside the while loop.
    int num = 0;
    int den = 1;
    while (true) { // This one will be broken at the limit or at a successful attempt.

        if (num+1 == den) {
            // Go to next denominator.
            num = 1;
            den *= 2;
            if (den > maxreso) return 4; // Unable to fit the right number of peaks.
        } else {
            // Go to next numerator.
            num += 2;
        }

        double threshold = (num*mx + (den-num)*mn) / den;

        // Peak flags.
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) peak[imeas] = residual[imeas] > threshold;

        if (insignificant_island > 0) {
            int unknown = 0;
            vector<int> cnt = {0,0};
            bool peak_cur = false;
            bool transition = false;
            for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
                bool p = peak[imeas]; // Or vector<bool>::reference p = peak[imeas]; bool &p = peak[imeas]; does not work.
                if (transition) {
                    cnt[!p] = 0;
                    cnt[p]++;
                    unknown++;
                    if (cnt[p] > insignificant_island) {
                        transition = false;
                        for (int iback=cnt[p] ; iback<unknown ; iback++) peak[imeas-iback] = p && peak_cur;
                        cnt[p] = 0;
                        unknown = 0;
                        peak_cur = p;
                    }
                } else {
                    if (p != peak_cur) {
                        transition = true;
                        cnt[p] = 1;
                        unknown = 1;
                    }
                }
            }
            if (transition) {
                for (int iback=0 ; iback<unknown ; iback++) peak[nmeas-1-iback] = false;
            }
        }

        // Count the number of peaks.
        size_t npeak_detected = 0;
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
            if (peak[imeas]) {
                if (imeas == 0) npeak_detected++;
                else if (!peak[imeas-1]) npeak_detected++;
            }
        }
        if (npeak_detected == ngauss) break;

    }

    // Set a-priori peak position right in the middle of the edges. Then, we will
    // fit some Gausses in an attempt to fit the entire system at once. That is
    // the cleanest way, but if it does not converge, we may face a challenge.
    // In such a case, we have to fall back to cutting our line in fragments.
    // But then, we still need the a-priori peak positions.

    size_t ipeak = 0;
    vector<double> peak_start(ngauss);
    vector<double> peak_end(ngauss);
    vector<size_t> idx_peak_start(ngauss);
    vector<size_t> idx_peak_end(ngauss);
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
        if (peak[imeas]) {
            if (imeas == 0) {
                peak_start[ipeak] = 1.5*x[imeas]-0.5*x[imeas+1];
                idx_peak_start[ipeak] = imeas;
            } else if (!peak[imeas-1]) {
                peak_start[ipeak] = 0.5*(x[imeas] + x[imeas-1]);
                idx_peak_start[ipeak] = imeas;
            }
            if (imeas == nmeas-1) {
                peak_end[ipeak] = 1.5*x[imeas]-0.5*x[imeas-1];
                idx_peak_end[ipeak] = imeas;
                ipeak++;
            } else if (!peak[imeas+1]) {
                peak_end[ipeak] = 0.5*(x[imeas] + x[imeas+1]);
                idx_peak_end[ipeak] = imeas;
                ipeak++;
            }
        }
    }

    vector<double> gausspart(3*ngauss);
    vector<double> gausspart_tol(3*ngauss);
    for (size_t ipeak=0 ; ipeak<ngauss ; ipeak++) {
        gausspart[3*ipeak] = residual[idx_peak_start[ipeak]];
        for (size_t imeas=idx_peak_start[ipeak]+1 ; imeas<=idx_peak_end[ipeak] ; imeas++) if (residual[imeas] > gausspart[3*ipeak]) gausspart[3*ipeak] = residual[imeas];
        gausspart[3*ipeak+1] = 0.5*(peak_start[ipeak]+peak_end[ipeak]);
        gausspart[3*ipeak+2] = (peak_end[ipeak]-peak_start[ipeak]) / (2.0*sqrt(2.0*log(2.0)));
        gausspart_tol[3*ipeak] = inv_magtol*gausspart[3*ipeak];
        gausspart_tol[3*ipeak+1] = inv_postol;
        gausspart_tol[3*ipeak+2] = inv_widthtol / (2.0*sqrt(2.0*log(2.0)));
    }

    if (!gaussfit) {
        memcpy(res,gausspart.data(),3*ngauss*sizeof(double));
        return 0;
    }

    // The polynomial part is still quite bad, because the peaks contribute to their
    // averages. The peak heights are also bad for the same reasons. The peak width
    // and positions are based on x-coordinates and should be fine.

    // A better apriori would be to lower the polynomial part and increase the
    // peak heights accordingly. This is possible with a linear fit on just the
    // peak intensities and the polynomial part, but that costs more computation
    // time than doing an additional iteration, because that linear fit is kind of
    // an iteration, but not using the splitting optimization.

    // Maybe, it is possible to evaulate the polynomial part and compare it to the
    // outermost points in each sub-region and transfer that difference from polynomial
    // part to peak strength. The most important issue is that it should stabilize
    // the inversion, preventing unlucky state changes for non-perfect Gausses.

    // So, some possible optimization implementation could be added.

    // Now, split these up into split gauss fits according to tail size.
    // Construct member multigaussfits.
    size_t ipeak_prv = 0;
    size_t imeas_start = 0;
    vector<double> tailsize(ngauss);
    for (size_t ipeak=0 ; ipeak<ngauss ; ipeak++) tailsize[ipeak] = domain_range * gausspart[3*ipeak+2] * (2.0*sqrt(2.0*log(2.0)));
    while (x[imeas_start] < gausspart[1]-tailsize[0]) imeas_start++;
    for (size_t ipeak=1 ; ipeak<ngauss+1 ; ipeak++) {
        // Construction to also close off the final instance.
        bool construct_new = ipeak == ngauss;
        if (!construct_new) construct_new = (gausspart[3*ipeak+1] - gausspart[3*(ipeak-1)+1] > tailsize[ipeak-1]+tailsize[ipeak]);
        if (construct_new) {
            // Append the multigaussfit that includes peaks ipeak_prv up to but excluding ipeak.
            size_t imeas_end = imeas_start;
            while (x[imeas_end] < gausspart[3*(ipeak-1)+1] + tailsize[ipeak-1]) {
                if (imeas_end == nmeas-1) break;
                imeas_end++;
            }
            size_t nmeas_cur = imeas_end-imeas_start+1;
            size_t npeak_cur = ipeak-ipeak_prv;
            Multigaussfit mgf_cur(this,npeak_cur,npoly,nmeas_cur);
            // Copy inversion settings.
            mgf_cur.setMaxiter(inv_maxiter);
            mgf_cur.setSteptolerance(inv_steptolerance);
            mgf_cur.setInitialStepReducer(inv_initial_step_reducer);
            mgf_cur.setReducerfactorFail(inv_reducerfactor_fail);
            mgf_cur.setReducerfactorSuccess(inv_reducerfactor_success);
            mgf_cur.setReducerLimit(inv_reducer_limit);
            // Cut out its meas, mat and x.
            mgf_cur.setMeas(&meas[imeas_start]);
            if (s_y.size() != 0) mgf_cur.setNoise(&s_y[imeas_start],true);
            vector<bool> mask_cur(nmeas_cur);
            for (size_t imeas_cur=0 ; imeas_cur<nmeas_cur ; imeas_cur++) {
                mask_cur[imeas_cur] = (*mask)[imeas_start+imeas_cur];
            }
            mgf_cur.setMask(&mask_cur);
            mgf_cur.setX(&x[imeas_start]);
            // Copy mat the easy way, rather than re-evaluating B-splines.
            for (size_t ipoly=0 ; ipoly<npoly ; ipoly++) memcpy(&mgf_cur.mat[ipoly*nmeas_cur],&mat[ipoly*nmeas+imeas_start],nmeas_cur*sizeof(double));
            // Convince multigaussfit that the mat is completed.
            mgf_cur.new_x = false;
            // Write the apriori.
            size_t nstate_cur = mgf_cur.getNstate();
            vector<double> apriori_cur(nstate_cur);
            memcpy(apriori_cur.data(),&gausspart[3*ipeak_prv],3*npeak_cur*sizeof(double));
            memcpy(&apriori_cur[3*npeak_cur],polypart.data(),npoly*sizeof(double));
            mgf_cur.setApriori(apriori_cur.data());
            vector<double> tolerance(nstate_cur);
            memcpy(tolerance.data(),&gausspart_tol[3*ipeak_prv],3*npeak_cur*sizeof(double));
            for (size_t istate=3*npeak_cur ; istate<nstate_cur ; istate++) tolerance[istate] = -inv_backgroundtol*med; // Minus sign because med is negative.
            mgf_cur.setTol(tolerance.data());
            // Progress running counters. Not if this was the final closure.
            if (ipeak != ngauss) {
                ipeak_prv = ipeak;
                imeas_start = imeas_end;
                while (x[imeas_start] < gausspart[3*ipeak+1] - tailsize[ipeak]) imeas_start++;
            }
            // Solve the inversion.
            vector<double> res_cur(nstate_cur);
            int conv = mgf_cur.inversion_calculate(res_cur.data(),debug);
            if (conv != 0) return conv;
            // Now, put the results in the correct places.
            // With some reverse engineering, we can memorize how many gauss parameters
            // included in this member.
            memcpy(res,res_cur.data(),3*npeak_cur*sizeof(double));
            res += 3*npeak_cur;
        }
    }

    return 0;

} // }}}

