// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "bspline.h"

// Constructor.
Bspline::Bspline( // {{{
    Logger *creator, // Object that creates the B-spline.
    size_t a_order, // Order of the B-spline in Wikipedia perspective (n).
    size_t a_nknot, // Number of knot points.
    double *a_knots // Knot points in ascending order.
) : Logger(creator)
{
    // 1. Copy the arguments.
    order = a_order;
    nknot = a_nknot;
    knots.resize(nknot);
    memcpy(knots.data(),a_knots,nknot*sizeof(double));

    error = false;
    // Assert that you at least have two knots.
    if (nknot < 2) {error = true; return;}

    // Assert that the knots are ascending.
    for (size_t iknot=1 ; iknot<nknot ; iknot++) {
        if (knots[iknot] < knots[iknot-1]) {error = true; return;}
    }
    nspline = nknot + order - 2;

    polyfactors.resize((nknot-1)*order*order);
    vector<double> other_polyfactors((nknot-1)*(order-1)*(order-1)); // Large enough.
    double *poly_new;
    double *poly_old;
    if ((order & 1) == 1) {
        poly_new = polyfactors.data();
        poly_old = other_polyfactors.data();
    } else {
        poly_old = polyfactors.data();
        poly_new = other_polyfactors.data();
    }
    for (size_t ipoly=0 ; ipoly<nknot-1 ; ipoly++) poly_new[ipoly] = 1.0;
    for (size_t iorder=1 ; iorder<order ; iorder++) {
        double *sw = poly_old;
        poly_old = poly_new;
        poly_new = sw;
        for (size_t ipoly=0 ; ipoly<(nknot-1)*(iorder+1)*(iorder+1) ; ipoly++) poly_new[ipoly] = 0.0;
        double *poly_new_run = poly_new;
        double *poly_old_run = poly_old;
        for (size_t iint=0 ; iint<nknot-1 ; iint++) {
            for (size_t ispline=0 ; ispline<iorder ; ispline++) {
                // Our interval starts at iint and ends at iint+1.
                // For iorder=1, that is the interval. For higher iorders, the leftmost
                // interval expands to the left, so the left knot of the ramp is.
                // iint+1-iorder and the right knot is iint+1.
                // But then, we move ispline to the right, so we are at iint+1-iorder+ispline
                // to iint+1+ispline.
                size_t ileft = iint+1+ispline<iorder?0:iint+1-iorder+ispline;
                size_t iright = nknot-1<iint+1+ispline?nknot-1:iint+1+ispline;
                double left_x = 1.0 / (knots[ileft]-knots[iright]);
                double right_x = -left_x;
                double left_c = knots[iright] * right_x;
                double right_c = knots[ileft] * left_x;
                for (size_t iterm=0 ; iterm<iorder ; iterm++) {
                    poly_new_run[iterm] += left_c * poly_old_run[iterm];
                    poly_new_run[iterm+1] += left_x * poly_old_run[iterm];
                }
                poly_new_run += iorder+1;
                for (size_t iterm=0 ; iterm<iorder ; iterm++) {
                    poly_new_run[iterm] += right_c * poly_old_run[iterm];
                    poly_new_run[iterm+1] += right_x * poly_old_run[iterm];
                }
                poly_old_run += iorder;
            }
            poly_new_run += iorder+1;
        }
    }

} // }}}
Bspline::~Bspline() {}

// Full Jacobian formula.
int Bspline::jaccalc( // {{{
    size_t nx, // Number of evaluation points.
    double *x, // Evaluation points, sorted in ascending order, or at least sorted in ascending knot intervals.
    size_t *ideriv, // Array that indicates for each point which derivative to take, zero is just the function itself.
    double *mat // Resulting Jacobian matrix (if valid pointer is given).
)
{

    check_error(error,"B-spline not properly constructed.");

    for (size_t imat=0 ; imat<nx*nspline ; imat++) mat[imat] = 0.0;

    // Assert that the requested derivatives exist.
    for (size_t ix=0 ; ix<nx ; ix++) {
        check_error(ideriv[ix] >= order,"Error: Too low order for requested derivative.");
    }

    // Pre-calculate factorial terms.
    vector<double> facterm(order*order);
    double* fac_cur = facterm.data();
    vector<double> facs(order);
    facs[0] = 1.0;
    for (size_t ifac=1 ; ifac<order ; ifac++) facs[ifac] = ifac*facs[ifac-1];
    for (size_t iterm=0 ; iterm<order ; iterm++) {
        for (size_t ider=0 ; ider<=iterm ; ider++) {
            fac_cur[ider] = facs[iterm] / facs[iterm-ider];
        }
        fac_cur += order;
    }

    // Evaluate per interval, this assumes ascending x.
    size_t iint = 0;
    double *polyfactors_cur = polyfactors.data();
    for (size_t ix=0 ; ix<nx ; ix++) {
        while (x[ix] > knots[iint+1] && iint+2 != nknot) {
            iint++;
            mat += nx;
            polyfactors_cur += order*order;
        }
        while (x[ix] < knots[iint] && iint != 0) {
            iint--;
            mat -= nx;
            polyfactors_cur -= order*order;
        }
        // Evaluate x.
        for (size_t iterm=ideriv[ix] ; iterm<order ; iterm++) {
            double preterm = facterm[iterm*order+ideriv[ix]] * pow(x[ix],iterm-ideriv[ix]);
            for (size_t ispline=0 ; ispline<order ; ispline++) mat[ispline*nx+ix] += preterm * polyfactors_cur[ispline*order+iterm];
        }
    }

    return 0;

} // }}}

// Jacobian formula without ideriv.
int Bspline::jaccalc( // {{{
    size_t nx, // Number of evaluation points.
    double *x, // Evaluation points, sorted in ascending order, or at least sorted in ascending knot intervals.
    double *mat // Resulting Jacobian matrix (if valid pointer is given).
)
{

    check_error(error,"B-spline not properly constructed.");

    for (size_t imat=0 ; imat<nx*nspline ; imat++) mat[imat] = 0.0;

    // Evaluate per interval, this assumes ascending x.
    size_t iint = 0;
    double *polyfactors_cur = polyfactors.data();
    for (size_t ix=0 ; ix<nx ; ix++) {
        while (iint+2 != nknot && x[ix] > knots[iint+1]) {
            iint++;
            mat += nx;
            polyfactors_cur += order*order;
        }
        while (iint != 0 && x[ix] < knots[iint]) {
            iint--;
            mat -= nx;
            polyfactors_cur -= order*order;
        }
        // Evaluate x.
        for (size_t iterm=0 ; iterm<order ; iterm++) {
            double preterm = pow(x[ix],iterm);
            for (size_t ispline=0 ; ispline<order ; ispline++) mat[ispline*nx+ix] += preterm * polyfactors_cur[ispline*order+iterm];
        }
    }

    return 0;

} // }}}

// Application of full Jacobian formula.
int Bspline::jacapply( // {{{
    size_t nx, // Number of evaluation points.
    double *x, // Evaluation points, sorted in ascending order, or at least sorted in ascending knot intervals.
    size_t *ideriv, // Array that indicates for each point which derivative to take, zero is just the function itself.
    double *terms, // Terms with which to apply the Jacobian, saving the effort to save the Jacobian. You will get res = mat*terms then.
    double *res // Array where the result is written if you apply the Jacobian with res = mat*terms.
)
{

    check_error(error,"B-spline not properly constructed.");

    for (size_t ix=0 ; ix<nx ; ix++) res[ix] = 0.0;

    // Assert that the requested derivatives exist.
    for (size_t ix=0 ; ix<nx ; ix++) {
        check_error(ideriv[ix] >= order,"Error: Too low order for requested derivative.");
    }

    // Matrix multiply the terms with the polyfactors over the spline dimension.
    vector<double> polyterms((nknot-1)*order);
    for (size_t imat=0 ; imat<(nknot-1)*order ; imat++) polyterms[imat] = 0.0;
    for (size_t iint=0 ; iint<nknot-1 ; iint++) {
        for (size_t ispline=0 ; ispline<order ; ispline++) {
            for (size_t iterm=0 ; iterm<order ; iterm++) {
                polyterms[iint*order+iterm] += terms[ispline+iint] * polyfactors[iint*order*order+ispline*order+iterm];
            }
        }
    }

    // Pre-calculate factorial terms.
    vector<double> facterm(order*order);
    double* fac_cur = facterm.data();
    vector<double> facs(order);
    facs[0] = 1.0;
    for (size_t ifac=1 ; ifac<order ; ifac++) facs[ifac] = ifac*facs[ifac-1];
    for (size_t iterm=0 ; iterm<order ; iterm++) {
        for (size_t ider=0 ; ider<=iterm ; ider++) {
            fac_cur[ider] = facs[iterm] / facs[iterm-ider];
        }
        fac_cur += order;
    }

    // Evaluate per interval, this assumes ascending x.
    double *polyterms_cur = polyterms.data();
    size_t iint = 0;
    for (size_t ix=0 ; ix<nx ; ix++) {
        while (x[ix] > knots[iint+1] && iint+2 != nknot) {
            iint++;
            polyterms_cur += order;
        }
        while (x[ix] < knots[iint] && iint != 0) {
            iint--;
            polyterms_cur -= order;
        }
        // Evaluate x.
        for (size_t iterm=ideriv[ix] ; iterm<order ; iterm++) {
            res[ix] += facterm[iterm*order+ideriv[ix]] * polyterms_cur[iterm] * pow(x[ix],iterm-ideriv[ix]);
        }
    }

    return 0;

} // }}}

// Application without ideriv.
int Bspline::jacapply( // {{{
    size_t nx, // Number of evaluation points.
    double *x, // Evaluation points, sorted in ascending order, or at least sorted in ascending knot intervals.
    double *terms, // Terms with which to apply the Jacobian, saving the effort to save the Jacobian. You will get res = mat*terms then.
    double *res // Array where the result is written if you apply the Jacobian with res = mat*terms.
)
{

    check_error(error,"B-spline not properly constructed.");

    for (size_t ix=0 ; ix<nx ; ix++) res[ix] = 0.0;

    // Matrix multiply the terms with the polyfactors over the spline dimension.
    vector<double> polyterms((nknot-1)*order);
    for (size_t imat=0 ; imat<(nknot-1)*order ; imat++) polyterms[imat] = 0.0;
    for (size_t iint=0 ; iint<nknot-1 ; iint++) {
        for (size_t ispline=0 ; ispline<order ; ispline++) {
            for (size_t iterm=0 ; iterm<order ; iterm++) {
                polyterms[iint*order+iterm] += terms[ispline+iint] * polyfactors[iint*order*order+ispline*order+iterm];
            }
        }
    }

    // Evaluate per interval, this assumes ascending x.
    size_t iint = 0;
    double *polyterms_cur = polyterms.data();
    for (size_t ix=0 ; ix<nx ; ix++) {
        while (x[ix] > knots[iint+1] && iint+2 != nknot) {
            iint++;
            polyterms_cur += order;
        }
        while (x[ix] < knots[iint] && iint != 0) {
            iint--;
            polyterms_cur -= order;
        }
        // Evaluate x.
        for (size_t iterm=0 ; iterm<order ; iterm++) {
            res[ix] += polyterms_cur[iterm] * pow(x[ix],iterm);
        }
    }

    return 0;

} // }}}

