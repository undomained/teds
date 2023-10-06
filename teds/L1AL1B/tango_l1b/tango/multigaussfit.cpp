// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "bspline.h"
#include "multigaussfit.h"

Multigaussfit::Multigaussfit( // {{{
    Logger *creator, // Object that creates the fitter.
    size_t a_ngauss, // Number of peaks (gausses).
    size_t a_npoly, // Number of polynomial parameters for the background.
    size_t a_nmeas // Number of elements in the measurement vector.
) : Inversion(creator,3*a_ngauss+a_npoly,a_nmeas)
{

    // Copy dimensions.
    ngauss = a_ngauss;
    npoly = a_npoly;
    // Shape arrays.
    mat.resize(npoly*nmeas);
    new_x = true;

} // }}}
Multigaussfit::~Multigaussfit() {}

int Multigaussfit::fwd( // {{{
    double *state, // State vector.
    double *signal, // Modelled signal (outout).
    double *jacobian // Modelled jacobian (outout).
)
{

    // Initialize with zeros.
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) signal[imeas] = 0.0;

    // Running pointers.
    double *st = state;
    double *jac = jacobian;

    // Gauss part.
    for (size_t igauss=0 ; igauss<ngauss ; igauss++) {
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
            jac[0*nmeas+imeas] = exp(-(pow(x[imeas]-st[1],2) / (2.0*pow(st[2],2))));
            signal[imeas] += st[0] * jac[0*nmeas+imeas];
            jac[1*nmeas+imeas] = st[0] * jac[0*nmeas+imeas] * (x[imeas]-st[1]) / pow(st[2],2);
            jac[2*nmeas+imeas] = st[0] * jac[0*nmeas+imeas] * (pow(x[imeas]-st[1],2) / (pow(st[2],3)));
        }
        st += 3;
        jac += 3*nmeas;
    }

    // Calculate mat if necessary.
    if (new_x) {
        vector<double> knots = {x[0],x[nmeas-1]};
        Bspline b(this,npoly,2,knots.data());
        handle(b.jaccalc(nmeas,x,mat.data()));
        new_x = false;
    }
    // Copy it to the Jacobian and apply it to the measurement.
    memcpy(jac,mat.data(),npoly*nmeas*sizeof(double));
    for (size_t ipoly=0 ; ipoly<npoly ; ipoly++) {
        for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
            signal[imeas] += st[ipoly] * mat[ipoly*nmeas+imeas];
        }
    }

    return 0;

} // }}}

// Sets independent axis.
void Multigaussfit::setX( // {{{
    double *a_x, // Independent axis array.
    bool copy
)
{

    if (copy) {
        own_x.resize(nmeas);
        x = own_x.data();
        memcpy(x,a_x,nmeas*sizeof(double));
    } else x = a_x;
    new_x = true;

} // }}}

