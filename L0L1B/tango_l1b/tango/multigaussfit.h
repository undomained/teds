// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef MULTIGAUSSFIT_H
#define MULTIGAUSSFIT_H

#include "header.h"
#include "inversion.h"

// This is an inversion routine that attempts to fit a set of Gaussians.
// This is used to recognize a number of peaks in a line.
class Multigaussfit : public Inversion { // {{{

    // Constructor.
    public:
    Multigaussfit(
        Logger *creator, // Object that creates the fitter.
        size_t a_ngauss, // Number of peaks (gausses).
        size_t a_npoly, // Number of polynomial parameters for the background.
        size_t a_nmeas // Number of elements in the measurement vector.
    );
    ~Multigaussfit(); // Destructor.
    // TODO: Do this somehow a bit differently.
    // Maybe, the multigaussfit should have a private constructor and have
    // friended classes construct instances and access new_x or completely
    // keep track of this optimization.
    public:
    bool new_x;
    vector<double> mat;
    private:
    size_t ngauss; // Number of peaks.
    size_t npoly; // Number of polynomial terms, one higher than the order.
    vector<double> own_x;
    double *x;

    // Forward model to override.
    int fwd(
        double *state, // State vector.
        double *signal, // Modelled signal (outout).
        double *jacobian // Modelled jacobian (outout).
    ) override;

    // Setter.
    public:
    void setX(
        double *a_x, // Independent parameter of the line.
        bool copy=false // Flag for copying the argument.
    );

}; // }}}

#endif
