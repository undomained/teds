#pragma once

#include "header.h"
#include "logger.h"

namespace tango {

class Bspline : public Logger {

    private:
    // Passed in constructor.
    size_t order; // Order in Wikipedia perspective (n), so n=1 is order-zero polynomial.
    size_t nknot; // Number of knots.
    vector<double> knots; // Position of the knots in ascending order.

    bool error = true; // Flag for proper construction.

    public:
    vector<double> polyfactors; // Pre-calculated prefactors for the Jacobian.

    public:
    Bspline(
        Logger *creator,
        size_t a_order,
        size_t a_nknot,
        double *knots
    ); // Constructor.
    ~Bspline(); // Destructor.
    size_t nspline; // Number of B-splines (N+n-2). This is the number of state parameters if you do an inversion.
    // Routines to directly apply a Jacobian. This does not return the Jacobian
    // matrix.
    int jacapply(size_t nx, double *x, double *terms, double *res); // Without ideriv, so ideriv is considered 0 everywhere.
    int jacapply(size_t nx, double *x, size_t *ideriv, double *terms, double *res); // Full implementation.
    // Routines for returninng the Jacobian, for instance, for inversion.
    int jaccalc(size_t nx, double *x, double *mat); // Without ideriv, so ideriv is considered 0 everywhere.
    int jaccalc(size_t nx, double *x, size_t *ideriv, double *mat); // Full implementation.

};

} // namespace tango
