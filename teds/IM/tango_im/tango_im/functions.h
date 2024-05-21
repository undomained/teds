#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "header.h"

// Header of normal functions.

string now_timestring();
int juliandate(
    int day,
    int month,
    int year
);
string format(const string& fmt, ...);
int chol_solve(
    size_t nstate, // Size of the matrix (square).
    double *mat // Matrix (will be turned into her inverse. It is represented as just one triangle, because the matrix must be symmetric for this algorithm.
);

int powi(
    int base,
    unsigned int exp
);

void median(
    size_t sz,
    const double *arr,
    double &res
);
void fancymedian(
    size_t sz,
    const double *arr,
    double &res
);

// Linear interpolation (not using B-splines, because it is no rocket science).
// The independent values can be ascending or descending per array.
void linear_interpol(
    size_t sz1, // Size of the original array.
    size_t sz2, // Size of the array of interpolated values.
    double *z1, // Independent variable of original array.
    double *z2, // Independent variable of target array.
    double *f1, // Dependent variable of original array.
    double *f2 // Dependent variable of target array (output).
);
void fill_holes(
    size_t sz,
    vector<bool> &mask,
    double *arr
);

// Do a binary search to locate an index n such that x falls in the
// range list[n]...list[n+1].
[[nodiscard]] auto binaryFindIdx(const std::vector<double>& list,
                                 const double x) -> int;

#endif
