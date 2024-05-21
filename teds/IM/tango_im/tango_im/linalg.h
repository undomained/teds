// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Determine which headers are included based on the linear library used.

#pragma once

#ifdef USE_MKL

#include <mkl_blas.h>
#include <mkl_lapack.h>

#else

extern "C" {
    auto dgttrf_(const int*,
                 double*,
                 double*,
                 double*,
                 double*,
                 int*,
                 int*) -> void;
    auto dgttrs_(const char*,
                 const int*,
                 const int*,
                 const double*,
                 const double*,
                 const double*,
                 const double*,
                 const int*,
                 double*,
                 const int*,
                 int*) -> void;
}

#endif
