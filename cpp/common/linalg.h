// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Determine which headers are included based on the linear library used.

#pragma once

#include <vector>

extern "C"
{
    auto dgttrf_(const int*, double*, double*, double*, double*, int*, int*)
      -> void;
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
    auto dgbsv_(const int*,
                const int*,
                const int*,
                const int*,
                double*,
                const int*,
                int*,
                double*,
                const int*,
                int*) -> void;
    auto dsyrk_(const char* uplo,
                const char* trans,
                const int* n,
                const int* k,
                const double* alpha,
                const double* a,
                const int* lda,
                const double* beta,
                double* c,
                const int* ldc) -> void;
    auto dgemm_(const char* transa,
                const char* transb,
                const int* m,
                const int* n,
                const int* k,
                const double* alpha,
                const double* a,
                const int* lda,
                const double* b,
                const int* ldb,
                const double* beta,
                double* c,
                const int* ldc) -> void;
}

namespace tango {

// Wrapper around the Fortran dgemm routine for matrix
// multiplication. The interface assumes row-major ordering or arrays.
// Inputs:
//   trans_a - whether matrix 'a' should be transposed ('n' or 't')
//   trans_b - whether matrix 'b' should be transposed
//   n_rows_a - number of rows of matrix A
//   n_cols_a - number of columns of matrix A
auto dgemm(const char trans_a,
           const char trans_b,
           const int n_rows_a,
           const int n_cols_a,
           const double* a,
           const std::vector<double>& b,
           std::vector<double>& c) -> void;
auto dgemm(const char trans_a,
           const char trans_b,
           const int n_rows_a,
           const int n_cols_a,
           const std::vector<double>& a,
           const std::vector<double>& b,
           std::vector<double>& c) -> void;

} // namespace tango
