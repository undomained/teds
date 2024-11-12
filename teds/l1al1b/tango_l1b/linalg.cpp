// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "linalg.h"

namespace tango {

auto dgemm(const char trans_a,
           const char trans_b,
           const int n_rows_a,
           const int n_cols_a,
           const std::vector<double>& a,
           const std::vector<double>& b,
           std::vector<double>& c) -> void
{
    // Throughout this function, a and b are swapped from the Fortran
    // perspective. Thus, trans_b is used where trans_a would normally
    // be used, n_rows_b instead of n_rows_a, and so on. If a and b
    // were not swapped, the function would yield c^T instead of c.
    constexpr double one { 1.0 };
    constexpr double zero { 0.0 };
    int n_rows_b { trans_a == 'n' ? n_cols_a : n_rows_a };
    int n_cols_b { static_cast<int>(b.size()) / n_rows_b };
    if (trans_b == 't') {
        std::swap(n_rows_b, n_cols_b);
    }
    const int M { trans_b == 'n' ? n_cols_b : n_rows_b };
    const int N { trans_a == 'n' ? n_rows_a : n_cols_a };
    const int K { trans_b == 'n' ? n_rows_b : n_cols_b };
    const int LDA { trans_b == 'n' ? M : K };
    const int LDB { trans_a == 'n' ? K : N };
    c.resize(M * N);
    dgemm_(&trans_b,
           &trans_a,
           &M,
           &N,
           &K,
           &one,
           b.data(),
           &LDA,
           a.data(),
           &LDB,
           &zero,
           c.data(),
           &M);
}

} // namespace tango
