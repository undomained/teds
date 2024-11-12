// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "b_spline_2d.h"

#include "linalg.h"

namespace tango {

// Compute the solution of a system of linear equations ATA * X = AT,
// where ATA is a band matrix of order N with KL subdiagonals and KU
// superdiagonals, and X and AT are N-by-NRHS matrices. For B-splines,
// A is the B-spline matrix and the number of diagonals is the
// B-spline order. On exit, X^T (because of Fortran ordering) will be
// stored in A.
static auto banded_ATA_X_AT(const int n_diagonals,
                            const int n_rows_A,
                            std::vector<double>& A) -> void
{
    constexpr char uplo { 'u' };
    constexpr char trans { 'n' };
    constexpr double one { 1.0 };
    constexpr double zero { 0.0 };
    const int n_cols_A { static_cast<int>(A.size()) / n_rows_A };
    std::vector<double> ATA(n_cols_A * n_cols_A);
    dsyrk_(&uplo,
           &trans,
           &n_cols_A,
           &n_rows_A,
           &one,
           A.data(),
           &n_cols_A,
           &zero,
           ATA.data(),
           &n_rows_A);
    // Complete the upper (lower in Fortran) triangle
    for (int i {}; i < n_rows_A; ++i) {
        for (int k {}; k <= n_diagonals; ++k) {
            const int j { std::min(i + k, n_cols_A - 1) };
            ATA[i * n_cols_A + j] = ATA[j * n_cols_A + i];
        }
    }
    // Convert (A^T A) into a banded matrix
    const int KL { n_diagonals };
    const int KU { n_diagonals };
    const int n_rows_band { 2 * KL + KU + 1 };
    std::vector<double> ATA_banded(n_rows_band * n_cols_A);
    for (int i {}; i < n_cols_A; ++i) {
        for (int j { std::max(0, i - n_diagonals) };
             j <= std::min(i + n_diagonals, n_cols_A - 1);
             ++j) {
            constexpr int first_idx { 0 }; // C++: 0, Fortran: 1
            ATA_banded[(KL + KU + first_idx + i - j) * n_cols_A + j] =
              ATA[i * n_cols_A + j];
        }
    }
    // Transpose ATA_banded for dgbsv re-using ATA as a work array
    ATA.resize(n_cols_A * n_rows_band);
    for (int i {}; i < n_cols_A; ++i) {
        for (int j {}; j < n_rows_band; ++j) {
            ATA[i * n_rows_band + j] = ATA_banded[j * n_cols_A + i];
        }
    }
    ATA_banded = std::move(ATA);
    std::vector<int> ipiv(n_cols_A);
    int info {};
    dgbsv_(&n_cols_A,
           &n_diagonals,
           &n_diagonals,
           &n_rows_A,
           ATA_banded.data(),
           &n_rows_band,
           ipiv.data(),
           A.data(),
           &n_rows_A,
           &info);
}

BSpline2D::BSpline2D(const int order,
                     const std::vector<double>& x_values_r,
                     const std::vector<double>& x_values_c,
                     const std::vector<double>& data_in)
{
    // Initialize 1D B-splines across rows and columns of the input grid
    b_spline_r = { order, x_values_r };
    b_spline_c = { order, x_values_c };
    // B-spline matrix across rows
    std::vector<double> B_mat_r {};
    b_spline_r.evalBasis(x_values_r, B_mat_r);
    // B-spline matrix across columns
    std::vector<double> B_mat_c {};
    b_spline_c.evalBasis(x_values_c, B_mat_c);
    // Find X from A^T A X = A^T where A is the B-spline matrix across
    // rows. The result X^T is stored in B_mat_r.
    banded_ATA_X_AT(order, static_cast<int>(x_values_r.size()), B_mat_r);
    const std::vector<double>& XT { B_mat_r };
    // Find Y from B^T B Y = B^T where B is the B-spline matrix across
    // columns. The result Y^T is stored in B_mat_c.
    banded_ATA_X_AT(order, static_cast<int>(x_values_c.size()), B_mat_c);
    const std::vector<double>& YT { B_mat_c };
    // Compute P Y^T
    std::vector<double> PYT {};
    dgemm('n', 'n', x_values_r.size(), x_values_c.size(), data_in, YT, PYT);
    // Compute control points from Q = X P Y^T
    dgemm('t',
          'n',
          x_values_r.size(),
          b_spline_r.nStates(),
          XT,
          PYT,
          control_points);
}

// Given a target point x, find the knot interval i_x ... i_x+1 and
// only evaluate x for basis functions i_x-N ... i_x where N is the
// B-spline order. The rest of the basis functions do not contribute
// to the spline value at x.
static auto nonzeroBSplineValues(const BSpline& b_spline,
                                 const double x,
                                 std::vector<double>& control_points,
                                 std::vector<double>& b_spline_values) -> void
{
    const int i_x { b_spline.findInterval(x) };
    for (int k {}; k <= b_spline.getOrder(); ++k) {
        control_points[i_x - k] = 1.0;
        b_spline_values[k] = b_spline.deBoor(control_points, x);
        control_points[i_x - k] = 0.0;
    }
}

auto BSpline2D::eval(const std::vector<double>& x,
                     const std::vector<double>& y,
                     std::vector<double>& z) const -> void
{
    z.assign(x.size(), 0.0);
    std::vector<double> b_spline_c_values(b_spline_r.getOrder() + 1);
    std::vector<double> b_spline_r_values(b_spline_r.getOrder() + 1);
    // Work array for evaluating B-spline at x for individual basis
    // functions, large enough for both 1D B-splines.
    std::vector<double> control_points_tmp(
      std::max(b_spline_r.nStates(), b_spline_c.nStates())
        + b_spline_r.getOrder() + 1,
      0.0);
    for (int i_point {}; i_point < static_cast<int>(x.size()); ++i_point) {
        const int i_r { b_spline_r.findInterval(x[i_point]) };
        const int i_c { b_spline_c.findInterval(y[i_point]) };
        nonzeroBSplineValues(
          b_spline_r, x[i_point], control_points_tmp, b_spline_r_values);
        nonzeroBSplineValues(
          b_spline_c, y[i_point], control_points_tmp, b_spline_c_values);
        for (int i {}; i <= b_spline_r.getOrder(); ++i) {
            for (int j {}; j <= b_spline_r.getOrder(); ++j) {
                z[i_point] +=
                  control_points[(i_r - i) * b_spline_c.nStates() + i_c - j]
                  * b_spline_r_values[i] * b_spline_c_values[j];
            }
        }
    }
}

} // namespace tango
