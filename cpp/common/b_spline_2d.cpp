// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "b_spline_2d.h"

#include "linalg.h"

namespace tango {

BSpline2D::BSpline2D(const int order,
                     const Eigen::ArrayXd& x_values_r,
                     const Eigen::ArrayXd& x_values_c,
                     const Eigen::Ref<const ArrayXXd> data_in)
{
    // Initialize 1D B-splines across rows and columns of the input grid
    b_spline_r = { order, x_values_r };
    b_spline_c = { order, x_values_c };
    // B-spline matrices across rows and columns
    Eigen::SparseMatrix<double> B_mat_r(b_spline_r.genBasis(x_values_r));
    Eigen::SparseMatrix<double> B_mat_c(b_spline_c.genBasis(x_values_c));
    // Find X from A^T A X = A^T where A is the B-spline matrix across
    // rows.
    Eigen::MatrixXd L(sparseToBanded(B_mat_r.transpose() * B_mat_r));
    Eigen::VectorXd D {};
    ldlt_decompose(L, D);
    Eigen::MatrixXd B_mat_r_full(B_mat_r.transpose());
    auto X = ldlt_solve(L, D, B_mat_r_full);
    // Find Y from A^T A Y = A^T where A is the B-spline matrix across
    // columns.
    L = sparseToBanded(B_mat_c.transpose() * B_mat_c);
    ldlt_decompose(L, D);
    Eigen::MatrixXd B_mat_c_full(B_mat_c.transpose());
    auto Y = ldlt_solve(L, D, B_mat_c_full);
    // Control points are given by X P Y^T where P are the datapoints
    auto P(data_in.matrix().reshaped<Eigen::RowMajor>(x_values_r.size(),
                                                      x_values_c.size()));
    control_points = X * P * Y.transpose();
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

auto BSpline2D::eval(const ArrayXXd& x,
                     const ArrayXXd& y,
                     Eigen::Ref<ArrayXXd> z) const -> void
{
    std::vector<double> b_spline_c_values(b_spline_r.getOrder() + 1);
    std::vector<double> b_spline_r_values(b_spline_r.getOrder() + 1);
    // Work array for evaluating B-spline at x for individual basis
    // functions, large enough for both 1D B-splines.
    std::vector<double> control_points_tmp(
      std::max(b_spline_r.nStates(), b_spline_c.nStates())
        + b_spline_r.getOrder() + 1,
      0.0);
    auto z_view(z.reshaped<Eigen::RowMajor>(x.rows(), x.cols()));
    z_view = 0.0;
    for (int i_row {}; i_row < static_cast<int>(x.rows()); ++i_row) {
        for (int i_col {}; i_col < static_cast<int>(x.cols()); ++i_col) {
            const int i_r { b_spline_r.findInterval(x(i_row, i_col)) };
            const int i_c { b_spline_c.findInterval(y(i_row, i_col)) };
            nonzeroBSplineValues(b_spline_r,
                                 x(i_row, i_col),
                                 control_points_tmp,
                                 b_spline_r_values);
            nonzeroBSplineValues(b_spline_c,
                                 y(i_row, i_col),
                                 control_points_tmp,
                                 b_spline_c_values);
            for (int i {}; i <= b_spline_r.getOrder(); ++i) {
                for (int j {}; j <= b_spline_r.getOrder(); ++j) {
                    z_view(i_row, i_col) += control_points(i_r - i, i_c - j)
                                            * b_spline_r_values[i]
                                            * b_spline_c_values[j];
                }
            }
        }
    }
}

} // namespace tango
