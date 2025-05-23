// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "b_spline_2d.h"

#include "linalg.h"

namespace tango {

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

BSpline2D::BSpline2D(const int order,
                     const Eigen::ArrayXd& y_values_in,
                     const Eigen::ArrayXd& x_values_in,
                     const ArrayXXd& y_values_out,
                     const ArrayXXd& x_values_out)
  : y_values_out { y_values_out }, x_values_out { x_values_out }
{
    // STEP 1 - Find the X and Y matrices used for generating the
    //          control points.

    // Initialize 1D B-splines across rows and columns of the input grid
    b_spline_r = { order, y_values_in };
    b_spline_c = { order, x_values_in };
    // B-spline matrices across rows and columns
    Eigen::SparseMatrix<double> B_mat_r = b_spline_r.genBasis(y_values_in);
    Eigen::SparseMatrix<double> B_mat_c = b_spline_c.genBasis(x_values_in);
    // Find X from A^T A X = A^T where A is the B-spline matrix across
    // rows.
    Eigen::MatrixXd L = sparseToBanded(B_mat_r.transpose() * B_mat_r);
    Eigen::VectorXd D {};
    ldlt_decompose(L, D);
    Eigen::MatrixXd B_mat_r_full = B_mat_r.transpose();
    X = ldlt_solve(L, D, B_mat_r_full);
    // Find Y from A^T A Y = A^T where A is the B-spline matrix across
    // columns.
    L = sparseToBanded(B_mat_c.transpose() * B_mat_c);
    ldlt_decompose(L, D);
    Eigen::MatrixXd B_mat_c_full = B_mat_c.transpose();
    Y = ldlt_solve(L, D, B_mat_c_full);

    // STEP 2 - Precompute 1D spline values for evaluating the 2D
    //          spline. This saves time if the 2D spline is called
    //          multiple times but the input and output grids don't
    //          change.

    std::vector<double> b_spline_c_values(b_spline_r.getOrder() + 1);
    std::vector<double> b_spline_r_values(b_spline_r.getOrder() + 1);
    // Work array for evaluating B-spline at x for individual basis
    // functions, large enough for both 1D B-splines.
    std::vector<double> control_points_tmp(
      std::max(b_spline_r.nStates(), b_spline_c.nStates())
        + b_spline_r.getOrder() + 1,
      0.0);
    spline_values.resize(y_values_out.rows() * b_spline_r_values.size(),
                         y_values_out.cols() * b_spline_c_values.size());
#pragma omp parallel for firstprivate(                                         \
    b_spline_c_values, b_spline_r_values, control_points_tmp)
    for (long int i_row = 0; i_row < y_values_out.rows(); ++i_row) {
        for (long int i_col {}; i_col < y_values_out.cols(); ++i_col) {
            nonzeroBSplineValues(b_spline_r,
                                 y_values_out(i_row, i_col),
                                 control_points_tmp,
                                 b_spline_r_values);
            nonzeroBSplineValues(b_spline_c,
                                 x_values_out(i_row, i_col),
                                 control_points_tmp,
                                 b_spline_c_values);
            for (int i {}; i <= b_spline_r.getOrder(); ++i) {
                for (int j {}; j <= b_spline_c.getOrder(); ++j) {
                    spline_values(i_row * b_spline_r_values.size() + i,
                                  i_col * b_spline_c_values.size() + j) =
                      b_spline_r_values[i] * b_spline_c_values[j];
                }
            }
        }
    }
}

auto BSpline2D::genControlPoints(const Eigen::Ref<const ArrayXXd> data) -> void
{
    auto P = data.matrix().reshaped<Eigen::RowMajor>(X.cols(), Y.cols());
    control_points = X * P * Y.transpose();
}

auto BSpline2D::eval(Eigen::Ref<ArrayXXd> z) const -> void
{
    const int n_values_r { b_spline_r.getOrder() + 1 };
    const int n_values_c { b_spline_c.getOrder() + 1 };
    auto z_view =
      z.reshaped<Eigen::RowMajor>(y_values_out.rows(), y_values_out.cols());
    z_view = 0.0;
    for (long int i_row {}; i_row < y_values_out.rows(); ++i_row) {
        for (int long i_col {}; i_col < y_values_out.cols(); ++i_col) {
            const int i_r { b_spline_r.findInterval(
              y_values_out(i_row, i_col)) };
            const int i_c { b_spline_c.findInterval(
              x_values_out(i_row, i_col)) };
            for (int i {}; i < n_values_r; ++i) {
                for (int j {}; j < n_values_c; ++j) {
                    z_view(i_row, i_col) +=
                      control_points(i_r - i, i_c - j)
                      * spline_values(i_row * n_values_r + i,
                                      i_col * n_values_c + j);
                }
            }
        }
    }
}

} // namespace tango
