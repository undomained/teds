// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for generating and evaluating 2D B-splines. Example usage:
//   // Initialize 3-order b-spline with two sets of knots, one in x
//   // and one in y direction, and a set of data points on a
//   // rectilinear grid defined by the knots.
//   BSpline2D spline { 3, knots_x, knots_y, data };
//   // Evaluate b-spline for some argument
//   spline.eval(arg_x, arg_y);

#pragma once

#include "b_spline.h"
#include "eigen.h"

namespace tango {

class BSpline2D
{
private:
    // 1D B-spline across rows of the input grid
    BSpline b_spline_r {};
    // 1D B-spline across columns of the input grid
    BSpline b_spline_c {};
    // Control points are given by C = X P Y^T where P are datapoints
    Eigen::MatrixXd X {};
    Eigen::MatrixXd Y {};
    // Precomputed 1D spline values (products)
    ArrayXXd spline_values {};
    // Control point matrix with dimensions Nr x Nc where Nr is the
    // number of B-spline states across rows (see BSpline::nStates).
    Eigen::MatrixXd control_points {};
    // Output grids (2D)
    ArrayXXd y_values_out {};
    ArrayXXd x_values_out {};

public:
    BSpline2D() = default;
    // Construct 1D B-splines and solve the 2D B-spline equation to
    // find the control points:
    //
    //   A^T A Q B^T B - A^T P B = 0,
    //
    // where A is the B-spline matrix across rows, B the B-spline
    // matrix across columns, P is a set of data points on the input
    // grid, and Q is a set of control points. The control points can
    // be expressed as
    //
    //   Q = X P Y^T,
    //   X = (A^T A)^-1 A^T,
    //   Y = (B^T B)^-1 B^T,
    //
    // where X and Y are found by solving
    //
    //   A^T A X = A^T,
    //   B^T B Y = B^T.
    //
    // Parameters
    // ----------
    // order
    //     B-spline order, same in both directions
    // x_values_in
    //     Input grid coordinates across rows, to be used as the
    //     B-spline knots.
    // y_values_in
    //     Input grid coordinates across columns, to be used as the
    //     B-spline knots.
    // x_values_out
    //     Output grid coordinates across rows
    // y_values_out
    //     Output grid coordinates across columns
    BSpline2D(const int order,
              const Eigen::ArrayXd& y_values_in,
              const Eigen::ArrayXd& x_values_in,
              const ArrayXXd& y_values_out,
              const ArrayXXd& x_values_out);
    // Computed control points using the X and Y matrices generated in
    // the constructor. The input data can change but not the
    // corresponding grids.
    //
    // Parameters
    // ----------
    // data
    //     Data values on input grid y_values_in and x_values_in used in
    //     the constructor.
    auto genControlPoints(const Eigen::Ref<const ArrayXXd> data) -> void;
    // Evaluate the 2D B-spline for points on a target grid (can be
    // irregular). The result for one target point can be expressed as
    //
    //   X(x,y) = Sum_i Sum_j N_i(x) N_j(y) Q_ij,
    //
    // where N_i(x) is the ith B-spline basis function across rows
    // evaluated at the row coordinate x, N_j(y) is the jth column
    // basis function evaluated at y, and Q_ij is a control point
    // corresponding to the grid point ij. While the sum formally runs
    // over all basis functions in both directions, in practice we
    // only evaluate a limited set of non-zero basis functions by
    // making use of de Boor's algorithm. The result array z should be
    // allocated and initialized outside this routine.
    auto eval(Eigen::Ref<ArrayXXd> z) const -> void;
};

} // namespace tango
