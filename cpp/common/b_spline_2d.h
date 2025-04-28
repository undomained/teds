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
    // Control point matrix with dimensions Nr x Nc where Nr is the
    // number of B-spline states across rows (see BSpline::nStates).
    Eigen::MatrixXd control_points {};

public:
    BSpline2D() = default;
    // Construct the 1D B-splines and solve the 2D B-spline equation
    // to find the control points:
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
    // Inputs:
    //   x_values_r - input grid coordinates across rows, to be used
    //                as the B-spline knots.
    //   x_values_c - the input grid coordinates across columns, to be
    //                used as the B-spline knots.
    //   data_in - the data values on that grid.
    BSpline2D(const int order,
              const Eigen::ArrayXd& x_values_r,
              const Eigen::ArrayXd& x_values_c,
              const Eigen::Ref<const ArrayXXd> data_in);
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
    auto eval(const ArrayXXd& x,
              const ArrayXXd& y,
              Eigen::Ref<ArrayXXd> z) const -> void;
};

} // namespace tango
