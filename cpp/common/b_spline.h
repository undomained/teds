// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for generating and evaluating 1D B-splines. Example usage:
//   // Initialize 3-order b-spline with a set of knots
//   BSpline spline { 3, knots };
//   // Construct the b-spline matrix with a set of data x-coordinates
//   Eigen::ArrayXd A {};
//   spline.evalBasis(x_values, A);
//   // Solve a linear system to find control points C from
//   // A C = P where P are the data y-coordinates (values)
//   ...
//   // Evaluate b-spline for some argument
//   spline.deBoor(C, 7.2);

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace tango {

class BSpline
{
private:
    // B-spline knots with both endpoints padded
    Eigen::ArrayXd knots {};
    // B-spline (polynomial) order
    int order {};
    // Work array for controlling which splines are evaluated. This is
    // used in the construction of the B-spline matrix.
    std::vector<double> control_points_tmp {};

public:
    BSpline() = default;
    BSpline(const int order, const Eigen::ArrayXd& knots);
    [[nodiscard]] auto getOrder() const -> int { return order; }
    // Return the size of the state vector for linear inversion. It is
    // (order+1) less than the number of knots because the endpoint
    // and the added knots should not be included.
    [[nodiscard]] auto nStates() const -> int;
    // Return index i such that x is in knots[i]..knots[i+1]
    auto findInterval(const double x) const -> int;
    // Evaluate the B-spline using de Boor algorithm for a given set
    // of control points and argument x.
    auto deBoor(const std::vector<double>& control_points,
                const double x) const -> double;
    // Construct B-spline matrix for a set of data x-values
    auto genBasis(const Eigen::ArrayXd& x_data) -> Eigen::SparseMatrix<double>;
    ~BSpline() = default;
};

} // namespace tango
