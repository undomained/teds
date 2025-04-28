// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// A simple class for generating and evaluating linear splines and
// their derivatives.
// Example usage:
//   LinearSpline spline { knots, values };
//   double val { spline.eval(3.7) };
// If the knots are equally spaced it has O(1) scaling, else O(log(N))
// scaling.

#pragma once

#include <Eigen/Dense>

namespace tango {

class LinearSpline
{
private:
    bool equal_spacing { true };
    double range {};
    Eigen::ArrayXd knots {};
    Eigen::ArrayXd values {};
    [[nodiscard]] auto lookupIdx(const double x) const -> int;

public:
    LinearSpline() = default;
    LinearSpline(const Eigen::ArrayXd& x_values,
                 const Eigen::ArrayXd& y_values);
    [[nodiscard]] auto eval(const double x) const -> double;
    [[nodiscard]] auto eval(const Eigen::ArrayXd& x) const -> Eigen::ArrayXd;
    [[nodiscard]] auto deriv(const double x) const -> double;
    [[nodiscard]] auto deriv(const Eigen::ArrayXd& x) const -> Eigen::ArrayXd;
    // Construct a new spline by inverting the x and y values
    [[nodiscard]] auto invert() const -> LinearSpline;
    ~LinearSpline() = default;
};

} // namespace tango
