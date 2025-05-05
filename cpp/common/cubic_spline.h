// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// A simple class for generating and evaluating cubic splines.
// Example usage:
//   CubicSpline spline { knots, values };
//   double val { spline.eval(3.7) };
// If the knots are equally spaced, the function evaluation has O(1)
// scaling, else O(log(N)) scaling.

#pragma once

#include <Eigen/Dense>

namespace tango {

class CubicSpline
{
private:
    bool equal_spacing { true };
    double range {};
    Eigen::ArrayXd knots {};
    Eigen::ArrayXd values {};
    Eigen::ArrayXd A {};
    Eigen::ArrayXd B {};
    Eigen::ArrayXd C {};
    [[nodiscard]] auto lookupIdx(const double x) const -> int;

public:
    CubicSpline(const Eigen::Ref<const Eigen::ArrayXd> x_values,
                const Eigen::Ref<const Eigen::ArrayXd> y_values);
    [[nodiscard]] auto eval(const double x) const -> double;
    [[nodiscard]] auto eval(const Eigen::ArrayXd& x) const -> Eigen::ArrayXd;
    [[nodiscard]] auto deriv(const double x) const -> double;
};

} // namespace tango
