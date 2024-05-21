// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// A simple class for generating and evaluating linear splines.
// Example usage:
//   LinearSpline spline { knots, values };
//   double val { spline.eval(3.7) };
// If the knots are equally spaced, the function evaluation has O(1)
// scaling, else O(log(N)) scaling.

#pragma once

#include <vector>

namespace tango {

class LinearSpline
{
private:
    bool equal_spacing { true };
    double range {};
    std::vector<double> knots {};
    std::vector<double> values {};
    // Determine if the knots are equally spaced
    auto setSpacing() -> void;
    auto lookupIdx(const double x) const -> int;

public:
    LinearSpline(const std::vector<double>& x_values,
                 const std::vector<double>& y_values);
    auto eval(const double x) const -> double;
};

} // namespace tango
