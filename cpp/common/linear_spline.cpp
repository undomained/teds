// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "linear_spline.h"

#include "algorithm.h"
#include "constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace tango {

// Threshold for determining if the knots are equally spaced
constexpr double equal_spacing_tol { 1e1
                                     * std::numeric_limits<double>::epsilon() };

LinearSpline::LinearSpline(const std::vector<double>& x_values,
                           const std::vector<double>& y_values)
  : knots { x_values }, values { y_values }
{
    // Determine if the knots are equally spaced
    double step0 { knots[1] - knots.front() };
    for (int i { 2 }; i < static_cast<int>(knots.size()); ++i) {
        const double step { knots[i] - knots[i - 1] };
        if (std::abs(step - step0) > equal_spacing_tol) {
            equal_spacing = false;
            break;
        }
    }
    if (equal_spacing) {
        range = knots.back() - knots.front();
    }
    if (knots[1] < knots[0]) {
        std::ranges::reverse(knots);
        std::ranges::reverse(values);
    }
}

[[nodiscard]] auto LinearSpline::lookupIdx(const double x) const -> int
{
    return static_cast<int>((x - knots.front()) / range
                            * static_cast<double>(knots.size() - 1));
}

[[nodiscard]] auto LinearSpline::eval(const double x) const -> double
{
    const int idx { equal_spacing ? lookupIdx(x) : binaryFindIdx(knots, x) };
    if (idx == fill::i) {
        if (x < knots.front()) {
            const double t { (x - knots.front()) / (knots[1] - knots.front()) };
            return std::lerp(values.front(), values[1], t);
        }
        const double& knots_2 { knots[knots.size() - 2] };
        const double t { (x - knots_2) / (knots.back() - knots_2) };
        return std::lerp(values[knots.size() - 2], values.back(), t);
    }
    const double t { (x - knots[idx]) / (knots[idx + 1] - knots[idx]) };
    return std::lerp(values[idx], values[idx + 1], t);
}

[[nodiscard]] auto LinearSpline::deriv(const double x) const -> double
{
    const int idx { equal_spacing ? lookupIdx(x) : binaryFindIdx(knots, x) };
    if (idx == fill::i) {
        if (x < knots.front()) {
            return (values[1] - values.front()) / (knots[1] - knots.front());
        }
        return (values.back() - values[knots.size() - 2])
               / (knots.back() - knots[knots.size() - 2]);
    }
    return (values[idx + 1] - values[idx]) / (knots[idx + 1] - knots[idx]);
}

[[nodiscard]] auto LinearSpline::invert() const -> LinearSpline
{
    return { values, knots };
}

} // namespace tango
