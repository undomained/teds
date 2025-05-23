// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "linear_spline.h"

#include "algorithm.h"

namespace tango {

// Threshold for determining if the knots are equally spaced
constexpr double equal_spacing_tol { 1e1
                                     * std::numeric_limits<double>::epsilon() };

LinearSpline::LinearSpline(const Eigen::ArrayXd& x_values,
                           const Eigen::ArrayXd& y_values)
  : knots { x_values }, values { y_values }
{
    // Determine if the knots are equally spaced
    double step0 { knots(1) - knots(0) };
    for (int i { 2 }; i < static_cast<int>(knots.size()); ++i) {
        const double step { knots(i) - knots(i - 1) };
        if (std::abs(step - step0) > equal_spacing_tol) {
            equal_spacing = false;
            break;
        }
    }
    if (equal_spacing) {
        range = knots(knots.size() - 1) - knots(0);
    }
    if (knots(1) < knots(0)) {
        std::ranges::reverse(knots);
        std::ranges::reverse(values);
    }
}

[[nodiscard]] auto LinearSpline::lookupIdx(const double x) const -> int
{
    return static_cast<int>((x - knots(0)) / range
                            * static_cast<double>(knots.size() - 1));
}

[[nodiscard]] auto LinearSpline::eval(const double x) const -> double
{
    int idx;
    if (x < knots(0)) {
        idx = 0;
    } else if (x > knots(knots.size() - 2)) {
        idx = static_cast<int>(knots.size() - 2);
    } else {
        idx = (equal_spacing ? lookupIdx(x) : binaryFindIdx(knots, x));
    }
    const double t { (x - knots(idx)) / (knots(idx + 1) - knots(idx)) };
    return std::lerp(values(idx), values(idx + 1), t);
}

[[nodiscard]] auto LinearSpline::eval(const Eigen::ArrayXd& x) const
  -> Eigen::ArrayXd
{
    Eigen::ArrayXd result(x.size());
    for (int i {}; i < static_cast<int>(x.size()); ++i) {
        result(i) = eval(x(i));
    }
    return result;
}

[[nodiscard]] auto LinearSpline::deriv(const double x) const -> double
{
    int idx;
    if (x < knots(0)) {
        idx = 0;
    } else if (x > knots(knots.size() - 2)) {
        idx = static_cast<int>(knots.size() - 2);
    } else {
        idx = (equal_spacing ? lookupIdx(x) : binaryFindIdx(knots, x));
    }
    return (values(idx + 1) - values(idx)) / (knots(idx + 1) - knots(idx));
}

[[nodiscard]] auto LinearSpline::deriv(const Eigen::ArrayXd& x) const
  -> Eigen::ArrayXd
{
    Eigen::ArrayXd result(x.size());
    for (int i {}; i < static_cast<int>(x.size()); ++i) {
        result(i) = deriv(x(i));
    }
    return result;
}

[[nodiscard]] auto LinearSpline::invert() const -> LinearSpline
{
    return LinearSpline(values, knots);
}

} // namespace tango
