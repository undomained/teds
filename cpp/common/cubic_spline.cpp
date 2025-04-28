// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "cubic_spline.h"

#include "algorithm.h"

namespace tango {

// Threshold for determining if the knots are equally spaced
constexpr double equal_spacing_tol { 1e1
                                     * std::numeric_limits<double>::epsilon() };

CubicSpline::CubicSpline(const Eigen::Ref<const Eigen::ArrayXd> x_values,
                         const Eigen::Ref<const Eigen::ArrayXd> y_values)
  : knots { x_values }, values { y_values }
{
    // Determine if the knots are equally spaced
    double step0 { knots(1) - knots(0) };
    for (int i { 2 }; i < knots.rows(); ++i) {
        const double step { knots(i) - knots(i - 1) };
        if (std::abs(step - step0) > equal_spacing_tol) {
            equal_spacing = false;
            break;
        }
    }
    if (equal_spacing) {
        range = knots(knots.size() - 1) - knots(0);
    }
    const Eigen::ArrayXd delta_x(knots(Eigen::seq(1, Eigen::last))
                                 - knots(Eigen::seqN(0, knots.size() - 1)));
    const Eigen::ArrayXd delta_y(y_values(Eigen::seq(1, Eigen::last))
                                 - y_values(Eigen::seqN(0, knots.size() - 1)));
    Eigen::ArrayXd diag(knots.size());
    diag(0) = 1.0;
    diag(diag.size() - 1) = 1.0;
    diag(Eigen::seq(1, Eigen::last - 1)) =
      2.0 / 3.0
      * (delta_x(Eigen::seqN(0, knots.size() - 2))
         + delta_x(Eigen::seqN(1, knots.size() - 2)));
    Eigen::ArrayXd diag_upper = delta_x / 3.0;
    Eigen::ArrayXd diag_lower = diag_upper;
    diag_upper(0) = 0.0;
    diag_lower(diag_lower.size() - 1) = 0.0;
    Eigen::ArrayXd rhs(knots.size());
    const auto n2 { rhs.size() - 2 };
    rhs(Eigen::seqN(1, n2)) =
      delta_y(Eigen::seqN(1, n2)) / delta_x(Eigen::seqN(1, n2))
      - delta_y(Eigen::seqN(0, n2)) / delta_x(Eigen::seqN(0, n2));
    rhs(0) = 0.0;
    rhs(rhs.size() - 1) = 0.0;
    // Solve the tridiagonal system that consists of the diagonals and
    // the right-hand side using the Thomas algoritm.
    for (int i { 1 }; i < rhs.size(); ++i) {
        const double w { diag_lower(i - 1) / diag(i - 1) };
        diag(i) -= w * diag_upper(i - 1);
        rhs(i) -= w * rhs(i - 1);
    }
    B.resize(rhs.size());
    B(B.size() - 1) = rhs(rhs.size() - 1) / diag(diag.size() - 1);
    for (int i { static_cast<int>(B.size()) - 2 }; i >= 0; --i) {
        B(i) = (rhs(i) - diag_upper(i) * B(i + 1)) / diag(i);
    }
    // Compute A and C from B
    const auto n1 { B.size() - 1 };
    A = delta_y / delta_x
        - delta_x * (2 * B(Eigen::seqN(0, n1)) + B(Eigen::seqN(1, n1))) / 3.0;
    C = (B(Eigen::seqN(1, n1)) - B(Eigen::seqN(0, n1))) / (3.0 * delta_x);
}

[[nodiscard]] auto CubicSpline::lookupIdx(const double x) const -> int
{
    return std::min(
      std::max(0,
               static_cast<int>((x - knots(0)) / range
                                * (static_cast<double>(knots.size() - 1)))),
      static_cast<int>(knots.size() - 1));
}

[[nodiscard]] auto CubicSpline::eval(const double x) const -> double
{
    // We're doing linear extrapolation in case one has to go outside
    // the interpolation range.
    if (x <= knots(0)) {
        return values(0) + A(0) * (x - knots(0));
    }
    if (x >= knots(knots.size() - 1)) {
        return values(values.size() - 1)
               + A(A.size() - 1) * (x - knots(knots.size() - 1));
    }
    const int idx { equal_spacing ? lookupIdx(x) : binaryFindIdx(knots, x) };
    const double Dx { x - knots(idx) };
    const double Dx2 { Dx * Dx };
    const double Dx3 { Dx2 * Dx };
    return values(idx) + A(idx) * Dx + B(idx) * Dx2 + C(idx) * Dx3;
}

[[nodiscard]] auto CubicSpline::eval(const Eigen::ArrayXd& x) const
  -> Eigen::ArrayXd
{
    Eigen::ArrayXd result(x.size());
    for (int i {}; i < static_cast<int>(x.size()); ++i) {
        result(i) = eval(x(i));
    }
    return result;
}

} // namespace tango
