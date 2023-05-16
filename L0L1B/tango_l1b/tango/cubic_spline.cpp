// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "cubic_spline.h"

#include "linalg.h"

#include <cmath>

// Threshold for determining if the knots are equally spaced
constexpr double equal_spacing_tol {
    1e1 * std::numeric_limits<double>::epsilon() };

CubicSpline::CubicSpline(const std::vector<double>& x_values,
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
    std::vector<double> delta_x(knots.size() - 1);
    for (int i {}; i < static_cast<int>(delta_x.size()); ++i) {
        delta_x[i] = knots[i + 1] - knots[i];
    }
    std::vector<double> delta_y(knots.size() - 1);
    for (int i {}; i < static_cast<int>(delta_y.size()); ++i) {
        delta_y[i] = values[i + 1] - values[i];
    }
    std::vector<double> diag(knots.size());
    diag.front() = diag.back() = 1.0;
    for (int i { 1 }; i < static_cast<int>(knots.size()) - 1; ++i) {
        diag[i] = 2.0 / 3.0 * (delta_x[i - 1] + delta_x[i]);
    }
    std::vector<double> diag_upper(diag.size() - 1);
    std::vector<double> diag_lower(diag.size() - 1);
    for (int i {}; i < static_cast<int>(diag_upper.size()) - 1; ++i) {
        diag_upper[i + 1] = delta_x[i + 1] / 3;
        diag_lower[i] = delta_x[i] / 3;
    }
    diag_upper.front() = 0.0;
    diag_lower.back() = 0.0;
    const int n { static_cast<int>(knots.size()) };
    std::vector<double> du2(knots.size() - 2);
    std::vector<int> ipiv(knots.size());
    int info {};
    dgttrf_(&n,
            diag_lower.data(),
            diag.data(),
            diag_upper.data(),
            du2.data(),
            ipiv.data(),
            &info);
    const char trans { 'n' };
    const int nrhs { 1 };
    B.resize(knots.size());
    B.front() = B.back() = 0.0;
    for (int i { 1 }; i < static_cast<int>(B.size()) - 1; ++i) {
        B[i] = delta_y[i] / delta_x[i] - delta_y[i - 1] / delta_x[i - 1];
    }
    dgttrs_(&trans,
            &n,
            &nrhs,
            diag_lower.data(),
            diag.data(),
            diag_upper.data(),
            du2.data(),
            ipiv.data(),
            B.data(),
            &n,
            &info);
    C.resize(B.size() - 1);
    for (int i {}; i < static_cast<int>(C.size()); ++i) {
        C[i] = (B[i + 1] - B[i]) / (3.0 * delta_x[i]);
    }
    A.resize(B.size() - 1);
    for (int i {}; i < static_cast<int>(A.size()); ++i) {
        A[i] =
          delta_y[i] / delta_x[i] - delta_x[i] * (2.0 * B[i] + B[i + 1]) / 3.0;
    }
}

auto CubicSpline::lookupIdx(const double x) const -> int
{
    return static_cast<int>((x - knots.front()) / range * (knots.size() - 1));
}

auto CubicSpline::binaryFindIdx(const double x) const -> int
{
    int i_begin {};
    int i_end { static_cast<int>(knots.size() - 1) };
    int i_mid {};
    while (true) {
        i_mid = (i_begin + i_end) / 2;
        if (x < knots[i_mid]) {
            i_end = i_mid - 1;
        } else if (x < knots[i_mid + 1]) {
            return i_mid;
        } else {
            i_begin = i_mid + 1;
        }
    }
}

auto CubicSpline::eval(const double x) const -> double
{
    // We're doing linear extrapolation in case one has to go outside
    // the interpolation range.
    if (x <= knots.front()) {
        return values.front() + A.front() * (x - knots.front());
    }
    if (x >= knots.back()) {
        return values.back() + A.back() * (knots.back() - x);
    }
    const int idx { equal_spacing ? lookupIdx(x) : binaryFindIdx(x) };
    const double Dx { x - knots[idx] };
    const double Dx2 { Dx * Dx };
    const double Dx3 { Dx2 * Dx };
    return values[idx] + A[idx] * Dx + B[idx] * Dx2 + C[idx] * Dx3;
}
