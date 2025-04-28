// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "b_spline.h"

#include "algorithm.h"

namespace tango {

BSpline::BSpline(const int order, const Eigen::ArrayXd& knots) : order { order }
{
    // Number of knots to skip for better smoothness
    const int i_start { order / 2 + 1 };
    this->knots.resize(knots.size() + order + 1);
    // Copy the B-spline knots and pad them from both ends by (order)
    // elements.
    for (int i {}; i < static_cast<int>(knots.size()) - order - 1; ++i) {
        this->knots(i + order + 1) = knots(i + i_start);
    }
    for (int i {}; i < order + 1; ++i) {
        this->knots(i) = knots(0);
        this->knots(this->knots.size() - 1 - i) = knots(knots.size() - 1);
    }
    control_points_tmp.resize(this->knots.size(), 0.0);
}

[[nodiscard]] auto BSpline::nStates() const -> int
{
    return static_cast<int>(knots.size()) - order - 1;
}

auto BSpline::findInterval(const double x) const -> int
{
    // The lowest and highest intervals are 'order' and
    // 'knots.size()-order-2', meaning knots added for padding are
    // excluded.
    return std::min(std::max(order, binaryFindIdx(knots, x)),
                    static_cast<int>(knots.size()) - order - 2);
}

auto BSpline::deBoor(const std::vector<double>& control_points,
                     const double x) const -> double
{
    const int p { findInterval(x) };
    constexpr int max_b_spline_order { 20 };
    std::array<double, max_b_spline_order> d_cur;
    std::copy(control_points.cbegin() + p - order,
              control_points.cbegin() + p + 1,
              d_cur.begin());
    for (int r { 1 }; r <= order; ++r) {
        for (int j { order }; j >= r; --j) {
            const double alpha { (x - knots[j + p - order])
                                 / (knots[j + 1 + p - r]
                                    - knots[j + p - order]) };
            d_cur[j] = (1 - alpha) * d_cur[j - 1] + alpha * d_cur[j];
        }
    }
    return d_cur[order];
}

auto BSpline::genBasis(const Eigen::ArrayXd& x_data)
  -> Eigen::SparseMatrix<double>
{
    std::vector<Eigen::Triplet<double>> triplets {};
    triplets.reserve(x_data.size() * (order + 1));
    for (int i_point {}; i_point < static_cast<int>(x_data.size()); ++i_point) {
        const int i_x { findInterval(x_data[i_point]) };
        for (int k {}; k <= order; ++k) {
            control_points_tmp[i_x - k] = 1.0;
            triplets.push_back({ i_point,
                                 i_x - k,
                                 deBoor(control_points_tmp, x_data[i_point]) });
            control_points_tmp[i_x - k] = 0.0;
        }
    }
    Eigen::SparseMatrix<double> B(x_data.size(), nStates());
    B.setFromTriplets(triplets.cbegin(), triplets.cend());
    return B;
}

} // namespace tango
