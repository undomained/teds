// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "quaternion.h"

#include "algorithm.h"

#include <cmath>

namespace tango {

Quaternion::Quaternion(const std::array<double, dims::vec * dims::vec>& rot)
{
    if (const double trace { rot[0 * dims::vec + 0] + rot[1 * dims::vec + 1]
                             + rot[2 * dims::vec + 2] };
        trace >= 0.0) {
        const double s { 0.5 / std::sqrt(trace + 1.0) };
        v = { (rot[2 * dims::vec + 1] - rot[1 * dims::vec + 2]) * s,
              (rot[0 * dims::vec + 2] - rot[2 * dims::vec + 0]) * s,
              (rot[1 * dims::vec + 0] - rot[0 * dims::vec + 1]) * s };
        r = 0.25 / s;
    } else {
        int i {};
        if (rot[1 * dims::vec + 1] > rot[0 * dims::vec + 0]) {
            i = 1;
        }
        if (rot[1 * dims::vec + 1] > rot[i * dims::vec + i]) {
            i = 2;
        }
        const int j { (i + 1) % dims::vec };
        const int k { (j + 1) % dims::vec };
        double s { std::sqrt(rot[i * dims::vec + i] - rot[j * dims::vec + j]
                             - rot[k * dims::vec + k] + 1.0) };
        v[i] = 0.5 * s;
        s = 0.5 / s;
        r = (rot[k * dims::vec + j] - rot[j * dims::vec + k]) * s;
        v[j] = (rot[j * dims::vec + i] + rot[i * dims::vec + j]) * s;
        v[k] = (rot[k * dims::vec + i] + rot[i * dims::vec + k]) * s;
    }
}

auto Quaternion::normalize() -> void
{
    const double norm_inv { 1 / std::sqrt(r * r + dotProduct(v, v)) };
    for (double& val : v) {
        val *= norm_inv;
    }
    r *= norm_inv;
}

auto Quaternion::conjugate() -> void
{
    for (double& el : v) {
        el = -el;
    }
}

auto Quaternion::rotate(std::array<double, dims::vec>& a) const -> void
{
    const double vv { dotProduct(v, v) };
    const double av { dotProduct(a, v) };
    std::array<double, dims::vec> va {};
    crossProduct(v, a, va);
    for (int i {}; i < dims::vec; ++i) {
        a[i] = (r * r - vv) * a[i] + 2.0 * av * v[i] + 2.0 * r * va[i];
    }
}

auto Quaternion::multiply(const Quaternion& rhs) -> void
{
    double r_new { r * rhs.r - dotProduct(v, rhs.v) };
    std::array<double, dims::vec> cross {};
    crossProduct(v, rhs.v, cross);
    for (int i {}; i < dims::vec; ++i) {
        v[i] = r * rhs.v[i] + rhs.r * v[i] + cross[i];
    }
    r = std::move(r_new);
}

auto Quaternion::slerp(const Quaternion& q, const double t) -> void
{
    const double theta { std::acos(dotProduct(v, q.v) + r * q.r) };
    const double f1 { std::sin((1 - t) * theta) / sin(theta) };
    const double f2 { std::sin(t * theta) / sin(theta) };
    for (int i {}; i < dims::vec; ++i) {
        v[i] = f1 * v[i] + f2 * q.v[i];
    }
    r = f1 * r + f2 * q.r;
}

} // namespace tango
