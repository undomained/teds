// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Quaternion algebra

#pragma once

#include "constants.h"

#include <array>
#include <vector>

namespace tango {

class Quaternion
{
private:
    // Vector part
    std::array<double, dims::vec> v {};
    // Real part
    double r {};

public:
    Quaternion() = default;
    Quaternion(const double x, const double y, const double z, const double r)
      : v { x, y, z }, r { r }
    {}
    // Initialize from a rotation matrix using the algorithm of
    // "Quaternion Calculus and Fast Animation", Ken Shoemake, 1987 to
    // convert the rotation matrix to a quaternion. Also, see the
    // implementation in Eigen: src/Geometry/Quaternion.h .
    Quaternion(const std::array<double, dims::vec * dims::vec>& rot);

    auto operator[](const int i) const -> const double& { return v[i]; }
    auto operator[](const int i) -> double& { return v[i]; }
    auto real() const -> const double& { return r; }
    auto real() -> double& { return r; }

    auto normalize() -> void;
    // Reverse the vector component
    auto conjugate() -> void;
    // Apply Q V Q^-1 to vector V where Q is this quaternion. The
    // formula is Q V Q^-1 = (r^2 - v^2)V + 2(Vv)v + 2r(v x V).
    auto rotate(std::array<double, dims::vec>& a) const -> void;
    // Multiply this quaternion with rhs from the right and store the
    // result in self:
    //   (self.r * rhs.r - self.v * rhs.v,
    //    self.r * rhs.v + rhs.r * self.v + self.v x rhs.v)
    auto multiply(const Quaternion& rhs) -> void;
    // Spherical linear interpolation (Ken Shoemake "Animating
    // rotation with quaternion curves", SIGGRAPH Computer Graphics,
    // 19:245 (1985).
    auto slerp(const Quaternion& q, const double t) -> void;
    ~Quaternion() = default;
};

} // namespace tango
