// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// General purpose math routines

#pragma once

#include "constants.h"

#include <array>
#include <vector>

namespace tango {

// Do a binary search to locate an index n such that x falls in the
// range list[n]...list[n+1].
[[nodiscard]] auto binaryFindIdx(const std::vector<double>& list,
                                 const double x) -> int;

[[nodiscard]] auto dotProduct(const std::array<double, dims::vec>& a,
                              const std::array<double, dims::vec>& b) -> double;

template <typename T>
auto crossProduct(const T a,
                  const T b,
                  std::array<double, dims::vec>& c) -> void
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

} // namespace tango
