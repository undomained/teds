// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// General purpose math routines

#pragma once

#include <vector>

namespace tango {

// Do a binary search to locate an index n such that x falls in the
// range list[n]...list[n+1].
[[nodiscard]] auto binaryFindIdx(const std::vector<double>& list,
                                 const double x) -> int;

} // namespace tango
