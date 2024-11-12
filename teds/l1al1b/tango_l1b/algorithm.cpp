// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "algorithm.h"

#include "constants.h"

namespace tango {

[[nodiscard]] auto binaryFindIdx(const std::vector<double>& list,
                                 const double x) -> int
{
    int i_begin {};
    int i_end { static_cast<int>(list.size() - 1) };
    int i_mid {};
    if (x <= list.front()) {
        return i_begin;
    }
    if (x >= list.back()) {
        return i_end;
    }
    while (true) {
        i_mid = (i_begin + i_end) / 2;
        if (x < list[i_mid]) {
            i_end = i_mid - 1;
        } else if (x < list[i_mid + 1]) {
            return i_mid;
        } else {
            i_begin = i_mid + 1;
        }
    }
}

} // namespace tango
