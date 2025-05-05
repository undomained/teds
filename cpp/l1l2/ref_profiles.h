// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class with information about reference profiles used in retrieval

#pragma once

#include <map>
#include <string>

namespace tango {

struct RefProfiles
{
    std::map<std::string, Eigen::ArrayXd> gases {};
    std::map<std::string, double> gas_sums {};
    std::map<std::string, double> initial {};
};

} // namespace tango {
