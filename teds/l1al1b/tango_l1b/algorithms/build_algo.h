// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
/// Build class for the algorithms.

#pragma once

#include <string>
#include "base_algo.h"
#include <unordered_map>
#include <functional>


namespace tango {

//class PRNUCorrection;
//class DummyCorrection;

class BuildAlgo {
public:
    BuildAlgo();
    ~BuildAlgo();

    BaseAlgo* CreateAlgo(const std::string& name);

private:
    std::unordered_map<std::string, std::function<BaseAlgo*()>> algo_map;
};

} // namespace tango


