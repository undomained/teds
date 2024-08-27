// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
/// Build class for the algorithms.
#include "build_algo.h"

#include "prnu_correction.h"
#include "dummy_correction.h"

namespace tango {

BuildAlgo::BuildAlgo() {
    algo_map["DummyCorrection"] = []() -> BaseAlgo* { return new DummyCorrection(); };
    algo_map["PRNUCorrection"] = []() -> BaseAlgo* { return new PRNUCorrection(); };
}

BuildAlgo::~BuildAlgo() = default;

BaseAlgo* BuildAlgo::CreateAlgo(const std::string& name) {
    auto it = algo_map.find(name);
    if (it != algo_map.end()) {
        return it->second();  // Call the lambda function to create the corresponding algorithm
    } else {
        return nullptr;  // Algorithm not found
    }
}

} // namespace tango


