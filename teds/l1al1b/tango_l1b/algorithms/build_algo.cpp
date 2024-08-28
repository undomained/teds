// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
/// Build class for the algorithms.
#include "build_algo.h"

// Include all correction algorithm header files below
#include "dummy_correction.h"
#include "prnu.h"
#include "dark_offset.h"
#include "dark_current.h"
#include "noise.h"
#include "straylight.h"
#include "radiometric.h"
#include "dac.h"
#include "coaddition.h"

namespace tango {


// Add algorithm names to dictionary in constructor
// Left: algo name defined in proctable
// Right: algo name defined in cpp file. 
BuildAlgo::BuildAlgo() {
    algo_map["DummyCorrection"] = []() -> BaseAlgo* { return new DummyCorrection(); };
    algo_map["PRNU"] = []() -> BaseAlgo* { return new PRNU(); };
    algo_map["Dark_Offset"] = []() -> BaseAlgo* { return new DarkOffset(); };
    algo_map["Dark_Current"] = []() -> BaseAlgo* { return new DarkCurrent(); };
    algo_map["Noise"] = []() -> BaseAlgo* { return new Noise(); };
    algo_map["Straylight"] = []() -> BaseAlgo* { return new Straylight(); };
    algo_map["Radiometric"] = []() -> BaseAlgo* { return new Radiometric(); };
    algo_map["DAC"] = []() -> BaseAlgo* { return new DAC(); };
    algo_map["Coaddition"] = []() -> BaseAlgo* { return new Coaddition(); };
}

BuildAlgo::~BuildAlgo() = default;

BaseAlgo* BuildAlgo::CreateAlgo(const std::string& name) {
    auto it = algo_map.find(name);
    if (it != algo_map.end()) {
        return it->second();  // Call the lambda function to create the corresponding algorithm
    } else {
        spdlog::warn("Algorithm {} not found, check BuildAlgo and proctable", name);
        return nullptr;  // Algorithm not found
    }
}

} // namespace tango


