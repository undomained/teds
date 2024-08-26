// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
#include "base_algo.h"

#include <spdlog/spdlog.h>

namespace tango {

BaseAlgo::BaseAlgo()
{
}

BaseAlgo::~BaseAlgo()
{
}

void BaseAlgo::algoCheckInput(const CKD& ckd, L1& l1) {
    spdlog::info("BaseAlgo algoCheckInput fct still to be filled in");
}

std::string BaseAlgo::getName() const {
    return std::string("BaseAlgo");
}

//void BaseAlgo::unloadData() {
//    spdlog::info("BaseAlgo unload fct still to be filled in");
//}

void BaseAlgo::algoExecute(const CKD& ckd, const bool enabled, L1& l1)
{
    spdlog::info("BaseAlgo algoExecute fct still to be filled in");
}

} // namespace tango

