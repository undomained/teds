// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
#include "base_algo.h"

namespace tango {

BaseAlgo::BaseAlgo()
{
}

BaseAlgo::~BaseAlgo()
{
}

bool BaseAlgo::algoCheckInput(L1& l1, const Dataset& input_data){
    spdlog::info("BaseAlgo algoCheckInput fct still to be filled in");
    return true;
}

std::string BaseAlgo::getName() const {
    return std::string("BaseAlgo");
}

//void BaseAlgo::unloadData() {
//    spdlog::info("BaseAlgo unload fct still to be filled in");
//}

void BaseAlgo::algoExecute(L1& l1, const Dataset& input_data)
{
    spdlog::info("BaseAlgo algoExecute fct still to be filled in");
}

} // namespace tango

