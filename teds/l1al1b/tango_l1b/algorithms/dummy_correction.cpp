// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dummy_correction.h"


namespace tango {

//DummyCorrection::DummyCorrection() {
//}

//DummyCorrection::~DummyCorrection() {
//}

std::string DummyCorrection::getName() const {
    return std::string("DummyCorrection");
}

bool DummyCorrection::algoCheckInput(const CKD& ckd, L1& l1) {
    spdlog::warn("DummyCorrection algoCheckInput fct still to be filled in");

    return true;
}

//void DummyCorrection::unloadData() {
//    spdlog::info("DummyCorrection unload fct still to be filled in");
//}

void DummyCorrection::algoExecute(L1& l1, const Dataset& input_data) {

    spdlog::info("DummyCorrection algoExecute fct still to be filled in");
}

} // namespace tango
