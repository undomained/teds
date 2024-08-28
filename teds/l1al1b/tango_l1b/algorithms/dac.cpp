// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dac.h"

namespace tango {

//DAC::DAC() {
//}

//DAC::~DAC() {
//}

std::string DAC::getName() const {
    return std::string("DAC");
}

bool DAC::algoCheckInput(const CKD& ckd, L1& l1) {
    spdlog::warn("DAC algoCheckInput fct still to be filled in");

    return true;
}

//void DAC::unloadData() {
//    spdlog::info("DAC unload fct still to be filled in");
//}

void DAC::algoExecute(const CKD& ckd, L1& l1) {
    spdlog::warn("Conversion factor (bits->voltage) not implemented");
    float f_DAC = 1.0;
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            l1.image[i] *= f_DAC;
        }
    }
}

} // namespace tango
