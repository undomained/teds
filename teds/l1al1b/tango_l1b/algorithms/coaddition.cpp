// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "coaddition.h"

namespace tango {

//Coaddition::Coaddition() {
//}

//Coaddition::~Coaddition() {
//}

std::string Coaddition::getName() const {
    return std::string("Coaddition");
}

bool Coaddition::algoCheckInput(const CKD& ckd, L1& l1) {
    if (l1.nr_coadditions == 0) {
        spdlog::warn("nr of coadditions = 0, skipping coaddition step");
        return false;
    } else {
        return true;
    }
    
}

//void Coaddition::unloadData() {
//    spdlog::info("Coaddition unload fct still to be filled in");
//}

void Coaddition::algoExecute(const CKD& ckd, L1& l1) {
    if (getModelType() == "L1B"){
        spdlog::warn("Coaddition correction for L1B not implemented yet");
    }
    
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            if (getModelType() == "L1B"){
                l1.stdev[i] *= std::sqrt(l1.nr_coadditions);
            } else if (getModelType() == "IM"){
                l1.image[i] *= l1.nr_coadditions;
            }
        }
    }
}

} // namespace tango
