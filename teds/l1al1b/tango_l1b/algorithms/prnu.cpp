// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "prnu.h"

namespace tango {

//PRNU::PRNU() {
//}

//PRNU::~PRNU() {
//}

std::string PRNU::getName() const {
    return std::string("PRNU");
}

bool PRNU::algoCheckInput(const CKD& ckd, L1& l1)
{
    // Check if image and ckd have the same dimensions
    if (l1.image.size() == ckd.prnu.prnu.size()) {
        return true;
    } else {
        spdlog::warn("Image and CKD dimensions do not match, skipping");
        return false;
    }
}

//void PRNU::unloadData() {
////    m_prnuBinned = NULL;
////    BaseAlgorithm::unloadData();
//    spdlog::info("PRNU unload fct still to be filled in");
//}

void PRNU::algoExecute(L1& l1, const Dataset& input_data) {
    CKD const& ckd = input_data.get_container<CKD>("ckd");
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            if (getModelType() == "L1B"){
                l1.image[i] /= ckd.prnu.prnu[i];
                l1.stdev[i] /= ckd.prnu.prnu[i];
            } else if (getModelType() == "IM"){
                l1.image[i] *= ckd.prnu.prnu[i];
//                l1.stdev[i] *= ckd.prnu.prnu[i];
            }
        }
    }
}

} // namespace tango
