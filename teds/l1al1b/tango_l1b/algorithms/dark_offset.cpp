// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dark_offset.h"

namespace tango {

//DarkOffset::DarkOffset() {
//}

//DarkOffset::~DarkOffset() {
//}

std::string DarkOffset::getName() const {
    return std::string("DarkOffset");
}

bool DarkOffset::algoCheckInput(const CKD& ckd, L1& l1) {
    // Check if image and ckd have the same dimensions
    if (l1.image.size() == ckd.dark.offset.size()) {
        return true;
    } else {
        spdlog::warn("Image and CKD dimensions do not match, skipping");
        return false;
    }
}

//void DarkOffset::unloadData() {
//    spdlog::info("DarkOffset unload fct still to be filled in");
//}

void DarkOffset::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            if (getModelType() == "L1B"){
                l1.image[i] -= ckd.dark.offset[i];
            } else if (getModelType() == "IM"){
                l1.image[i] += ckd.dark.offset[i];
            }
        }
    }
}

} // namespace tango
