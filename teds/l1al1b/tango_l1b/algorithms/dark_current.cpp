// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dark_current.h"

namespace tango {

//DarkCurrent::DarkCurrent() {
//}

//DarkCurrent::~DarkCurrent() {
//}

std::string DarkCurrent::getName() const {
    return std::string("DarkCurrent");
}

bool DarkCurrent::algoCheckInput(L1& l1, const Dataset& input_data){
    // Check if image and ckd have the same dimensions
    CKD const& ckd = input_data.get_container<CKD>("ckd");
    if (l1.image.size() != ckd.dark.current.size()) {
        spdlog::warn("Dark Current: Image and CKD dimensions do not match, skipping");
        return false;
    } else if (l1.exposure_time <= 0) {
        spdlog::warn("Exposure time = 0, skipping");
        return false;
    } else {
        return true;
    }
}

//void DarkCurrent::unloadData() {
//    spdlog::info("DarkCurrent unload fct still to be filled in");
//}

void DarkCurrent::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            if (getModelType() == "L1B"){
                l1.image[i] -= ckd.dark.current[i] * l1.exposure_time;
            } else if (getModelType() == "IM"){
                l1.image[i] += ckd.dark.current[i] * l1.exposure_time;
            }
        }
    }
}

} // namespace tango
