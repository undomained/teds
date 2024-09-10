// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "radiometric.h"

namespace tango {

//Radiometric::Radiometric() {
//}

//Radiometric::~Radiometric() {
//}

std::string Radiometric::getName() const {
    return std::string("Radiometric");
}

bool Radiometric::algoCheckInput(const CKD& ckd, L1& l1) {
    // Check if image and ckd have the same dimensions
    if (l1.image.size() != ckd.rad.rad.size()) {
        spdlog::warn("Image and CKD dimensions do not match, skipping");
        return false;
    } else {
        return true;
    }
}

//void Radiometric::unloadData() {
//    spdlog::info("Radiometric unload fct still to be filled in");
//}

void Radiometric::algoExecute(const CKD& ckd, L1& l1) {
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            if (getModelType() == "L1B"){
                l1.image[i] *= (ckd.rad.rad[i] / l1.exposure_time);
                l1.units = "ph nm-1 s-1 sr-1 m-2"; // check if this is correct
            } else if (getModelType() == "IM"){
                l1.image[i] /= (ckd.rad.rad[i] / l1.exposure_time);
                l1.units = "counts";
            }
        }
    }
}

} // namespace tango
