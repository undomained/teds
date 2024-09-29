// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dac.h"

namespace tango {

//DAC::DAC() {
//}

//DAC::~DAC() {
//}

std::string DAC::getName() const {
    if (getModelType() == "IM") {
        return std::string("ADC");
    } else {
        return std::string("DAC");
    }
}

bool DAC::algoCheckInput(const CKD& ckd, L1& l1) {
    // Check if image and ckd have the same dimensions
    if (l1.image.size() > 0) {
        return true;
    } else {
        spdlog::warn("No input image");
        return false;
    }

    return true;
}

//void DAC::unloadData() {
//    spdlog::info("DAC unload fct still to be filled in");
//}

void DAC::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");

    spdlog::warn("Conversion factor (bits or counts -> voltage) not implemented");
    float f_DAC = 1.0; 

    if (getModelType() == "L1B") {
        for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
            if (!l1.pixel_mask[i]) {
                l1.image[i] *= f_DAC;
                l1.stdev[i] *= f_DAC;
            }
        }
        l1.units = 'V';
    } else if (getModelType() == "IM") {
        for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
            if (!l1.pixel_mask[i]) {
                l1.image[i] /= f_DAC;
                l1.stdev[i] /= f_DAC; 
            }
        }
        l1.units = "counts";
    }
}

} // namespace tango
