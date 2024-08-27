// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dark_offset.h"
#include <cstdint>
#include <spdlog/spdlog.h>


#include <concepts>
#include <iostream>

namespace tango {

//DarkOffset::DarkOffset() {
//}

//DarkOffset::~DarkOffset() {
//}

std::string DarkOffset::getName() const {
    return std::string("DarkOffset");
}

void DarkOffset::algoCheckInput(const CKD& ckd, L1& l1) {
    spdlog::info("DarkOffset algoCheckInput fct still to be filled in");
}

//void DarkOffset::unloadData() {
//    spdlog::info("DarkOffset unload fct still to be filled in");
//}

void DarkOffset::algoExecute(const CKD& ckd, L1& l1) {
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            l1.image[i] -= ckd.dark.offset[i];
        }
    }
}

} // namespace tango
