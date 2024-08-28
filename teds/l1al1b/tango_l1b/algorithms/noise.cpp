// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "noise.h"
#include <cstdint>
#include <spdlog/spdlog.h>


#include <concepts>
#include <iostream>

namespace tango {

//Noise::Noise() {
//}

//Noise::~Noise() {
//}

std::string Noise::getName() const {
    return std::string("Noise");
}

bool Noise::algoCheckInput(const CKD& ckd, L1& l1) {
    // Check if image and ckd have the same dimensions
    if ((l1.image.size() != ckd.noise.g.size()) && (l1.image.size() != ckd.noise.n2.size())) {
        spdlog::warn("Image and CKD dimensions do not match, skipping");
        return false;
    } else {
        return true;
    }
}

//void Noise::unloadData() {
//    spdlog::info("Noise unload fct still to be filled in");
//}

void Noise::algoExecute(const CKD& ckd, L1& l1) {
    spdlog::warn("Noise algoExecute fct still to be filled in");
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (l1.pixel_mask[i]) {
            continue;
        }
        const double noise2 { ckd.noise.g[i] * l1.image[i] + ckd.noise.n2[i] };
        
        //if (noise2 <= 0.0) {
        //    l1.pixel_mask[i] = true;
        //} else if (binned_detector_image) {
        //    l1.stdev[i] +=
        //      std::sqrt(noise2 / (l1.nr_coadditions * binning_table.binSize(i)));
        //} else {
        //    l1.stdev[i] = std::sqrt(noise2 / l1.nr_coadditions);
        //}
    }
}

} // namespace tango
