// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "noise.h"
#include <random>

namespace tango {

//Noise::Noise() {
//}

//Noise::~Noise() {
//}

std::string Noise::getName() const {
    return std::string("Noise");
}

bool Noise::algoCheckInput(const CKD& ckd, L1& l1) {
    // Fill noise2 if needed
    if ((l1.stdev.size()) != l1.image.size()) {
        l1.stdev.resize(l1.image.size());
    }
    if ((l1.noise2.size()) != l1.image.size()) {
        l1.noise2.resize(l1.image.size());
    }
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
    spdlog::warn("Noise, binning not yet taken into account");
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (l1.pixel_mask[i]) {
            continue;
        }
        const double noise2 { ckd.noise.g[i] * l1.image[i] + ckd.noise.n2[i] };

        if (getModelType()=="L1B"){
            
            // Not sure if this is the right order of calculations
            if (noise2 <= 0.0) {
                l1.pixel_mask[i] = true;
            } else {
                l1.noise2[i] = noise2;
            }

        } else if (getModelType()=="IM"){
            const int seed = 100;
            static std::mt19937 gen { static_cast<unsigned long>(seed) };
            std::normal_distribution<> d { 0.0, std::sqrt(noise2) };
            l1.image[i] += d(gen);
        }
        // another option for binned images has to be added here 
    }
}

} // namespace tango
