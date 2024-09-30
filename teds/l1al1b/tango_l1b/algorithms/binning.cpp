// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "binning.h"

namespace tango {

//Binning::Binning() {
//}

//Binning::~Binning() {
//}

std::string Binning::getName() const {
    return std::string("Binning");
}

bool Binning::algoCheckInput(const CKD& ckd, L1& l1) {
    if (l1.nr_coadditions == 0) {
        spdlog::warn("nr of coadditions = 0, skipping coaddition step");
        return false;
    } else {
        return true;
    }
    
}

//void Binning::unloadData() {
//    spdlog::info("Binning unload fct still to be filled in");
//}

void Binning::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");
    BinningTable const& binning = input_data.get_container<BinningTable>("binning");

    binning.bin(l1.image);   
    // For some reason the bin fct devides by binSize. Need to multiply by binsize here.
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        l1.image[i] = l1.image[i] * binning.binSize(i);
    }
}

} // namespace tango
