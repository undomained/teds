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

void Binning::binWavelength(L1& l1, BinningTable const& binning) {

    const auto n_binned_pixels = l1.image.size();

    const auto& wavelengths = *l1.wavelength;
    const auto n_detector_cols { wavelengths.front().size() };
    const auto n_detector_rows {(*l1.wavelength).size()};

    const auto n_binned_rows = static_cast <int>(round(n_binned_pixels/n_detector_cols));

    auto wavelength_binned = std::make_shared<std::vector<std::vector<double>>>(n_binned_rows, std::vector<double>(n_detector_cols));

    // Fill binned data with zeroes
    for (size_t i_bin {}; i_bin < n_binned_rows; ++i_bin) {
        for (size_t i {}; i < n_detector_cols; ++i) {
            (*wavelength_binned)[i_bin][i] = 0.0;
        }
    }

    for (size_t i_row {}; i_row < n_detector_rows; ++i_row) {
        for (size_t i_col {}; i_col < n_detector_cols; ++i_col) {
            int idx = i_row *n_detector_cols + i_col ;
            int i_row_bin = int(binning.binIndex(idx)/n_detector_cols);
            (*wavelength_binned)[i_row_bin][i_col] += wavelengths[i_row][i_col];
        }
    }

// todo THIS ONLY works in case there is no binning in col direction. Needs to be fixed.
    for (size_t i_col {}; i_col < n_detector_cols; ++i_col) {
        for (int i_row {}; i_row < static_cast<int>((*wavelength_binned).size()); ++i_row) {
            int idx = i_row*n_detector_cols + i_col;
            (*wavelength_binned)[i_row][i_col] /= binning.binSize(idx);
        }
    }

    l1.wavelength_binned = wavelength_binned;
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

    // wavelength also needs to be binned.
    binWavelength(l1, binning);
}

} // namespace tango
