// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "draw_on_detector.h"
#include "../cubic_spline.h"
#include <numeric>

namespace tango {

//DrawOnDetector::DrawOnDetector() {
//}

//DrawOnDetector::~DrawOnDetector() {
//}

std::string DrawOnDetector::getName() const {
    return std::string("Draw On Detector");
}

bool DrawOnDetector::algoCheckInput(const CKD& ckd, L1& l1) {
    bool sig_input = l1.observation_sig.size() > 0;
    bool wl_input = (*l1.observation_wl).size() > 0;
    if (sig_input && wl_input){
        return true;
    } else {
        return false;
    }
}

//void DrawOnDetector::unloadData() {
//    spdlog::info("DrawOnDetector unload fct still to be filled in");
//}

void DrawOnDetector::algoExecute(const CKD& ckd, L1& l1) {
    int n_pixels = ckd.n_detector_rows * ckd.n_detector_cols;
    l1.image.assign(n_pixels, 0.0);

    // Create a range for columns and rows
    std::vector<double> cols(ckd.n_detector_cols);
    for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col) {
        cols[i_col] = i_col;
    }

    // First convert wavelengths to columns
    std::vector<std::vector<double>> lbl_in_cols(
        ckd.n_act, std::vector<double>(ckd.n_detector_cols, 0.0)); 
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        
        std::vector<double> wl = ckd.wave.wavelength[i_act]; // wavelength per column (ckd)
        std::vector<double> lbl = (l1.observation_sig[i_act]); // spectrum signal
        std::vector<double> observation_wls =(*l1.observation_wl)[i_act]; // spectrum wavelength

        const CubicSpline spl_wl_to_sig { observation_wls, lbl}; // spline of signal
        // Interpolate signal to wavelength of column
        for (int i_col :cols){
            lbl_in_cols[i_act][i_col] = spl_wl_to_sig.eval(wl[i_col]);
        }
    }

    // Now per column interpolate act_pos to rows
    for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col){
        std::vector<double> lbl_this_col(ckd.n_act, 0.0);
        std::vector<double> row_this_col(ckd.n_act, 0.0);
        for (int i_act {}; i_act < ckd.n_act; ++i_act){
            lbl_this_col[i_act] = lbl_in_cols[i_act][i_col];
            row_this_col[i_act] = ckd.swath.row_indices[i_act][i_col];
        }

        const CubicSpline spl_lbl_vs_row {row_this_col, lbl_this_col};
        // Interpolate lbl spectrum to integer (detector) rows
        for (int i_row{}; i_row< ckd.n_detector_rows; ++i_row){
            l1.image[i_row * ckd.n_detector_cols + i_col] =
                spl_lbl_vs_row.eval(i_row);
        }
    }
    std::vector<double> image_zeros(l1.image.size());
    l1.stdev = image_zeros;
    *l1.wavelength = ckd.wave.wave_map; // Set wavelength
    l1.units = "counts";
}

} // namespace tango
