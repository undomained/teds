// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "extract_spectra.h"
#include "../cubic_spline.h"
#include <numeric>

namespace tango {

//ExtractSpectra::ExtractSpectra() {
//}

//ExtractSpectra::~ExtractSpectra() {
//}

std::string ExtractSpectra::getName() const {
    return std::string("ExtractSpectra");
}

bool ExtractSpectra::algoCheckInput(const CKD& ckd, L1& l1) {
//    bool sig_input = l1.observation_sig.size() > 0;
//    bool wl_input = (*l1.observation_wl).size() > 0;
//    if (sig_input && wl_input){
//        return true;
//    } else {
//        return false;
//    }
    return true;
}

void ExtractSpectra::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");

   // Fill observation data with zero
    l1.observation_sig.resize(ckd.n_act);
    l1.observation_std.resize(ckd.n_act);
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        l1.observation_sig[i_act].assign(ckd.n_detector_cols, 0.0);
        l1.observation_std[i_act].assign(ckd.n_detector_cols, 0.0);
    }

    for (size_t i_col {}; i_col < ckd.n_detector_cols; ++i_col){
        for (size_t i_act {}; i_act < ckd.n_act; ++i_act) {

            const double row_mid = ckd.swath.row_indices[i_act][i_col];
            double row_first = row_mid;
            double row_last = row_mid;

            if (i_act > 0){
                row_first = ckd.swath.row_indices[i_act-1][i_col];
            }
            if (i_act < (ckd.n_act-1)){
                row_last = ckd.swath.row_indices[i_act+1][i_col];
            }
            int row_first_i {static_cast<int>(row_first) };
            int row_last_i {static_cast<int>(row_last) };
            int row_mid_i {static_cast<int>(row_mid) };

            // Check boundaries
            if (row_first_i < 0){
                row_first_i = 0;
            }
            if (row_mid_i < 0){
                row_mid_i = 0;
            }
            if (row_last_i < 0){
                row_last_i = 0;
            }
            const int n_row = ckd.n_detector_rows;
            if (row_first_i >= n_row){
                row_first_i = n_row-1;
            }
            if (row_mid_i >= n_row){
                row_mid_i = n_row-1;
            }
            if (row_last_i >= n_row){
                row_last_i = n_row-1;
            }

            // Determine weight factors for rows below row_mid and above row_mid
            int n_rows_below = row_mid_i - row_first_i;
            int n_rows_above = row_last_i - row_mid_i;
            double weight_factor_below = 0.;
            if (n_rows_below > 0){
                weight_factor_below = 1. / (n_rows_below) ;
            }
            double weight_factor_above = 0.;
            if (n_rows_above > 0){
                weight_factor_above = 1. / (n_rows_above) ;
            }

            // Add signal of rows below row_mid to signal according to weight
            double total_weight = 0.;
            double weight = 0.0;
            double n = 0.;
            for (int i_row {row_first_i}; i_row <= row_mid_i; ++i_row){
                total_weight += weight;
                l1.observation_sig[i_act][i_col] += weight * l1.image[i_row * ckd.n_detector_cols + i_col];
                n = {l1.stdev[i_row * ckd.n_detector_cols + i_col]*weight};
                l1.observation_std[i_act][i_col] += n*n;
                weight += weight_factor_below;
            }
            // Add signal of rows above row_mid to signal according to weight
            weight = 1.0;
            for (int i_row {row_mid_i}; i_row <= row_last_i; ++i_row){
                total_weight += weight;
                l1.observation_sig[i_act][i_col] += weight * l1.image[i_row * ckd.n_detector_cols + i_col];
                n = {l1.stdev[i_row * ckd.n_detector_cols + i_col]*weight};
                l1.observation_std[i_act][i_col] += n*n;
                weight -= weight_factor_above;
            }
            // Divide by total_weight
            l1.observation_sig[i_act][i_col] /= total_weight;
            l1.observation_std[i_act][i_col] = std::sqrt(l1.observation_std[i_act][i_col] / (total_weight * total_weight));
        }
    }

// TODO To be able to compare to input, should we not also go to original wavelength grid of spectra??????
//
}

} // namespace tango
