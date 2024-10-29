// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "ISRF.h"
#include "../cubic_spline.h"
#include <numeric>
#include <iostream>
#include <cmath> // for exp, abs and pow

namespace tango {

//ISRF::ISRF() {
//}

//ISRF::~ISRF() {
//}

std::string ISRF::getName() const {
    return std::string("ISRF");
}

bool ISRF::algoCheckInput(const CKD& ckd, L1& l1) {
    // Check if isrf wavelengths match with input wavelengths. If so, wavelength
    // dependend ISRF can be applied
    int n_input_samples = ckd.n_act * ckd.n_lbl;
    int n_x0 = ckd.isrf.x0.size();
    int w = ckd.isrf.w.size();
    int n = ckd.isrf.n.size();

    auto& wl_input = (*l1.observation_wl).front(); 
    int n_wl_input = wl_input.size(); // number of input wavelengths

    if (n_x0 != n_input_samples || w != n_input_samples || n != n_input_samples ){
        // CKD size not correct
        return false;
    } else if (ckd.n_lbl != n_wl_input){
        spdlog::error("Number of line-by-line samples given in CKD does not match input spectrum");
        return false;
    } else {
        return true;
    }
}

//void ISRF::unloadData() {
//    spdlog::info("ISRF unload fct still to be filled in");
//}

void ISRF::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");

    // ISRF with highest width is below 0.001 at ±0.6193nm
    // We compute ISRF to ±0.8nm
    float wl_differ = 0.0; // wavelength increments from input spectra
    auto& wl_input = (*l1.observation_wl).front(); 
    std::vector<float> inc_input(ckd.n_lbl - 1);
    for (int i_wl {}; i_wl < ckd.n_lbl - 1; ++i_wl) {
            inc_input[i_wl] = wl_input[i_wl + 1]- wl_input[i_wl];
        }
    float mean_inc_input =
            std::accumulate(inc_input.begin(), inc_input.end(), 0.0)
            / (inc_input.size());
    
    // number of indices to compute ISRF
    int max_ix_wing = (0.8 / mean_inc_input) + 1; 
    std::vector<float> dwl(2 * max_ix_wing + 1);
    int n_dwl = dwl.size();
    for (int i {}; i < dwl.size(); i++) {
        dwl[i] = (i - max_ix_wing) * mean_inc_input; // wavelength difference. Used as x in supergaussian
    }

    for (int i_act {}; i_act < ckd.n_act; i_act++) {
        std::vector<double> signal_conv(ckd.n_lbl, 0.0);
        for (int i_wl {}; i_wl < ckd.n_lbl; i_wl++) {
            float x0 = ckd.isrf.x0[i_act * ckd.n_lbl + i_wl];
            float w = ckd.isrf.w[i_act * ckd.n_lbl + i_wl];
            float n = ckd.isrf.n[i_act * ckd.n_lbl + i_wl];

            // left bound of dwl
            int i_dwl_0 = std::max(max_ix_wing - i_wl, 0); 
            int i_dwl_1 = std::min(n_dwl - 1 , ckd.n_lbl - i_wl + max_ix_wing - 1);
            
            // Calculate in-bound part of ISRF
            std::vector<float> this_dwl(dwl.begin() + i_dwl_0, dwl.begin() + i_dwl_1 + 1);
            std::vector<float> this_isrf(i_dwl_1 - i_dwl_0 + 1, 0);
            for (int i {}; i < this_dwl.size(); i++) {
                float x = this_dwl[i];
                this_isrf[i] = std::exp(-std::pow(std::abs(x-x0) / w, 2 * n));
            }

            // Normalization factor 
            double norm_inv = 1 / (std::accumulate(this_isrf.begin(), this_isrf.end(), 0.0)); 

            // Left and right bounds of input spectrum,
            int i_input_0 = std::max(i_wl - max_ix_wing, 0);
            int i_input_1 = std::min(i_wl + max_ix_wing, ckd.n_lbl - 1);

            // Calculate in-bound part of input spectrum
            std::vector<double> this_signal(
                l1.observation_sig[i_act].begin() + i_input_0,
                l1.observation_sig[i_act].begin() + i_input_1 + 1);


            // Carry out convolution
            for (int k {}; k < this_isrf.size(); ++k) {
                signal_conv[i_wl] += norm_inv * this_isrf[k] * this_signal[k];
            }
        }

        l1.observation_sig[i_act] = signal_conv;
        if (i_act ==0){
            for (int i_wl {}; i_wl < ckd.n_lbl; i_wl++) {
            spdlog::info("{}", l1.observation_sig[i_act][i_wl]);
            }
        }
    }

    /* OLD VERSION
    int center_ix = (ckd.n_isrf_samples - 1) / 2;
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        std::vector<double> signal_conv(ckd.n_lbl, 0.0);
        for (int i_wl {}; i_wl < ckd.n_lbl; ++i_wl) {

            // Left and right bounds of isrf samples, convolution kernel
            int i_isrf_0 = std::max(center_ix - i_wl, 0);
            int i_isrf_1 = std::min(ckd.n_lbl + center_ix - i_wl - 1,
                                    ckd.n_isrf_samples - 1);
            // Left and right bounds of wavelengths,
            int i_input_0 = std::max(i_wl - center_ix, 0);
            int i_input_1 = std::min(i_wl + center_ix, ckd.n_lbl - 1);

            std::vector<double> this_isrf(
                ckd.rad.isrf[i_act][i_wl].begin() + i_isrf_0,
                ckd.rad.isrf[i_act][i_wl].begin() + i_isrf_1);
            double norm_inv =
                1
                / (std::accumulate(this_isrf.begin(), this_isrf.end(), 0.0));

            std::vector<double> this_signal(
                l1.observation_sig[i_act].begin() + i_input_0,
                l1.observation_sig[i_act].begin() + i_input_1);
            for (int k {}; k < this_isrf.size(); ++k) {
                signal_conv[i_wl] +=
                    norm_inv * this_isrf[k] * this_signal[k];
            }
        }
        
        l1.observation_sig[i_act] = signal_conv;
    }
    */
}

} // namespace tango
