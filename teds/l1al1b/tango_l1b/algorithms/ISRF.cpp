// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "ISRF.h"
#include "../cubic_spline.h"
#include <numeric>

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
    float wl_differ = 0.0;
    int n_wl_input = (*l1.wavelength).front().size();
    int n_wl_isrf = (ckd.n_lbl);
    int center_ix = (ckd.n_isrf_samples - 1) / 2;
    if (n_wl_input == n_wl_isrf) {
        // sizes match, now look for matching wavelengths
        if ((*l1.wavelength).front()[0]
            == ckd.rad.isrf_wl.front().front()[center_ix]) {
            // first wavelength matches as well, now look at increments
            for (int i_act {}; i_act < ckd.n_act; ++i_act) {
                std::vector<double> inc_isrf(n_wl_isrf - 1);
                std::vector<double> inc_input(n_wl_input - 1);
                for (int i_wl {}; i_wl < n_wl_input - 1; ++i_wl) {
                    inc_isrf[i_wl] = ckd.rad.isrf_wl[0][i_wl + 1][center_ix]
                                     - ckd.rad.isrf_wl[0][i_wl][center_ix];
                    inc_input[i_wl] = (*l1.wavelength)[i_act][i_wl + 1]
                                      - (*l1.wavelength)[i_act][i_wl];
                }
                double mean_inc_isrf =
                  std::accumulate(inc_isrf.begin(), inc_isrf.end(), 0.0)
                  / (inc_isrf.size());
                double mean_inc_input =
                  std::accumulate(inc_input.begin(), inc_input.end(), 0.0)
                  / (inc_input.size());
                wl_differ += std::abs(mean_inc_input - mean_inc_isrf);
            }
        } else {
            wl_differ = 1;
        }
    }

    if (wl_differ == 0) {
        return true;
    } else {
        return false;
    }
    
}

//void ISRF::unloadData() {
//    spdlog::info("ISRF unload fct still to be filled in");
//}

void ISRF::algoExecute(const CKD& ckd, L1& l1) {
    spdlog::info("Applying wavelength dependend ISRF");
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
        l1.observation_sig[i_act] = std::move(signal_conv);
    }
}

} // namespace tango
