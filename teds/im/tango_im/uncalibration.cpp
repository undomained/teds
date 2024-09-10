// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "uncalibration.h"

#include <algorithm>
#include <iostream>
#include <random>
#include <ranges>
#include <spdlog/spdlog.h>
#include <tango_l1b/binning_table.h>
#include <tango_l1b/ckd.h>
#include <tango_l1b/cubic_spline.h>
#include <tango_l1b/fourier.h>
#include <tango_l1b/l1.h>

namespace tango {

auto applyISRF(const CKD& ckd,
               const bool enabled,
               const double fwhm_gauss,
               L1& l1) -> void
{
    
    // If this process is disabled then linearly interpolate the
    // line-by-line spectra onto the CKD wavelength grids. We cannot
    // simply return like the other processes.
    if (!enabled) {
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            const LinearSpline spline { (*l1.wavelength)[i_act],
                                        l1.spectra[i_act].signal };
            l1.spectra[i_act].signal.resize(ckd.n_detector_cols);
            for (int i {}; i < ckd.n_detector_cols; ++i) {
                l1.spectra[i_act].signal[i] =
                  spline.eval(ckd.wave.wavelength[i_act][i]);
            }
        }
        return;
    }

    // Check if isrf wavelengths match with input wavelengths. If so, wavelength
    // dependend ISRF can be applied
    double wl_differ = 0.0;
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

    // If the wavelength dependend ISRF ckd matches the input scene, we apply
    // the wavelength dependend ISRF Otherwise a constant gaussian distribution
    // will be used
    if (!wl_differ) {
        spdlog::info("Applying wavelength dependend ISRF");
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
                  l1.spectra[i_act].signal.begin() + i_input_0,
                  l1.spectra[i_act].signal.begin() + i_input_1);
                for (int k {}; k < this_isrf.size(); ++k) {
                    signal_conv[i_wl] +=
                      norm_inv * this_isrf[k] * this_signal[k];
                }
            }
            l1.spectra[i_act].signal = std::move(signal_conv);
        }
    } else {
        spdlog::warn("Input scene and ISRF CKD have different wavelength "
                     "grids, using a constant ISRF");
        // Convolve the radiances with an ISRF from the line-by-line grid
        const double sigma { fwhm_gauss / (2.0 * sqrt(2.0 * log(2.0))) };
        const double sigma_inv { 1.0
                                 / (2.0 * sigma * sigma) }; // norm of exponent
        const double norm_inv { 1.0 / (sigma * std::sqrt(2 * M_PI)) }; // norm
        std::vector<double> inc_input(n_wl_input - 1);
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            std::vector<double> signal_conv(ckd.n_lbl,
                                            0.0); // convolution result
            // Create normal distribution array with same increments as lbl
            // spectrum
            for (int i_wl {}; i_wl < n_wl_input - 1; ++i_wl) {
                inc_input[i_wl] = (*l1.wavelength)[i_act][i_wl + 1]
                                  - (*l1.wavelength)[i_act][i_wl];
            }
            double inc =
              std::accumulate(inc_input.begin(), inc_input.end(), 0.0)
              / (inc_input.size());
            // set mu = 0, always center at current evaluated wavelength
            // set limits of the normal distribution at 3 times FWHM
            int n_wings =
              std::ceil(3.0 * fwhm_gauss / inc) - 1; // samples in wings
            int n_gauss = 2 * n_wings + 1;           // total samples in gauss

            // now create distribution
            std::vector<double> gauss(n_gauss);
            for (int i {}; i < n_gauss; ++i) {
                double x = (i - n_wings) * inc;
                // make sure to multiply with the stepsize because its not a
                // continuous function
                gauss[i] = norm_inv * std::exp(-sigma_inv * x * x) * inc;
            }
            // Carry out convolution
            for (int i_wl {}; i_wl < ckd.n_lbl; ++i_wl) {
                // set left and right bounds for input (lbl) and gauss to
                // prevent out-of-bounds
                int i_lbl_0 = std::max(i_wl - n_wings, 0);
                int i_lbl_1 = std::min(i_wl + n_wings, ckd.n_lbl - 1);
                int i_gau_0 = std::max(n_wings - i_wl, 0);
                int i_gau_1 =
                  std::min(ckd.n_lbl + n_wings - i_wl - 1, n_gauss - 1);
                std::vector<double> this_gauss(gauss.begin() + i_gau_0,
                                               gauss.begin() + i_gau_1);
                double norm_gauss = // to make the sum unity for left and right
                                    // bounds
                  1
                  / (std::accumulate(
                    this_gauss.begin(), this_gauss.end(), 0.0));
                std::vector<double> this_signal(
                  l1.spectra[i_act].signal.begin() + i_lbl_0,
                  l1.spectra[i_act].signal.begin() + i_lbl_1);
                for (int k {}; k < this_gauss.size(); ++k) {
                    signal_conv[i_wl] +=
                      norm_gauss * this_gauss[k] * this_signal[k];
                }
            }
            l1.spectra[i_act].signal = std::move(signal_conv);
            
        }
    }
    l1.level = ProcLevel::swath;
}


auto drawOnDetector(const CKD& ckd, L1& l1) -> void
{
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
        std::vector<double> wl = ckd.wave.wavelength[i_act];
        const CubicSpline spl_wl_to_col { wl, cols };

        std::vector<double> lbl = (l1.spectra[i_act].signal);
        std::vector<double> observation_wls =( *l1.wavelength)[i_act];
        std::vector<double> col_ix(observation_wls.size(), 0.0); 
        // Calculate decimal col indices of line-by-line spectrum
        for (int i_wave {}; i_wave < (observation_wls.size()); ++i_wave) {
            col_ix[i_wave] = spl_wl_to_col.eval(observation_wls[i_wave]);
        }

        // Interpolate line-by-line spectrum to integer (detector) column range
        const CubicSpline spl_lbl_vs_col {col_ix, lbl};
        for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col){
            lbl_in_cols[i_act][i_col] = spl_lbl_vs_col.eval(i_col);
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

    l1.level = ProcLevel::l1b;
}


auto radiometric(const CKD& ckd, const bool enabled, L1& l1) -> void
{   
    if (!enabled) {
        return;
    }
    const double exposure_time_inv { 1.0 / l1.exposure_time };

    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] /= ckd.rad.rad[i] * exposure_time_inv;
        }
    }
    
    l1.level = ProcLevel::stray;
}


auto strayLight(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (std::isnan(l1.image[i])) {
            l1.image[i] = l1.image[i - 1];
        }
    }
    std::vector<double> conv_result(ckd.npix, 0.0);
    std::vector<std::complex<double>> image_fft(
      *std::ranges::max_element(ckd.stray.kernel_fft_sizes));
    for (int i_kernel {}; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
        std::vector<double> image_weighted { l1.image };
        for (int i {}; i < ckd.npix; ++i) {
            image_weighted[i] *= ckd.stray.weights[i_kernel][i];
        }
        const int image_n_rows {
            ckd.stray.edges[i_kernel * box::n + box::t]
            - ckd.stray.edges[i_kernel * box::n + box::b]
        };
        const int image_n_cols {
            ckd.stray.edges[i_kernel * box::n + box::r]
            - ckd.stray.edges[i_kernel * box::n + box::l]
        };
        std::vector<double> sub_image(image_n_rows * image_n_cols);
        for (int i {}; i < image_n_rows; ++i) {
            for (int j {}; j < image_n_cols; ++j) {
                sub_image[i * image_n_cols + j] = image_weighted
                  [(i + ckd.stray.edges[i_kernel * box::n + box::b])
                     * ckd.n_detector_cols
                   + j + ckd.stray.edges[i_kernel * box::n + box::l]];
            }
        }
        std::vector<double> conv_result_sub {};
        convolve(image_n_rows,
                 image_n_cols,
                 sub_image,
                 ckd.stray.kernel_rows[i_kernel],
                 ckd.stray.kernel_cols[i_kernel],
                 ckd.stray.kernels_fft[i_kernel],
                 image_fft,
                 conv_result_sub);
        for (int i {}; i < image_n_rows; ++i) {
            for (int j {}; j < image_n_cols; ++j) {
                conv_result[(i + ckd.stray.edges[i_kernel * box::n + box::b])
                              * ckd.n_detector_cols
                            + j
                            + ckd.stray.edges[i_kernel * box::n + box::l]] +=
                  conv_result_sub[i * image_n_cols + j];
            }
        }
    }
    std::vector<double> image_conv(ckd.npix);
    for (int i {}; i < ckd.npix; ++i) {
        image_conv[i] = (1.0 - ckd.stray.eta[i]) * l1.image[i] + conv_result[i];
    }
    l1.image = std::move(image_conv);
    l1.level = ProcLevel::prnu;
}

auto prnu(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] *= ckd.prnu.prnu[i];
        }
    }
    l1.level = ProcLevel::nonlin;
}

auto nonlinearity(const CKD& ckd,
                  const bool enabled,
                  const LinearSpline& nonlin_spline,
                  L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] = nonlin_spline.eval(l1.image[i]);
        }
    }
    l1.level = ProcLevel::noise;
}

auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] += ckd.dark.current[i] * l1.exposure_time;
        }
    }
}

auto noise(const CKD& ckd, const bool enabled, const int seed, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    static std::mt19937 gen { static_cast<unsigned long>(seed) };
    for (int i {}; i < ckd.npix; ++i) {
        const double noise_value { std::sqrt(ckd.noise.n2[i]
                                             + l1.image[i] * ckd.noise.g[i]) };
        std::normal_distribution<> d { 0.0, noise_value };
        l1.image[i] += d(gen);
    }
    l1.level = ProcLevel::dark;
}

auto darkOffset(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] += ckd.dark.offset[i];
        }
    }
    l1.level = ProcLevel::raw;
}

auto digitalToAnalog(const BinningTable& binning_table, L1& l1) -> void
{
    for (double& val : l1.image) {
        val *= l1.nr_coadditions;
    }
    // Note: The next line of code seems to be taking the information of binning
    // table into account, but at the end it divides by count_table (BF) That is
    // why in netxt steps there is a loop where multiplication by binSize (BF)
    // takes place
    binning_table.bin(l1.image);
    l1.image_i32.resize(binning_table.nBins());
    for (int i {}; i < binning_table.nBins(); ++i) {
        l1.image_i32[i] =
          static_cast<int>(std::round(l1.image[i] * binning_table.binSize(i)));
    }
    l1.image = std::vector<double> {};
    l1.level = ProcLevel::l1a;
}

} // namespace tango
