// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "uncalibration.h"

#include <algorithm>
#include <random>
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
    // Convolve the radiances with an ISRF from the line-by-line grid
    // onto the target grids taken from the CKD.
    const double sigma { fwhm_gauss / (2.0 * sqrt(2.0 * log(2.0))) };
    const double sigma_inv { 1.0 / (2.0 * sigma * sigma) };
    const double norm_inv { 1.0 / (sigma * std::sqrt(2 * std::numbers::pi)) };
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        const auto& lbl_wave { (*l1.wavelength)[i_act] };
        std::vector<double> signal_conv(ckd.n_detector_cols, 0.0);
        for (int i_wave {}; i_wave < ckd.n_detector_cols; ++i_wave) {
            for (int i_lbl {}; i_lbl < static_cast<int>(lbl_wave.size());
                 ++i_lbl) {
                // Convolution result for this CKD wavelength
                double conv { ckd.wave.wavelength[i_act][i_wave]
                              - lbl_wave[i_lbl] };
                // Reduce the computational cost by considering the
                // limited extent of the Gaussian.
                constexpr double wave_rel_threshold { 3.0 };
                if (std::abs(conv) > wave_rel_threshold * fwhm_gauss) {
                    continue;
                }
                conv *= conv * sigma_inv;
                conv = norm_inv * std::exp(-conv);
                double delta_lambda;
                if (i_lbl == 0) {
                    delta_lambda = lbl_wave[i_lbl + 1] - lbl_wave[i_lbl];
                } else if (i_lbl == static_cast<int>(lbl_wave.size()) - 1) {
                    delta_lambda =
                      lbl_wave.back() - lbl_wave[lbl_wave.size() - 2];
                } else {
                    delta_lambda =
                      0.5 * (lbl_wave[i_lbl + 1] - lbl_wave[i_lbl - 1]);
                }
                conv *= delta_lambda;
                signal_conv[i_wave] += conv * l1.spectra[i_act].signal[i_lbl];
            }
        }
        l1.spectra[i_act].signal = std::move(signal_conv);
    }
    l1.level = ProcLevel::rad;
}

auto radiometric(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    const double exposure_time_inv { 1.0 / l1.exposure_time };
    for (int i_act {}; i_act < ckd.n_act; i_act++) {
        for (int i {}; i < ckd.n_detector_cols; ++i) {
            l1.spectra[i_act].signal[i] /=
              ckd.rad.rad[i_act][i] * exposure_time_inv;
        }
    }
    l1.level = ProcLevel::swath;
}

auto drawOnDetector(const CKD& ckd, L1& l1) -> void
{
    l1.image.assign(ckd.npix, 0.0);
    std::vector<double> x_values(ckd.n_act);
    std::vector<double> y_values(ckd.n_act);
    for (int i_wave {}; i_wave < ckd.n_detector_cols; ++i_wave) {
        bool needReverse = false;
        if (ckd.swath.row_indices[ckd.n_act - 1 ][i_wave] < ckd.swath.row_indices[0][i_wave]){
            needReverse = true;
        }
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            const int act_idx { static_cast<int>(i_act) };
            if (needReverse){
                const int act_idx { static_cast<int>(ckd.n_act - 1 - i_act) };
            }

            x_values[i_act] = ckd.swath.row_indices[act_idx][i_wave];
            y_values[i_act] = l1.spectra[act_idx].signal[i_wave];
        }
        CubicSpline spline { x_values, y_values };
        for (int i_spat {}; i_spat < ckd.n_detector_rows; ++i_spat) {
            l1.image[i_spat * ckd.n_detector_cols + i_wave] =
              spline.eval(i_spat);
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
    // Note: The next line of code seems to be taking the information of binning table into account,
    // but at the end it divides by count_table (BF)
    // That is why in netxt steps there is a loop where multiplication by binSize (BF) takes place
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
