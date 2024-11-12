// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "forward_models.h"

#include <algorithm>
#include <random>
#include <tango_l1b/b_spline_2d.h>
#include <tango_l1b/binning_table.h>
#include <tango_l1b/ckd.h>
#include <tango_l1b/fourier.h>
#include <tango_l1b/l1.h>

namespace tango {

auto applyISRF(const CKD& ckd,
               const bool enabled,
               const Eigen::SparseMatrix<double>& isrf,
               L1& l1) -> void
{
    l1.level = ProcLevel::l1b;
    // If this process is disabled then linearly interpolate the
    // line-by-line spectra onto the CKD wavelength grids. We cannot
    // simply return like the other processes.
    const int n_waves_in { static_cast<int>((*l1.wavelengths).front().size()) };
    const int n_waves_out { static_cast<int>(ckd.swath.wavelengths.size()) };
    std::vector<double> spectra_out(ckd.n_act * n_waves_out);
    if (!enabled) {
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            const std::vector<double> spectrum {
                l1.spectra.begin() + i_act * n_waves_in,
                l1.spectra.begin() + (i_act + 1) * n_waves_in
            };
            const LinearSpline spline { (*l1.wavelengths)[i_act], spectrum };
            for (int i {}; i < n_waves_out; ++i) {
                spectra_out[i_act * n_waves_out + i] =
                  spline.eval(ckd.swath.wavelengths[i]);
            }
        }
        l1.spectra = std::move(spectra_out);
        return;
    }
    Eigen::VectorXd result(ckd.n_detector_cols);
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        result = isrf
                 * Eigen::Map<Eigen::VectorXd>(&l1.spectra[i_act * n_waves_in],
                                               n_waves_in);
        for (int i {}; i < n_waves_out; ++i) {
            spectra_out[i_act * n_waves_out + i] = result[i];
        }
    }
    l1.spectra = std::move(spectra_out);
}

auto radiometric(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::swath;
    if (!enabled) {
        return;
    }
    const double exposure_time_inv { 1.0 / l1.exposure_time };
    for (int i_act {}; i_act < ckd.n_act; i_act++) {
        for (int i {}; i < static_cast<int>(ckd.swath.wavelengths.size());
             ++i) {
            l1.spectra[i_act * ckd.swath.wavelengths.size() + i] /=
              ckd.rad.rad[i_act][0] * exposure_time_inv;
        }
    }
}

auto mapToDetector(const CKD& ckd, const int b_spline_order, L1& l1) -> void
{
    l1.level = ProcLevel::stray;
    const BSpline2D bspline_2d {
        b_spline_order, ckd.swath.act_angles, ckd.swath.wavelengths, l1.spectra
    };
    bspline_2d.eval(ckd.swath.act_map, ckd.swath.wavelength_map, l1.image);
    l1.stdev.resize(ckd.npix, 1.0);
}

auto strayLight(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::prnu;
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
}

auto prnu(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::nonlin;
    if (!enabled) {
        return;
    }
    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] *= ckd.prnu.prnu[i];
        }
    }
}

auto nonlinearity(const CKD& ckd,
                  const bool enabled,
                  const LinearSpline& nonlin_spline,
                  L1& l1) -> void
{
    l1.level = ProcLevel::dark_current;
    if (!enabled) {
        return;
    }
    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] = nonlin_spline.eval(l1.image[i]);
        }
    }
}

auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::noise;
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
    l1.level = ProcLevel::dark_offset;
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
}

auto darkOffset(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::raw;
    if (!enabled) {
        return;
    }
    for (int i {}; i < ckd.npix; ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] += ckd.dark.offset[i];
        }
    }
}

auto digitalToAnalog(const BinningTable& binning_table, L1& l1) -> void
{
    l1.level = ProcLevel::l1a;
    for (double& val : l1.image) {
        val *= l1.nr_coadditions;
    }
    binning_table.bin(l1.image);
    l1.image_i32.resize(binning_table.nBins());
    for (int i {}; i < binning_table.nBins(); ++i) {
        l1.image_i32[i] =
          static_cast<int>(std::round(l1.image[i] * binning_table.binSize(i)));
    }
    l1.image = std::vector<double> {};
}

} // namespace tango
