// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "calibration.h"

#include "b_spline_2d.h"
#include "binning_table.h"
#include "ckd.h"
#include "cubic_spline.h"
#include "fourier.h"
#include "l1.h"

#include <algorithm>
#include <numeric>

namespace tango {

auto binningTable(const CKD& ckd,
                  const BinningTable& binning_table,
                  L1& l1) -> void
{
    l1.level = ProcLevel::raw;
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] /= binning_table.binSize(i);
        }
    }
}

auto darkOffset(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::dark_offset;
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] -= ckd.dark.offset[i];
        }
    }
}

auto noise(const CKD& ckd,
           const bool enabled,
           const BinningTable& binning_table,
           L1& l1) -> void
{
    l1.level = ProcLevel::noise;
    if (!enabled) {
        l1.stdev.assign(l1.image.size(), 1.0);
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!ckd.pixel_mask[i]) {
            const double var { ckd.noise.g[i] * std::abs(l1.image[i])
                               + ckd.noise.n2[i] };
            l1.stdev[i] =
              std::sqrt(var / (l1.nr_coadditions * binning_table.binSize(i)));
        }
    }
}

auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::dark_current;
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] -= ckd.dark.current[i] * l1.exposure_time;
        }
    }
}

auto nonlinearity(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::nonlin;
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.stdev[i] *= ckd.nonlin.spline.deriv(l1.image[i]);
            l1.image[i] = ckd.nonlin.spline.eval(l1.image[i]);
        }
    }
}

auto prnu(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.level = ProcLevel::prnu;
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!ckd.pixel_mask[i]) {
            l1.image[i] /= ckd.prnu.prnu[i];
            l1.stdev[i] /= ckd.prnu.prnu[i];
        }
    }
}

auto removeBadValues(const CKD& ckd, L1& l1) -> void
{
    std::vector<double> x_values {};
    std::vector<double> y_values_image {};
    std::vector<double> y_values_stdev {};
    for (int i_row {}; i_row < ckd.n_detector_rows; ++i_row) {
        x_values.clear();
        y_values_image.clear();
        y_values_stdev.clear();
        for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col) {
            const int idx { i_row * ckd.n_detector_cols + i_col };
            if (!ckd.pixel_mask[idx]) {
                x_values.push_back(static_cast<double>(i_col));
                y_values_image.push_back(l1.image[idx]);
                y_values_stdev.push_back(l1.stdev[idx]);
            }
        }
        const CubicSpline spline_image { x_values, y_values_image };
        for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col) {
            l1.image[i_row * ckd.n_detector_cols + i_col] =
              spline_image.eval(static_cast<double>(i_col));
        }
        const CubicSpline spline_stdev { x_values, y_values_stdev };
        for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col) {
            l1.stdev[i_row * ckd.n_detector_cols + i_col] =
              spline_stdev.eval(static_cast<double>(i_col));
        }
    }
}

auto strayLight(const CKD& ckd,
                const BinningTable& binning_table,
                const int n_van_cittert,
                L1& l1) -> void
{
    l1.level = ProcLevel::stray;
    if (n_van_cittert == 0) {
        return;
    }
    // Unbin the image using the stray light binning table
    std::vector<double> image_unbin(ckd.npix);
    binning_table.unbin(l1.image, image_unbin);
    // "Ideal" image, i.e. the one without stray light
    std::vector<double> image_ideal { image_unbin };
    // Part of the unbinned detector image
    std::vector<double> sub_image {};
    // Result of convoling one kernel with a subimage
    std::vector<double> conv_result_sub {};
    // Allocate the image FFT once to save time
    std::vector<std::complex<double>> image_fft(
      *std::ranges::max_element(ckd.stray.kernel_fft_sizes));
    // Result of taking the whole convolution
    std::vector<double> conv_result(ckd.npix);
    // Van Cittert algorithm
    for (int i_vc {}; i_vc < n_van_cittert; ++i_vc) {
        std::ranges::fill(conv_result, 0.0);
        for (int i_kernel {}; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
            // Number of rows and column in this subimage
            const int image_n_rows {
                ckd.stray.edges[i_kernel * box::n + box::t]
                - ckd.stray.edges[i_kernel * box::n + box::b]
            };
            const int image_n_cols {
                ckd.stray.edges[i_kernel * box::n + box::r]
                - ckd.stray.edges[i_kernel * box::n + box::l]
            };
            // Each iteration starts with image_ideal as the best
            // current estimate.
            sub_image.resize(image_n_rows * image_n_cols);
            for (int i {}; i < image_n_rows; ++i) {
                for (int j {}; j < image_n_cols; ++j) {
                    const int idx {
                        (i + ckd.stray.edges[i_kernel * box::n + box::b])
                          * ckd.n_detector_cols
                        + j + ckd.stray.edges[i_kernel * box::n + box::l]
                    };
                    sub_image[i * image_n_cols + j] =
                      image_ideal[idx] * ckd.stray.weights[i_kernel][idx];
                }
            }
            // Result of taking a convolution using one of the
            // subimages and kernels
            convolve(image_n_rows,
                     image_n_cols,
                     sub_image,
                     ckd.stray.kernel_rows[i_kernel],
                     ckd.stray.kernel_cols[i_kernel],
                     ckd.stray.kernels_fft[i_kernel],
                     image_fft,
                     conv_result_sub);
            // The full convolution is a sum over all convolutions
            for (int i {}; i < image_n_rows; ++i) {
                for (int j {}; j < image_n_cols; ++j) {
                    conv_result
                      [(i + ckd.stray.edges[i_kernel * box::n + box::b])
                         * ckd.n_detector_cols
                       + j + ckd.stray.edges[i_kernel * box::n + box::l]] +=
                      conv_result_sub[i * image_n_cols + j];
                }
            }
        }
        for (int i {}; i < ckd.npix; ++i) {
            image_ideal[i] =
              (image_unbin[i] - conv_result[i]) / (1 - ckd.stray.eta[i]);
        }
    }
    binning_table.bin(image_ideal, l1.image);
}

auto mapFromDetector(const CKD& ckd, const int b_spline_order, L1& l1) -> void
{
    l1.level = ProcLevel::swath;
    std::vector<double> rows(ckd.n_detector_rows);
    std::vector<double> cols(ckd.n_detector_cols);
    std::iota(rows.begin(), rows.end(), 0.0);
    std::iota(cols.begin(), cols.end(), 0.0);
    const BSpline2D bspline_2d_image { b_spline_order, rows, cols, l1.image };
    bspline_2d_image.eval(ckd.swath.row_map, ckd.swath.col_map, l1.spectra);
    const BSpline2D bspline_2d_noise { b_spline_order, rows, cols, l1.stdev };
    bspline_2d_noise.eval(
      ckd.swath.row_map, ckd.swath.col_map, l1.spectra_stdev);
    l1.image = std::vector<double> {};
    l1.stdev = std::vector<double> {};
}

auto changeWavelengthGrid(const CKD& ckd, L1& l1) -> void
{
    const int n_waves_in { static_cast<int>(ckd.swath.wavelengths.size()) };
    std::vector<double> spectra_out(ckd.n_act * ckd.n_detector_cols);
    std::vector<double> spectra_stdev_out(ckd.n_act * ckd.n_detector_cols);
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        const std::vector<double> spectrum {
            l1.spectra.begin() + i_act * n_waves_in,
            l1.spectra.begin() + (i_act + 1) * n_waves_in
        };
        const CubicSpline spline { ckd.swath.wavelengths, spectrum };
        for (int i {}; i < ckd.n_detector_cols;++i) {
            spectra_out[i_act * ckd.n_detector_cols + i] =
              spline.eval(ckd.wave.wavelengths[i_act][i]);
        }
    }
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        const std::vector<double> spectrum {
            l1.spectra_stdev.begin() + i_act * n_waves_in,
            l1.spectra_stdev.begin() + (i_act + 1) * n_waves_in
        };
        const CubicSpline spline { ckd.swath.wavelengths, spectrum };
        for (int i {}; i < ckd.n_detector_cols;++i) {
            spectra_stdev_out[i_act * ckd.n_detector_cols + i] =
              spline.eval(ckd.wave.wavelengths[i_act][i]);
        }
    }
    l1.spectra = std::move(spectra_out);
    l1.spectra_stdev = std::move(spectra_stdev_out);
}

auto radiometric(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    // Radiometric calibration is the last step to get to L1B
    l1.level = ProcLevel::l1b;
    if (!enabled) {
        return;
    }
    const double exposure_time_inv { 1.0 / l1.exposure_time };
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        for (int i {}; i < ckd.n_detector_cols; ++i) {
            l1.spectra[i_act * ckd.n_detector_cols + i] *=
              ckd.rad.rad[i_act][i] * exposure_time_inv;
            l1.spectra_stdev[i_act * ckd.n_detector_cols + i] *=
              ckd.rad.rad[i_act][i] * exposure_time_inv;
        }
    }
}

} // namespace tango
