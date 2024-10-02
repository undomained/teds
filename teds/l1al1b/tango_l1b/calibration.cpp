// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "calibration.h"

#include "binning_table.h"
#include "ckd.h"
#include "cubic_spline.h"
#include "fourier.h"
#include "l1.h"

#include <algorithm>

namespace tango {

// Smooth over bad values. This is necessary for algorithms such as
// stray light correction which use all pixels.
static auto fillHoles(const std::vector<bool>& pixel_mask,
                      std::vector<double>& image) -> void
{
    const int i_last { static_cast<int>(image.size() - 1) };
    for (int i {}; i < static_cast<int>(image.size()); ++i) {
        if (pixel_mask[i]) {
            // Unless we are at either end of the image array or one
            // or both of the neighboring pixels are bad, take the
            // average of the neighboring values.
            int n_good {};
            image[i] = 0.0;
            if (i > 0 && !pixel_mask[i - 1]) {
                image[i] += image[i - 1];
                ++n_good;
            }
            if (i < i_last && !pixel_mask[i + 1]) {
                image[i] += image[i + 1];
                ++n_good;
            }
            image[i] /= std::max(1, n_good);
        }
    }
}

auto binningTable(const BinningTable& binning_table,
                  const Unbin unbin,
                  const int n_rows,
                  const int n_cols,
                  L1& l1) -> void
{
    auto pixel_mask { l1.pixel_mask };
    if (unbin != Unbin::none) {
        binning_table.bin(pixel_mask);
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!pixel_mask[i]) {
            l1.image[i] /= binning_table.binSize(i);
        }
    }
    l1.level = ProcLevel::raw;
    if (unbin == Unbin::none) {
        return;
    }
    // This algorithm assumes no binning across columns
    std::vector<double> row_indices {};
    for (int i_row {}; i_row < n_rows;) {
        const int bin_size { binning_table.binSize(
          binning_table.binIndex(i_row * n_cols)) };
        row_indices.push_back(i_row + 0.5 * bin_size - 0.5);
        i_row += bin_size;
    }
    const int n_rows_binned { static_cast<int>(row_indices.size()) };
    std::vector<double> row_values(n_rows_binned);
    std::vector<double> image_unbinned(n_rows * n_cols);
    if (unbin == Unbin::nearest) {
        binning_table.unbin(l1.image, image_unbinned);
    } else {
        fillHoles(pixel_mask, l1.image);
        for (int i_col {}; i_col < n_cols; ++i_col) {
            for (int i_row {}; i_row < n_rows_binned; ++i_row) {
                row_values[i_row] = l1.image[i_row * n_cols + i_col];
            }
            if (unbin == Unbin::linear) {
                const LinearSpline spline { row_indices, row_values };
                for (int i_row {}; i_row < n_rows; ++i_row) {
                    image_unbinned[i_row * n_cols + i_col] =
                      spline.eval(static_cast<double>(i_row));
                }
            } else if (unbin == Unbin::cubic) {
                const CubicSpline spline { row_indices, row_values };
                for (int i_row {}; i_row < n_rows; ++i_row) {
                    image_unbinned[i_row * n_cols + i_col] =
                      spline.eval(static_cast<double>(i_row));
                }
            }
        }
    }
    l1.image = std::move(image_unbinned);
    l1.stdev.resize(l1.image.size());
    // Image dimensions have changed and this should be reflected in
    // the binning_table NetCDF attribute in case such a product is
    // produced.
    l1.binning_table_id = 0;
}

auto darkOffset(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            l1.image[i] -= ckd.dark.offset[i];
        }
    }
}

auto noise(const CKD& ckd,
           const bool enabled,
           const BinningTable& binning_table,
           L1& l1) -> void
{
    if (!enabled) {
        l1.stdev.assign(l1.image.size(), 1.0);
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (l1.pixel_mask[i]) {
            continue;
        }
        const double var { ckd.noise.g[i] * l1.image[i] + ckd.noise.n2[i] };
        if (var <= 0.0) {
            l1.pixel_mask[i] = true;
        } else if (l1.binning_table_id != 0) {
            l1.stdev[i] =
              std::sqrt(var / (l1.nr_coadditions * binning_table.binSize(i)));
        } else {
            l1.stdev[i] = std::sqrt(var / l1.nr_coadditions);
        }
    }
}

auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            l1.image[i] -= ckd.dark.current[i] * l1.exposure_time;
        }
    }
    l1.level = ProcLevel::dark;
}

auto nonlinearity(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            l1.stdev[i] *= ckd.nonlin.spline.deriv(l1.image[i]);
            l1.image[i] = ckd.nonlin.spline.eval(l1.image[i]);
        }
    }
    l1.level = ProcLevel::nonlin;
}

auto prnu(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
        if (!l1.pixel_mask[i]) {
            l1.image[i] /= ckd.prnu.prnu[i];
            l1.stdev[i] /= ckd.prnu.prnu[i];
        }
    }
    l1.level = ProcLevel::prnu;
}

auto strayLight(const CKD& ckd,
                const BinningTable& binning_table,
                const int n_van_cittert,
                L1& l1) -> void
{
    if (n_van_cittert == 0) {
        return;
    }
    fillHoles(l1.pixel_mask, l1.image);
    // Unbin the image using the stray light binning table
    std::vector<double> image_unbin(ckd.npix);
    if (l1.binning_table_id == 0) {
        image_unbin = std::move(l1.image);
    } else {
        binning_table.unbin(l1.image, image_unbin);
    }
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
    if (l1.binning_table_id == 0) {
        l1.image = std::move(image_ideal);
    } else {
        binning_table.bin(image_ideal, l1.image);
    }
    l1.level = ProcLevel::stray;
}

auto extractSpectra(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    l1.spectra.resize(ckd.n_act);
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        l1.spectra[i_act].extract(
          ckd, l1.image, l1.stdev, l1.pixel_mask, i_act);
    }
    l1.image = std::vector<double> {};
    l1.stdev = std::vector<double> {};
    l1.level = ProcLevel::swath;
}

auto radiometric(const CKD& ckd, const bool enabled, L1& l1) -> void
{
    if (!enabled) {
        return;
    }
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        l1.spectra[i_act].calibrate(ckd, l1.exposure_time, i_act);
    }
    // Radiometric calibration is the last step to get to L1B
    l1.level = ProcLevel::l1b;
}

} // namespace tango
