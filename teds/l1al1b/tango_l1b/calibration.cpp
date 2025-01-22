// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "calibration.h"

#include "b_spline_2d.h"
#include "binning_table.h"
#include "ckd.h"
#include "cubic_spline.h"
#include "fourier.h"
#include "geometry.h"
#include "l1.h"
#include "solar_model.h"

#include <algorithm>
#include <numeric>
#include <span>

namespace tango {

auto binScaling(const CKD& ckd,
                const BinningTable& binning_table,
                L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::raw;
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i {}; i < ckd.npix_binned; ++i) {
            if (!ckd.pixel_mask[i]) {
                l1_prod.signal[i_alt * ckd.npix_binned + i] /=
                  l1_prod.nr_coadditions * binning_table.binSize(i);
            }
        }
    }
    l1_prod.noise.assign(l1_prod.signal.size(),
                         static_cast<double>(l1_prod.nr_coadditions));
    l1_prod.nr_coadditions = 1;
}

auto darkOffset(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::dark_offset;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i {}; i < ckd.npix_binned; ++i) {
            if (!ckd.pixel_mask[i]) {
                l1_prod.signal[i_alt * ckd.npix_binned + i] -=
                  ckd.dark.offset[i];
            }
        }
    }
}

auto noise(const CKD& ckd,
           const bool enabled,
           const BinningTable& binning_table,
           L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::noise;
    if (!enabled) {
        std::ranges::fill(l1_prod.noise, 1.0);
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i {}; i < ckd.npix_binned; ++i) {
            if (!ckd.pixel_mask[i]) {
                const int idx { i_alt * ckd.npix_binned + i };
                const double var { ckd.noise.g[i]
                                     * std::abs(l1_prod.signal[idx])
                                   + ckd.noise.n2[i] };
                const double coad { l1_prod.noise[idx] };
                l1_prod.noise[idx] =
                  std::sqrt(var / (coad * binning_table.binSize(i)));
            }
        }
    }
}

auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::dark_current;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i {}; i < ckd.npix_binned; ++i) {
            if (!ckd.pixel_mask[i]) {
                l1_prod.signal[i_alt * ckd.npix_binned + i] -=
                  ckd.dark.current[i] * l1_prod.exposure_time;
            }
        }
    }
}

auto nonlinearity(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::nonlin;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i {}; i < ckd.npix_binned; ++i) {
            if (!ckd.pixel_mask[i]) {
                const int idx { i_alt * ckd.npix_binned + i };
                l1_prod.noise[idx] *=
                  ckd.nonlin.spline.deriv(l1_prod.signal[idx]);
                l1_prod.signal[idx] =
                  ckd.nonlin.spline.eval(l1_prod.signal[idx]);
            }
        }
    }
}

auto prnu(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::prnu;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i {}; i < ckd.npix_binned; ++i) {
            if (!ckd.pixel_mask[i]) {
                const int idx { i_alt * ckd.npix_binned + i };
                l1_prod.signal[idx] /= ckd.prnu.prnu[i];
                l1_prod.noise[idx] /= ckd.prnu.prnu[i];
            }
        }
    }
}

auto removeBadValues(const CKD& ckd, L1& l1_prod) -> void
{
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        std::vector<double> x_values {};
        std::vector<double> y_values_signal {};
        std::vector<double> y_values_noise {};
        for (int i_row {}; i_row < ckd.n_detector_rows_binned; ++i_row) {
            x_values.clear();
            y_values_signal.clear();
            y_values_noise.clear();
            for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col) {
                const int idx { i_row * ckd.n_detector_cols + i_col };
                if (!ckd.pixel_mask[idx]) {
                    const int idx_full { i_alt * ckd.npix_binned + idx };
                    x_values.push_back(static_cast<double>(i_col));
                    y_values_signal.push_back(l1_prod.signal[idx_full]);
                    y_values_noise.push_back(l1_prod.noise[idx_full]);
                }
            }
            const CubicSpline spline_signal { x_values, y_values_signal };
            for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col) {
                l1_prod.signal[i_alt * ckd.npix_binned
                               + i_row * ckd.n_detector_cols + i_col] =
                  spline_signal.eval(static_cast<double>(i_col));
            }
            const CubicSpline spline_noise { x_values, y_values_noise };
            for (int i_col {}; i_col < ckd.n_detector_cols; ++i_col) {
                l1_prod.noise[i_alt * ckd.npix_binned
                              + i_row * ckd.n_detector_cols + i_col] =
                  spline_noise.eval(static_cast<double>(i_col));
            }
        }
    }
}

static auto vanCittertMulti(const CKD& ckd,
                            std::vector<std::complex<double>>& signal_fft,
                            std::vector<double>& signal_ideal,
                            std::vector<double>& sub_signal,
                            std::vector<double>& conv_result_sub,
                            std::vector<double>& conv_result) -> void
{
    conv_result.assign(ckd.npix, 0.0);
    for (int i_kernel {}; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
        // Number of rows and column in this subsignal
        const int signal_n_rows {
            ckd.stray.edges[i_kernel * box::n + box::t]
            - ckd.stray.edges[i_kernel * box::n + box::b]
        };
        const int signal_n_cols {
            ckd.stray.edges[i_kernel * box::n + box::r]
            - ckd.stray.edges[i_kernel * box::n + box::l]
        };
        // Each iteration starts with signal_ideal as the best
        // current estimate.
        sub_signal.resize(signal_n_rows * signal_n_cols);
        for (int i {}; i < signal_n_rows; ++i) {
            for (int j {}; j < signal_n_cols; ++j) {
                const int idx {
                    (i + ckd.stray.edges[i_kernel * box::n + box::b])
                      * ckd.n_detector_cols
                    + j + ckd.stray.edges[i_kernel * box::n + box::l]
                };
                sub_signal[i * signal_n_cols + j] =
                  signal_ideal[idx] * ckd.stray.weights[i_kernel][idx];
            }
        }
        // Result of taking a convolution using one of the
        // subsignals and kernels
        convolve(signal_n_rows,
                 signal_n_cols,
                 sub_signal,
                 ckd.stray.kernel_rows[i_kernel],
                 ckd.stray.kernel_cols[i_kernel],
                 ckd.stray.kernels_fft[i_kernel],
                 signal_fft,
                 conv_result_sub);
        // The full convolution is a sum over all convolutions
        for (int i {}; i < signal_n_rows; ++i) {
            for (int j {}; j < signal_n_cols; ++j) {
                conv_result[(i + ckd.stray.edges[i_kernel * box::n + box::b])
                              * ckd.n_detector_cols
                            + j
                            + ckd.stray.edges[i_kernel * box::n + box::l]] +=
                  conv_result_sub[i * signal_n_cols + j];
            }
        }
    }
}

static auto vanCittertSingle(const CKD& ckd,
                             std::vector<std::complex<double>>& signal_fft,
                             std::vector<double>& signal_ideal,
                             std::vector<double>& conv_result) -> void
{
    convolve(ckd.n_detector_rows,
             ckd.n_detector_cols,
             signal_ideal,
             ckd.stray.kernel_rows.front(),
             ckd.stray.kernel_cols.front(),
             ckd.stray.kernels_fft.front(),
             signal_fft,
             conv_result);
}

auto strayLight(const CKD& ckd,
                const BinningTable& binning_table,
                const int n_van_cittert,
                L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::stray;
    if (n_van_cittert == 0) {
        return;
    }
    // Pre-allocate the signal FFT for efficiency
    std::vector<std::complex<double>> signal_fft(
      *std::ranges::max_element(ckd.stray.kernel_fft_sizes));
    // Result of convolution
    std::vector<double> conv_result {};
    // With multiple kernels, part of the unbinned detector signal
    std::vector<double> sub_signal {};
    // With multiple kernels, result of convolving one kernel with
    // a subsignal
    std::vector<double> conv_result_sub {};
#pragma omp parallel for firstprivate(signal_fft) private(                     \
    conv_result, sub_signal, conv_result_sub)
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        // Unbin the detector image
        const std::span signal(&l1_prod.signal[i_alt * ckd.npix_binned],
                               ckd.npix_binned);
        std::vector<double> signal_unbin(ckd.npix);
        binning_table.unbin(signal, signal_unbin);
        // "Ideal" signal, i.e. the one without stray light
        std::vector<double> signal_ideal { signal_unbin };
        // Van Cittert algorithm
        for (int i_vc {}; i_vc < n_van_cittert; ++i_vc) {
            if (ckd.stray.n_kernels > 1) {
                vanCittertMulti(ckd,
                                signal_fft,
                                signal_ideal,
                                sub_signal,
                                conv_result_sub,
                                conv_result);
            } else {
                vanCittertSingle(ckd, signal_fft, signal_ideal, conv_result);
            }
            for (int i {}; i < ckd.npix; ++i) {
                signal_ideal[i] =
                  (signal_unbin[i] - conv_result[i]) / (1 - ckd.stray.eta[i]);
            }
        }
        std::vector<double> signal_binned(ckd.npix_binned);
        binning_table.bin(signal_ideal, signal_binned);
        std::copy(signal_binned.cbegin(), signal_binned.cend(), signal.begin());
    }
}

auto mapFromDetector(const CKD& ckd,
                     const BinningTable& binning_table,
                     const int b_spline_order,
                     const bool exact_drawing,
                     L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::swath;
    // Size of spectra per ALT position
    const size_t spec_size { ckd.n_act * ckd.swath.wavelengths.size() };
    l1_prod.spectra.assign(l1_prod.n_alt * spec_size, 0.0);
    l1_prod.spectra_noise.assign(l1_prod.n_alt * spec_size, 0.0);
    // Assume there is binning only across rows. From that determine
    // the number of rows of the binned detector image.
    std::vector<double> rows(ckd.n_detector_rows_binned, 0.0);
    for (int i {}; i < ckd.n_detector_rows_binned; ++i) {
        const double cur_bin_size { static_cast<double>(
          binning_table.binSize(i * ckd.n_detector_cols)) };
        if (i == 0) {
            rows[i] = 0.5 * cur_bin_size - 0.5;
        } else {
            const double prev_bin_size { static_cast<double>(
              binning_table.binSize((i - 1) * ckd.n_detector_cols)) };
            const double bin_size_half { 0.5 * (prev_bin_size - cur_bin_size) };
            rows[i] = rows[i - 1] + cur_bin_size + bin_size_half;
        }
    }
    std::vector<double> cols(ckd.n_detector_cols);
    std::iota(cols.begin(), cols.end(), 0.0);
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        const BSpline2D bspline_2d_signal {
            b_spline_order, rows, cols, &l1_prod.signal[i_alt * ckd.npix_binned]
        };
        bspline_2d_signal.eval(ckd.swath.row_map,
                               ckd.swath.col_map,
                               &l1_prod.spectra[i_alt * spec_size]);
        const BSpline2D bspline_2d_noise {
            b_spline_order, rows, cols, &l1_prod.noise[i_alt * ckd.npix_binned]
        };
        bspline_2d_noise.eval(ckd.swath.row_map,
                              ckd.swath.col_map,
                              &l1_prod.spectra_noise[i_alt * spec_size]);
    }
    if (!exact_drawing) {
        l1_prod.signal = std::vector<double> {};
        l1_prod.noise = std::vector<double> {};
        l1_prod.wavelengths =
          std::vector<std::vector<double>>(ckd.n_act, ckd.swath.wavelengths);
        return;
    }
    // Regrid row_map from intermediate wavelengths to wavelengths
    // corresponding to detector columns.
    std::vector<double> row_map(ckd.n_act * ckd.n_detector_cols);
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        const CubicSpline spline {
            ckd.swath.wavelengths,
            { ckd.swath.row_map.cbegin() + i_act * ckd.swath.wavelengths.size(),
              ckd.swath.row_map.cbegin()
                + (i_act + 1) * ckd.swath.wavelengths.size() }
        };
        for (int i {}; i < ckd.n_detector_cols; ++i) {
            row_map[i_act * ckd.n_detector_cols + i] =
              spline.eval(ckd.wave.wavelengths[i_act][i]);
        }
    }
    l1_prod.spectra.resize(l1_prod.n_alt * ckd.n_act * ckd.n_detector_cols);
    l1_prod.spectra_noise.resize(l1_prod.n_alt * ckd.n_act
                                 * ckd.n_detector_cols);
#pragma omp parallel for
    for (size_t i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            for (int i_wave {}; i_wave < ckd.n_detector_cols; ++i_wave) {
                const int row_dn { static_cast<int>(
                  row_map[i_act * ckd.n_detector_cols + i_wave]) };
                const int row_up { row_dn + 1 };
                const double weight {
                    row_map[i_act * ckd.n_detector_cols + i_wave] - row_dn
                };
                const int pix_dn { row_dn * ckd.n_detector_cols + i_wave };
                const int pix_up { row_up * ckd.n_detector_cols + i_wave };
                const size_t idx {
                    (i_alt * ckd.n_act + i_act) * ckd.n_detector_cols + i_wave
                };
                l1_prod.spectra[idx] =
                  weight * l1_prod.signal[i_alt * ckd.npix + pix_dn]
                  + (1 - weight) * l1_prod.signal[i_alt * ckd.npix + pix_up];

                const double a { weight
                                 * l1_prod.noise[i_alt * ckd.npix + pix_dn] };
                const double b { (1 - weight)
                                 * l1_prod.noise[i_alt * ckd.npix + pix_up] };
                l1_prod.spectra_noise[idx] = std::sqrt(a * a + b * b);
            }
        }
    }
    l1_prod.signal = std::vector<double> {};
    l1_prod.noise = std::vector<double> {};
    l1_prod.wavelengths = ckd.wave.wavelengths;
}

auto changeWavelengthGrid(const CKD& ckd, L1& l1_prod) -> void
{
    const int n_waves_in { static_cast<int>(ckd.swath.wavelengths.size()) };
    std::vector<double> spectra_out(l1_prod.n_alt * ckd.n_act
                                    * ckd.n_detector_cols);
    std::vector<double> spectra_noise_out(spectra_out.size());
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            const int idx_beg_in { (i_alt * ckd.n_act + i_act) * n_waves_in };
            const int idx_beg_out { (i_alt * ckd.n_act + i_act)
                                    * ckd.n_detector_cols };
            const CubicSpline spline { ckd.swath.wavelengths,
                                       { l1_prod.spectra.begin() + idx_beg_in,
                                         l1_prod.spectra.begin() + idx_beg_in
                                           + n_waves_in } };
            for (int i {}; i < ckd.n_detector_cols; ++i) {
                spectra_out[idx_beg_out + i] =
                  spline.eval(ckd.wave.wavelengths[i_act][i]);
            }
            const CubicSpline spline_noise {
                ckd.swath.wavelengths,
                { l1_prod.spectra_noise.begin() + idx_beg_in,
                  l1_prod.spectra_noise.begin() + idx_beg_in + n_waves_in }
            };
            for (int i {}; i < ckd.n_detector_cols; ++i) {
                spectra_noise_out[idx_beg_out + i] =
                  spline_noise.eval(ckd.wave.wavelengths[i_act][i]);
            }
        }
    }
    l1_prod.spectra = std::move(spectra_out);
    l1_prod.spectra_noise = std::move(spectra_noise_out);
    l1_prod.wavelengths = ckd.wave.wavelengths;
}

auto radiometric(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    // Radiometric calibration is the last step to get to L1B
    l1_prod.level = ProcLevel::l1b;
    if (!enabled) {
        return;
    }
    const double exposure_time_inv { 1.0 / l1_prod.exposure_time };
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            for (int i {}; i < ckd.n_detector_cols; ++i) {
                const int idx {
                    (i_alt * ckd.n_act + i_act) * ckd.n_detector_cols + i
                };
                l1_prod.spectra[idx] *=
                  ckd.rad.rad[i_act][i] * exposure_time_inv;
                l1_prod.spectra_noise[idx] *=
                  ckd.rad.rad[i_act][i] * exposure_time_inv;
            }
        }
    }
}

// Bin a multi-dimensional array. Stride is the distance between
// consecutive elements that should be averaged into one bin.
static auto binNDArray(const int bin,
                       const int stride,
                       std::vector<double>& array) -> void
{
    const int n_rows { static_cast<int>(array.size() / stride) };
    const int n_cols { stride };
    const int n_rows_binned { n_rows / bin };
    std::vector<double> array_binned(n_rows_binned * n_cols, 0.0);
    for (int i {}; i < n_rows_binned; ++i) {
        for (int j {}; j < n_cols; ++j) {
            for (int k {}; k < bin; ++k) {
                array_binned[i * n_cols + j] +=
                  array[(i * bin + k) * n_cols + j];
            }
        }
    }
    for (double& x : array_binned) {
        x /= bin;
    }
    array = std::move(array_binned);
}

auto binL1B(const int bin, L1& l1_prod) -> void
{
    const int n_waves { static_cast<int>(l1_prod.wavelengths.front().size()) };
    binNDArray(bin, n_waves, l1_prod.spectra);
    binNDArray(bin, n_waves, l1_prod.spectra_noise);
    binNDArray(bin, 1, l1_prod.geo.lat);
    binNDArray(bin, 1, l1_prod.geo.lon);
    binNDArray(bin, 1, l1_prod.geo.height);
    binNDArray(bin, 1, l1_prod.geo.sza);
    binNDArray(bin, 1, l1_prod.geo.saa);
    binNDArray(bin, 1, l1_prod.geo.vza);
    binNDArray(bin, 1, l1_prod.geo.vaa);
    // Wavelength array needs to be flattened first
    std::vector<double> wavelengths_flat(l1_prod.wavelengths.size() * n_waves);
    for (int i {}; i < static_cast<int>(l1_prod.wavelengths.size()); ++i) {
        for (int j {}; j < n_waves; ++j) {
            wavelengths_flat[i * n_waves + j] = l1_prod.wavelengths[i][j];
        }
    }
    binNDArray(bin, n_waves, wavelengths_flat);
    l1_prod.wavelengths.resize(l1_prod.wavelengths.size() / bin);
    for (int i {}; i < static_cast<int>(l1_prod.wavelengths.size()); ++i) {
        for (int j {}; j < n_waves; ++j) {
            l1_prod.wavelengths[i][j] = wavelengths_flat[i * n_waves + j];
        }
    }
}

} // namespace tango
