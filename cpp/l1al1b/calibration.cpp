// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "calibration.h"

#include <common/b_spline_2d.h>
#include <common/binning_table.h>
#include <common/ckd.h>
#include <common/cubic_spline.h>
#include <common/fourier.h>
#include <common/l1.h>

namespace tango {

auto binScaling(const BinningTable& binning_table, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::raw;
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        l1_prod.signal.row(i_alt) /=
          l1_prod.nr_coadditions * binning_table.count_table;
    }
    l1_prod.noise =
      ArrayXXd::Constant(l1_prod.signal.rows(),
                         l1_prod.signal.cols(),
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
        l1_prod.signal.row(i_alt) -= ckd.dark.offset;
    }
}

auto noise(const CKD& ckd,
           const bool enabled,
           const BinningTable& binning_table,
           L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::noise;
    if (!enabled) {
        l1_prod.noise = 1.0;
        return;
    }
    // Element-wise multiplication does not work with column-major
    // vectors ArrayXXd.
    typedef Eigen::Array<double, 1, Eigen::Dynamic, Eigen::RowMajor> ArrayXd;
    const Eigen::Ref<const ArrayXd> g(ckd.noise.g);
    const Eigen::Ref<const ArrayXd> n2(ckd.noise.n2);
    const Eigen::Ref<const ArrayXd> count_table(binning_table.count_table);
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        const auto& var(g * Eigen::abs(l1_prod.signal.row(i_alt)) + n2);
        const auto& coad(l1_prod.noise.row(i_alt));
        l1_prod.noise.row(i_alt) = Eigen::sqrt(var / (coad * count_table));
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
        l1_prod.signal.row(i_alt) -= ckd.dark.current * l1_prod.exposure_time;
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
        l1_prod.noise.row(i_alt) *=
          ckd.nonlin.spline.deriv(l1_prod.signal.row(i_alt));
        l1_prod.signal.row(i_alt) =
          ckd.nonlin.spline.eval(l1_prod.signal.row(i_alt));
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
        l1_prod.signal.row(i_alt) /= ckd.prnu.prnu;
        l1_prod.noise.row(i_alt) /= ckd.prnu.prnu;
    }
}

auto removeBadValues(const CKD& ckd, L1& l1_prod) -> void
{
    const int n_cols { ckd.n_detector_cols };
    const Eigen::ArrayXd col_indices =
      Eigen::ArrayXd::LinSpaced(n_cols, 0, n_cols - 1);
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        Eigen::ArrayXd x_values(n_cols);
        Eigen::ArrayXd y_values_signal(n_cols);
        Eigen::ArrayXd y_values_noise(n_cols);
        for (int i_row {}; i_row < ckd.n_detector_rows_binned; ++i_row) {
            auto mask_row(ckd.pixel_mask(Eigen::seqN(i_row * n_cols, n_cols)));
            auto signal_row(
              l1_prod.signal(i_alt, Eigen::seqN(i_row * n_cols, n_cols)));
            auto noise_row(
              l1_prod.noise(i_alt, Eigen::seqN(i_row * n_cols, n_cols)));
            int n_good {};
            for (int i_col {}; i_col < n_cols; ++i_col) {
                if (!mask_row(i_col)) {
                    x_values(n_good) = i_col;
                    y_values_signal(n_good) = signal_row(i_col);
                    y_values_noise(n_good) = noise_row(i_col);
                    ++n_good;
                }
            }
            const CubicSpline signal_spline { x_values.head(n_good),
                                              y_values_signal.head(n_good) };
            signal_row = signal_spline.eval(col_indices);
            const CubicSpline noise_spline { x_values.head(n_good),
                                             y_values_noise.head(n_good) };
            noise_row = noise_spline.eval(col_indices);
        }
    }
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
    Eigen::ArrayXcd signal_fft(ckd.stray.kernel_fft_sizes.maxCoeff());
    // Result of convolution
    ArrayXXd conv_result(ckd.n_detector_rows, ckd.n_detector_cols);
#pragma omp parallel for firstprivate(signal_fft, conv_result)
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        // Unbin the detector image
        const ArrayXXd signal_unbin =
          binning_table.unbin(l1_prod.signal.row(i_alt))
            .reshaped<Eigen::RowMajor>(ckd.n_detector_rows,
                                       ckd.n_detector_cols);
        // "Ideal" signal, i.e. the one without stray light
        ArrayXXd signal_ideal = signal_unbin;
        // Van Cittert algorithm
        for (int i_vc {}; i_vc < n_van_cittert; ++i_vc) {
            if (ckd.stray.n_kernels > 1) {
                convolveMulti(ckd.stray.edges,
                              ckd.stray.weights,
                              signal_ideal,
                              ckd.stray.kernel_rows,
                              ckd.stray.kernel_cols,
                              ckd.stray.kernels_fft,
                              signal_fft,
                              conv_result);
            } else {
                convolve(signal_ideal,
                         ckd.stray.kernel_rows(0),
                         ckd.stray.kernel_cols(0),
                         ckd.stray.kernels_fft.front(),
                         signal_fft,
                         conv_result);
            }
            signal_ideal = (signal_unbin - conv_result) / (1 - ckd.stray.eta);
        }
        l1_prod.signal.row(i_alt) = binning_table.bin(signal_ideal);
    }
}

auto mapFromDetector(const CKD& ckd,
                     const BinningTable& binning_table,
                     const int b_spline_order,
                     L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::swath;
    // Size of spectra per ALT position
    l1_prod.spectra.resize(l1_prod.n_alt * ckd.n_act,
                           ckd.wave.wavelengths.size());
    l1_prod.spectra_noise.resize(l1_prod.n_alt * ckd.n_act,
                                 ckd.wave.wavelengths.size());
    // Assume there is binning only across rows. From that determine
    // the number of rows of the binned detector image.
    Eigen::ArrayXd rows = Eigen::ArrayXd::Zero(ckd.n_detector_rows_binned);
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
    const Eigen::ArrayXd cols = Eigen::ArrayXd::LinSpaced(
      ckd.n_detector_cols, 0.0, ckd.n_detector_cols - 1);
    BSpline2D bspline_2d_signal {
        b_spline_order, rows, cols, ckd.swath.row_map, ckd.swath.col_map
    };
    BSpline2D bspline_2d_noise { bspline_2d_signal };
#pragma omp parallel for firstprivate(bspline_2d_signal, bspline_2d_noise)
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        bspline_2d_signal.genControlPoints(l1_prod.signal.row(i_alt));
        bspline_2d_signal.eval(l1_prod.spectra(
          Eigen::seqN(i_alt * ckd.n_act, ckd.n_act), Eigen::all));
        bspline_2d_noise.genControlPoints(l1_prod.noise.row(i_alt));
        bspline_2d_noise.eval(l1_prod.spectra_noise(
          Eigen::seqN(i_alt * ckd.n_act, ckd.n_act), Eigen::all));
    }
    l1_prod.signal.resize(0, 0);
    l1_prod.noise.resize(0, 0);
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
        l1_prod.spectra(Eigen::seqN(i_alt * ckd.n_act, ckd.n_act),
                        Eigen::all) *= ckd.rad.rad * exposure_time_inv;
        l1_prod.spectra_noise(Eigen::seqN(i_alt * ckd.n_act, ckd.n_act),
                              Eigen::all) *= ckd.rad.rad * exposure_time_inv;
    }
}

// Bin array across rows and/or oclumns
static auto binRowCol(const int row_bin,
                      const int col_bin,
                      const ArrayXXd& array) -> ArrayXXd
{
    ArrayXXd array_binned =
      ArrayXXd::Zero(array.rows() / row_bin, array.cols() / col_bin);
    for (int r {}; r < row_bin; ++r) {
        for (int c {}; c < col_bin; ++c) {
            array_binned += array(Eigen::seq(r, Eigen::last, row_bin),
                                  Eigen::seq(c, Eigen::last, col_bin));
        }
    }
    return array_binned /= row_bin * col_bin;
}

auto binL1B(const int bin, L1& l1_prod) -> void
{
    l1_prod.spectra = binRowCol(bin, 1, l1_prod.spectra);
    l1_prod.spectra_noise = binRowCol(bin, 1, l1_prod.spectra_noise);
    l1_prod.spectra_noise /= std::sqrt(static_cast<double>(bin));
    l1_prod.geo.lat = binRowCol(1, bin, l1_prod.geo.lat);
    l1_prod.geo.lon = binRowCol(1, bin, l1_prod.geo.lon);
    l1_prod.geo.height = binRowCol(1, bin, l1_prod.geo.height);
    l1_prod.geo.sza = binRowCol(1, bin, l1_prod.geo.sza);
    l1_prod.geo.saa = binRowCol(1, bin, l1_prod.geo.saa);
    l1_prod.geo.vza = binRowCol(1, bin, l1_prod.geo.vza);
    l1_prod.geo.vaa = binRowCol(1, bin, l1_prod.geo.vaa);
}

} // namespace tango
