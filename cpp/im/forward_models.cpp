// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "forward_models.h"

#include <common/b_spline_2d.h>
#include <common/binning_table.h>
#include <common/ckd.h>
#include <common/fourier.h>
#include <common/isrf.h>
#include <common/l1.h>
#include <netcdf>
#include <random>
#include <spdlog/spdlog.h>

namespace tango {

auto applyISRF(const CKD& ckd,
               const bool enabled,
               const ISRF& isrf,
               const std::string& sgm_filename,
               L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::l1b;
    // If this process is disabled then linearly interpolate the
    // line-by-line spectra onto the CKD wavelength grids. We cannot
    // simply return like the other processes.
    const size_t n_waves_in { static_cast<size_t>(l1_prod.wavelengths.size()) };
    const size_t n_waves_out { static_cast<size_t>(
      ckd.wave.wavelengths.size()) };
    ArrayXXd spectra_out(l1_prod.n_alt * ckd.n_act, n_waves_out);
    if (enabled) {
        // If the spectra are not in memory then read them one by one
        // for the convolution.
        if (l1_prod.spectra.size() == 0) {
            const netCDF::NcFile nc { sgm_filename, netCDF::NcFile::read };
            const netCDF::NcVar nc_var { nc.getVar("radiance") };
            ArrayXXd buf(ckd.n_act, n_waves_in);
            for (int i_alt {}; i_alt < l1_prod.n_alt; ++i_alt) {
                nc_var.getVar({ static_cast<size_t>(i_alt), 0, 0 },
                              { 1, static_cast<size_t>(ckd.n_act), n_waves_in },
                              buf.data());
#pragma omp parallel for
                for (int i_act = 0; i_act < ckd.n_act; ++i_act) {
                    const int act_idx { i_alt * ckd.n_act + i_act };
                    spectra_out.row(act_idx) =
                      isrf.convolve(buf.row(i_act), act_idx);
                }
            }
        } else {
#pragma omp parallel for
            for (int i = 0; i < l1_prod.n_alt * ckd.n_act; ++i) {
                spectra_out.row(i) = isrf.convolve(l1_prod.spectra.row(i), i);
            }
        }
    } else {
#pragma omp parallel for
        for (int i_act = 0; i_act < l1_prod.n_alt * ckd.n_act; ++i_act) {
            const LinearSpline spline { l1_prod.wavelengths,
                                        l1_prod.spectra.row(i_act) };
            for (size_t i_wave {}; i_wave < n_waves_out; ++i_wave) {
                spectra_out.row(i_act)(i_wave) =
                  spline.eval(ckd.wave.wavelengths(i_wave));
            }
        }
    }
    l1_prod.spectra = std::move(spectra_out);
    l1_prod.wavelengths = ckd.wave.wavelengths;
}

auto radiometric(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::swath;
    if (!enabled) {
        return;
    }
    const double exposure_time_inv { 1.0 / l1_prod.exposure_time };
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        l1_prod.spectra(Eigen::seqN(i_alt * ckd.n_act, ckd.n_act),
                        Eigen::all) /= ckd.rad.rad * exposure_time_inv;
    }
}

auto mapToDetector(const CKD& ckd,
                   const int b_spline_order,
                   L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::stray;
    l1_prod.signal.resize(l1_prod.n_alt, ckd.npix);
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        const BSpline2D bspline_2d {
            b_spline_order,
            ckd.swath.act_angles,
            ckd.wave.wavelengths,
            l1_prod.spectra(Eigen::seqN(i_alt * ckd.n_act, ckd.n_act),
                            Eigen::all)
        };
        bspline_2d.eval(ckd.swath.act_map,
                        ckd.swath.wavelength_map,
                        l1_prod.signal.row(i_alt));
        for (int i {}; i < ckd.npix; ++i) {
            l1_prod.signal(i_alt, i) = std::max(0.0, l1_prod.signal(i_alt, i));
        }
    }
    l1_prod.spectra = ArrayXXd {};
}

auto strayLight(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::prnu;
    if (!enabled) {
        return;
    }
    // Pre-allocate the signal FFT for efficiency
    Eigen::ArrayXcd signal_fft(ckd.stray.kernel_fft_sizes.maxCoeff());
    // Result of convolution
    ArrayXXd conv_result(ckd.n_detector_rows, ckd.n_detector_cols);
#pragma omp parallel for firstprivate(signal_fft, conv_result)
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        auto signal = l1_prod.signal.row(i_alt).reshaped<Eigen::RowMajor>(
          ckd.n_detector_rows, ckd.n_detector_cols);
        if (ckd.stray.n_kernels > 1) {
            convolveMulti(ckd.stray.edges,
                          ckd.stray.weights,
                          signal,
                          ckd.stray.kernel_rows,
                          ckd.stray.kernel_cols,
                          ckd.stray.kernels_fft,
                          signal_fft,
                          conv_result);
        } else {
            convolve(signal,
                     ckd.stray.kernel_rows(0),
                     ckd.stray.kernel_cols(0),
                     ckd.stray.kernels_fft.front(),
                     signal_fft,
                     conv_result);
        }
        signal = (1.0 - ckd.stray.eta) * signal + conv_result;
        l1_prod.signal.row(i_alt) = signal.reshaped<Eigen::RowMajor>();
    }
}

auto prnu(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::nonlin;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        l1_prod.signal.row(i_alt) *= ckd.prnu.prnu;
    }
}

auto nonlinearity(const bool enabled,
                  const LinearSpline& nonlin_spline,
                  L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::dark_current;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        l1_prod.signal.row(i_alt) =
          nonlin_spline.eval(l1_prod.signal.row(i_alt));
    }
}

auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::noise;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        l1_prod.signal.row(i_alt) += ckd.dark.current * l1_prod.exposure_time;
    }
}

auto noise(const CKD& ckd,
           const bool enabled,
           const int seed,
           const int n_coadditions,
           L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::dark_offset;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        static std::mt19937 gen { static_cast<unsigned long>(seed) };
        for (int i {}; i < ckd.npix; ++i) {
            const double noise_value { std::sqrt(
              (ckd.noise.n2(i)
               + std::abs(l1_prod.signal(i_alt, i)) * ckd.noise.g(i))
              / n_coadditions) };
            std::normal_distribution<> d { 0.0, noise_value };
            l1_prod.signal(i_alt, i) += d(gen);
        }
    }
}

auto darkOffset(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::raw;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        l1_prod.signal.row(i_alt) += ckd.dark.offset;
    }
}

auto binDetectorImages(const int n_rows,
                       const int n_cols,
                       const std::string& binning_filename,
                       const int binning_table_id,
                       const bool scale_by_binsize,
                       L1& l1_prod) -> void
{
    const BinningTable binning_table {
        n_rows, n_cols, binning_filename, binning_table_id
    };
    l1_prod.signal = binning_table.binMulti(l1_prod.signal, scale_by_binsize);
    l1_prod.binning_table_id = static_cast<uint8_t>(binning_table_id);
}

auto digitalToAnalog(const int nr_coadditions, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::l1a;
    l1_prod.nr_coadditions = static_cast<uint16_t>(nr_coadditions);
    l1_prod.signal *= l1_prod.nr_coadditions;
    l1_prod.signal = Eigen::round(l1_prod.signal);
}

auto estimateOptimalCoadd(const CKD& ckd,
                          const int FMC,
                          const double exposure_time,
                          const double f_sat,
                          const double full_well,
                          const double t_dead,
                          const L1& l1_prod) -> void
{
    constexpr double dwell { 24.0 };
    const double t_dwell { FMC / dwell };
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        double max_val {};
        for (int i {}; i < ckd.npix; ++i) {
            if (!ckd.pixel_mask(i)) {
                max_val = std::max(max_val, l1_prod.signal(i_alt, i));
            }
        }
        const double I_sig { max_val / exposure_time };
        const double coadd_raw { I_sig * t_dwell
                                 / (f_sat * full_well + I_sig * t_dead) };
        const int n_coadd { static_cast<int>(std::ceil(coadd_raw)) };
        spdlog::info("  {:5.2f} ({}) {:6.4} ms",
                     coadd_raw,
                     n_coadd,
                     1e3 * (t_dwell / n_coadd) - t_dead);
    }
}

} // namespace tango
