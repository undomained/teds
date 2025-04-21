// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "forward_models.h"

#include <common/isrf.h>

#include <Eigen/Sparse>
#include <algorithm>
#include <common/b_spline_2d.h>
#include <common/binning_table.h>
#include <common/ckd.h>
#include <common/cubic_spline.h>
#include <common/fourier.h>
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
    const size_t n_waves_in { l1_prod.wavelengths.size() };
    const int n_waves_out { static_cast<int>(ckd.wave.wavelengths.size()) };
    std::vector<double> spectra_out(l1_prod.n_alt * ckd.n_act * n_waves_out);
    if (enabled) {
        // If the spectra are not in memory then read them one by one
        // for the convolution.
        if (l1_prod.spectra.empty()) {
            const netCDF::NcFile nc { sgm_filename, netCDF::NcFile::read };
            const netCDF::NcVar nc_var { nc.getVar("radiance") };
            std::vector<double> buf(ckd.n_act * n_waves_in);
            const size_t n_act { static_cast<size_t>(ckd.n_act) };
            for (size_t i_alt {}; i_alt < static_cast<size_t>(l1_prod.n_alt);
                 ++i_alt) {
                nc_var.getVar(
                  { i_alt, 0, 0 }, { 1, n_act, n_waves_in }, buf.data());
#pragma omp parallel for
                for (int i_act = 0; i_act < ckd.n_act; ++i_act) {
                    const size_t act_idx { i_alt * ckd.n_act + i_act };
                    isrf.convolve(act_idx,
                                  &buf[i_act * n_waves_in],
                                  &spectra_out[act_idx * n_waves_out]);
                }
            }
        } else {
#pragma omp parallel for
            for (size_t i_alt = 0; i_alt < static_cast<size_t>(l1_prod.n_alt);
                 ++i_alt) {
                for (int i_act {}; i_act < ckd.n_act; ++i_act) {
                    const size_t act_idx { i_alt * ckd.n_act + i_act };
                    isrf.convolve(act_idx,
                                  &l1_prod.spectra[act_idx * n_waves_in],
                                  &spectra_out[act_idx * n_waves_out]);
                }
            }
        }
    } else {
#pragma omp parallel for
        for (size_t i_alt = 0; i_alt < static_cast<size_t>(l1_prod.n_alt);
             ++i_alt) {
            for (int i_act {}; i_act < ckd.n_act; ++i_act) {
                const size_t act_idx { i_alt * ckd.n_act + i_act };
                const std::vector<double> spectrum {
                    l1_prod.spectra.begin() + act_idx * n_waves_in,
                    l1_prod.spectra.begin() + (act_idx + 1) * n_waves_in
                };
                const LinearSpline spline { l1_prod.wavelengths, spectrum };
                for (int i {}; i < n_waves_out; ++i) {
                    spectra_out[act_idx * n_waves_out + i] =
                      spline.eval(ckd.wave.wavelengths[i]);
                }
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
    for (size_t i_alt = 0; i_alt < static_cast<size_t>(l1_prod.n_alt);
         ++i_alt) {
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            for (int i {}; i < ckd.n_wavelengths; ++i) {
                l1_prod.spectra[(i_alt * ckd.n_act + i_act) * ckd.n_wavelengths
                                + i] /=
                  ckd.rad.rad[i_act][i] * exposure_time_inv;
            }
        }
    }
}

auto mapToDetector(const CKD& ckd,
                   const int b_spline_order,
                   L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::stray;
    l1_prod.signal.assign(l1_prod.n_alt * ckd.npix, 0.0);
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        const BSpline2D bspline_2d {
            b_spline_order,
            ckd.swath.act_angles,
            ckd.wave.wavelengths,
            &l1_prod.spectra[i_alt * ckd.n_act * ckd.wave.wavelengths.size()]
        };
        bspline_2d.eval(ckd.swath.act_map,
                        ckd.swath.wavelength_map,
                        &l1_prod.signal[i_alt * ckd.npix]);
        for (int i {}; i < ckd.npix; ++i) {
            l1_prod.signal[i_alt * ckd.npix + i] =
              std::max(0.0, l1_prod.signal[i_alt * ckd.npix + i]);
        }
    }
    l1_prod.spectra = std::vector<double> {};
}

auto strayLight(const CKD& ckd, const bool enabled, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::prnu;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        std::vector<double> conv_result(ckd.npix, 0.0);
        std::vector<std::complex<double>> image_fft(
          *std::ranges::max_element(ckd.stray.kernel_fft_sizes));
        for (int i_kernel {}; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
            std::vector<double> image_weighted {
                l1_prod.signal.begin() + i_alt * ckd.npix,
                l1_prod.signal.begin() + (i_alt + 1) * ckd.npix
            };
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
                    conv_result
                      [(i + ckd.stray.edges[i_kernel * box::n + box::b])
                         * ckd.n_detector_cols
                       + j + ckd.stray.edges[i_kernel * box::n + box::l]] +=
                      conv_result_sub[i * image_n_cols + j];
                }
            }
        }
        for (int i {}; i < ckd.npix; ++i) {
            const int idx { i_alt * ckd.npix + i };
            l1_prod.signal[idx] =
              (1.0 - ckd.stray.eta[i]) * l1_prod.signal[idx] + conv_result[i];
        }
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
        for (int i {}; i < ckd.npix; ++i) {
            if (!ckd.pixel_mask[i]) {
                l1_prod.signal[i_alt * ckd.npix + i] *= ckd.prnu.prnu[i];
            }
        }
    }
}

auto nonlinearity(const CKD& ckd,
                  const bool enabled,
                  const LinearSpline& nonlin_spline,
                  L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::dark_current;
    if (!enabled) {
        return;
    }
#pragma omp parallel for
    for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
        for (int i {}; i < ckd.npix; ++i) {
            if (!ckd.pixel_mask[i]) {
                const int idx { i_alt * ckd.npix + i };
                l1_prod.signal[idx] = nonlin_spline.eval(l1_prod.signal[idx]);
            }
        }
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
        for (int i {}; i < ckd.npix; ++i) {
            if (!ckd.pixel_mask[i]) {
                l1_prod.signal[i_alt * ckd.npix + i] +=
                  ckd.dark.current[i] * l1_prod.exposure_time;
            }
        }
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
            const int idx { i_alt * ckd.npix + i };
            const double noise_value { std::sqrt(
              (ckd.noise.n2[i] + std::abs(l1_prod.signal[idx]) * ckd.noise.g[i])
              / n_coadditions) };
            std::normal_distribution<> d { 0.0, noise_value };
            l1_prod.signal[idx] += d(gen);
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
        for (int i {}; i < ckd.npix; ++i) {
            if (!ckd.pixel_mask[i]) {
                l1_prod.signal[i_alt * ckd.npix + i] += ckd.dark.offset[i];
            }
        }
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
    binning_table.bin(l1_prod.signal, scale_by_binsize);
    l1_prod.binning_table_id = static_cast<uint8_t>(binning_table_id);
}

auto digitalToAnalog(const int nr_coadditions, L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::l1a;
    l1_prod.nr_coadditions = static_cast<uint16_t>(nr_coadditions);
    for (double& val : l1_prod.signal) {
        val *= l1_prod.nr_coadditions;
        val = std::round(val);
    }
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
            if (!ckd.pixel_mask[i]) {
                max_val =
                  std::max(max_val, l1_prod.signal[i_alt * ckd.npix + i]);
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
