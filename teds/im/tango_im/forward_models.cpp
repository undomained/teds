// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "forward_models.h"

#include <Eigen/Sparse>
#include <algorithm>
#include <random>
#include <tango_l1b/b_spline_2d.h>
#include <tango_l1b/binning_table.h>
#include <tango_l1b/ckd.h>
#include <tango_l1b/fourier.h>
#include <tango_l1b/l1.h>

namespace tango {

// Generate an ISRF kernel for convolving the radiances from the
// line-by-line grid onto the intermediate grids from the CKD swath
// section.
static auto genISRFKernel(const double fwhm_gauss,
                          const double shape,
                          const std::vector<double>& lbl_wavelengths,
                          const std::vector<double>& wavelengths,
                          Eigen::SparseMatrix<double>& isrf) -> void
{
    isrf =
      Eigen::SparseMatrix<double>(wavelengths.size(), lbl_wavelengths.size());
    const bool do_gaussian { std::abs(shape - 2.0) < 1e-30 };
    const double sigma { fwhm_gauss / (2.0 * std::sqrt(2.0 * std::log(2.0))) };
    const double sigma_2inv { 1.0 / (2.0 * sigma * sigma) };
    std::vector<Eigen::Triplet<double>> triplets {};
    std::vector<double> row_sums(isrf.rows(), 0.0); // For normalization
    for (int i_wave {}; i_wave < isrf.rows(); ++i_wave) {
        for (int i_lbl {}; i_lbl < isrf.cols(); ++i_lbl) {
            // Convolution result for this CKD wavelength
            double conv { wavelengths[i_wave] - lbl_wavelengths[i_lbl] };
            // Reduce the computational cost by considering the
            // limited extent of the Gaussian.
            constexpr double wave_rel_threshold { 1.5 };
            if (std::abs(conv) > wave_rel_threshold * fwhm_gauss) {
                continue;
            }
            if (do_gaussian) {
                conv = std::exp(-(conv * conv * sigma_2inv));
            } else {
                conv = std::pow(
                  2.0, -std::pow(2 * std::abs(conv) / fwhm_gauss, shape));
            }
            row_sums[i_wave] += conv;
            triplets.push_back({ i_wave, i_lbl, conv });
        }
    }
    isrf.setFromTriplets(triplets.begin(), triplets.end());
    // Normalize
    for (int k {}; k < isrf.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it { isrf, k }; it;
             ++it) {
            it.valueRef() /= row_sums[it.row()];
        }
    }
}

auto applyISRF(const CKD& ckd,
               const bool enabled,
               const double fwhm_gauss,
               const double shape,
               L1& l1_prod) -> void
{
    l1_prod.level = ProcLevel::l1b;
    // If this process is disabled then linearly interpolate the
    // line-by-line spectra onto the CKD wavelength grids. We cannot
    // simply return like the other processes.
    const int n_waves_in { static_cast<int>(
      l1_prod.wavelengths.front().size()) };
    const int n_waves_out { static_cast<int>(ckd.swath.wavelengths.size()) };
    std::vector<double> spectra_out(l1_prod.n_alt * ckd.n_act * n_waves_out);
    if (enabled) {
        Eigen::SparseMatrix<double> isrf {};
        genISRFKernel(fwhm_gauss,
                      shape,
                      l1_prod.wavelengths.front(),
                      ckd.swath.wavelengths,
                      isrf);
#pragma omp parallel for
        for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
            Eigen::VectorXd result(ckd.n_detector_cols);
            for (int i_act {}; i_act < ckd.n_act; ++i_act) {
                const int act_idx { i_alt * ckd.n_act + i_act };
                result = isrf
                         * Eigen::Map<Eigen::VectorXd>(
                           &l1_prod.spectra[act_idx * n_waves_in], n_waves_in);
                for (int i {}; i < n_waves_out; ++i) {
                    spectra_out[act_idx * n_waves_out + i] = result[i];
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int i_alt = 0; i_alt < l1_prod.n_alt; ++i_alt) {
            for (int i_act {}; i_act < ckd.n_act; ++i_act) {
                const int act_idx { i_alt * ckd.n_act + i_act };
                const std::vector<double> spectrum {
                    l1_prod.spectra.begin() + act_idx * n_waves_in,
                    l1_prod.spectra.begin() + (act_idx + 1) * n_waves_in
                };
                const LinearSpline spline { l1_prod.wavelengths[i_act],
                                            spectrum };
                for (int i {}; i < n_waves_out; ++i) {
                    spectra_out[act_idx * n_waves_out + i] =
                      spline.eval(ckd.swath.wavelengths[i]);
                }
            }
        }
    }
    l1_prod.spectra = std::move(spectra_out);
    for (int i_act {}; i_act < ckd.n_act; ++i_act) {
        l1_prod.wavelengths[i_act] = ckd.swath.wavelengths;
    }
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
        for (int i_act {}; i_act < ckd.n_act; ++i_act) {
            for (int i {}; i < static_cast<int>(ckd.swath.wavelengths.size());
                 ++i) {
                l1_prod.spectra[(i_alt * ckd.n_act + i_act)
                                  * ckd.swath.wavelengths.size()
                                + i] /=
                  ckd.rad.rad[i_act][0] * exposure_time_inv;
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
            ckd.swath.wavelengths,
            &l1_prod.spectra[i_alt * ckd.n_act * ckd.swath.wavelengths.size()]
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
            const double noise_value { std::sqrt(ckd.noise.n2[i]
                                                 + std::abs(l1_prod.signal[idx])
                                                     * ckd.noise.g[i]) };
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
    if (scale_by_binsize) {
        binning_table.bin(l1_prod.signal);
    } else {
        binning_table.binUnscaled(l1_prod.signal);
    }
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

} // namespace tango
