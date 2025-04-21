// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "isrf.h"

#include "cubic_spline.h"

#include <cmath>
#include <netcdf>
#include <numeric>
#include <spdlog/spdlog.h>

namespace tango {

ISRF::ISRF(const std::vector<double>& wavelength_diffs,
           const std::vector<double>& isrf_data,
           const std::vector<double>& wavelengths_in,
           const std::vector<double>& wavelengths_out,
           const double wave_cutoff)
  : wavelengths_in { wavelengths_in }, wavelengths_out { wavelengths_out }
{
    step = wavelengths_in[1] - wavelengths_in[0];
    wavelength_in_range = wavelengths_in.back() - wavelengths_in.front();
    // This many values of the kernel need to be precomputed and
    // stored. Outside the range the kernel is 0.
    n_vals_half = static_cast<int>(wave_cutoff / step);
    n_vals = 2 * n_vals_half - 1;
    // Array for storing interpolated convolution values. There is one
    // array for wavelengths but possibly many sets of ISRF data
    // corresponding to different configurations. For instance, we
    // could have a different ISRF for each ALT/ACT position. This
    // class, however, does not explicitly know about the geometry. It
    // only stores the kernels in order that they are read.
    const auto n_kernels { isrf_data.size() / wavelength_diffs.size() };
    kernels.resize(n_kernels);
    for (int i_kernel {}; i_kernel < static_cast<int>(n_kernels); ++i_kernel) {
        auto& kernel { kernels[i_kernel] };
        kernel = Eigen::VectorXd(n_vals);
        const CubicSpline kernel_spline {
            wavelength_diffs,
            { isrf_data.cbegin() + i_kernel * wavelength_diffs.size(),
              isrf_data.cbegin() + (i_kernel + 1) * wavelength_diffs.size() }
        };
        for (int i {}; i < n_vals; ++i) {
            kernel[i] = kernel_spline.eval(-n_vals_half * step + i * step);
        }
        kernel /= kernel.sum();
    }
}

auto ISRF::fromFile(const std::string& filename,
                    const std::vector<double>& wavelengths_in,
                    const std::vector<double>& wavelengths_out) -> void

{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    std::vector<double> wavelengths {};
    std::vector<double> isrf {};
    if (const auto grp { nc.getGroup("isrf") }; grp.isNull()) {
        spdlog::info("  Reading homogeneous ISRF");
        const auto n_wavelengths { nc.getDim("wavelength").getSize() };
        wavelengths.resize(n_wavelengths);
        isrf.resize(n_wavelengths);
        nc.getVar("wavelength").getVar(wavelengths.data());
        nc.getVar("isrf").getVar(isrf.data());
    } else {
        spdlog::info("  Reading heterogenous ISRF");
        const auto n_alt { nc.getDim("along_track_sample").getSize() };
        const auto n_act { nc.getDim("across_track_sample").getSize() };
        const auto n_wavelengths { grp.getDim("wavelength").getSize() };
        wavelengths.resize(n_wavelengths);
        isrf.resize(n_alt * n_act * n_wavelengths);
        grp.getVar("wavelength").getVar(wavelengths.data());
        grp.getVar("isrf").getVar(isrf.data());
    }
    *this = ISRF(wavelengths, isrf, wavelengths_in, wavelengths_out);
}

auto ISRF::fromGauss(const std::vector<double>& wavelengths_in,
                     const std::vector<double>& wavelengths_out,
                     const double fwhm,
                     const double shape) -> void
{
    constexpr double wave_cutoff { 0.7 }; // nm
    const double step { (wavelengths_in[1] - wavelengths_in[0]) };
    const int n_vals_half { static_cast<int>(wave_cutoff / step) };
    const int n_vals { 2 * n_vals_half + 1 };
    std::vector<double> wavelength_diffs(n_vals);
    std::vector<double> isrf(n_vals);
    for (int i {}; i < n_vals; ++i) {
        wavelength_diffs[i] = -n_vals_half * step + i * step;
        isrf[i] = std::pow(
          2.0, -std::pow(2 * std::abs(wavelength_diffs[i]) / fwhm, shape));
    }
    *this = ISRF(wavelength_diffs, isrf, wavelengths_in, wavelengths_out);
}

auto ISRF::lookupIdx(const double lambda) const -> size_t
{
    return static_cast<size_t>((lambda - wavelengths_in.front())
                               / wavelength_in_range
                               * (wavelengths_in.size() - 1));
}

auto ISRF::convolve(const size_t i_kernel,
                    double* const data_in,
                    double* data_out) const -> void
{
    const auto kernel { kernels[std::min(kernels.size() - 1, i_kernel)] };
    for (int i_conv {}; i_conv < static_cast<int>(wavelengths_out.size());
         ++i_conv) {
        const double first_wavelength { wavelengths_out[i_conv]
                                        - n_vals_half * step };
        const auto first_idx { lookupIdx(first_wavelength) };
        data_out[i_conv] = kernel.dot(
          Eigen::Map<Eigen::VectorXd>(data_in + first_idx, kernel.size()));
    }
}

} // namespace tango
