// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "isrf.h"

#include "cubic_spline.h"

#include <netcdf>
#include <spdlog/spdlog.h>

namespace tango {

ISRF::ISRF(const Eigen::ArrayXd& wavelength_diffs,
           const ArrayXXd& isrf_data,
           const Eigen::ArrayXd& wavelengths_in,
           const Eigen::ArrayXd& wavelengths_out,
           const double wave_cutoff)
  : wavelengths_in { wavelengths_in }, wavelengths_out { wavelengths_out }
{
    step = wavelengths_in(1) - wavelengths_in(0);
    wavelength_in_range =
      wavelengths_in(wavelengths_in.size() - 1) - wavelengths_in(0);
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
    kernels.resize(isrf_data.rows(), n_vals);
    for (int i_kernel {}; i_kernel < static_cast<int>(isrf_data.rows());
         ++i_kernel) {
        auto kernel = kernels.row(i_kernel);
        const CubicSpline kernel_spline { wavelength_diffs,
                                          isrf_data.row(i_kernel) };
        for (int i {}; i < n_vals; ++i) {
            kernel(i) = kernel_spline.eval(-n_vals_half * step + i * step);
        }
        kernel /= kernel.sum();
    }
}

auto ISRF::fromFile(const std::string& filename,
                    const Eigen::ArrayXd& wavelengths_in,
                    const Eigen::ArrayXd& wavelengths_out) -> void

{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    Eigen::ArrayXd wavelengths {};
    ArrayXXd isrf {};
    if (const auto grp { nc.getGroup("isrf") }; grp.isNull()) {
        spdlog::info("  Reading homogeneous ISRF");
        const auto n_wavelengths { nc.getDim("wavelength").getSize() };
        wavelengths.resize(n_wavelengths);
        isrf.resize(1, n_wavelengths);
        nc.getVar("wavelength").getVar(wavelengths.data());
        nc.getVar("isrf").getVar(isrf.data());
    } else {
        spdlog::info("  Reading heterogenous ISRF");
        const auto n_alt { nc.getDim("along_track_sample").getSize() };
        const auto n_act { nc.getDim("across_track_sample").getSize() };
        const auto n_wavelengths { grp.getDim("wavelength").getSize() };
        wavelengths.resize(n_wavelengths);
        isrf.resize(n_alt * n_act, n_wavelengths);
        grp.getVar("wavelength").getVar(wavelengths.data());
        grp.getVar("isrf").getVar(isrf.data());
    }
    *this = ISRF(wavelengths, isrf, wavelengths_in, wavelengths_out);
}

auto ISRF::fromGauss(const Eigen::ArrayXd& wavelengths_in,
                     const Eigen::ArrayXd& wavelengths_out,
                     const double fwhm,
                     const double shape) -> void
{
    constexpr double wave_cutoff { 0.7 }; // nm
    const double step { (wavelengths_in(1) - wavelengths_in(0)) };
    const int n_vals_half { static_cast<int>(wave_cutoff / step) };
    const int n_vals { 2 * n_vals_half + 1 };
    auto wavelength_diffs = Eigen::ArrayXd::LinSpaced(
      n_vals, -n_vals_half * step, n_vals_half * step);
    ArrayXXd isrf(1, n_vals);
    isrf.row(0) =
      Eigen::pow(2.0, -(2 * wavelength_diffs.abs() / fwhm).pow(shape));
    *this = ISRF(wavelength_diffs, isrf, wavelengths_in, wavelengths_out);
}

auto ISRF::lookupIdx(const double lambda) const -> int
{
    return static_cast<int>((lambda - wavelengths_in(0)) / wavelength_in_range
                            * (wavelengths_in.size() - 1));
}

auto ISRF::convolve(const Eigen::Ref<const Eigen::VectorXd> data_in,
                    const int i_kernel) const -> Eigen::ArrayXd
{
    const auto& kernel = kernels.matrix().row(
      std::min(i_kernel, static_cast<int>(kernels.rows() - 1)));
    Eigen::ArrayXd data_out(wavelengths_out.size());
    for (int i_conv {}; i_conv < static_cast<int>(data_out.size()); ++i_conv) {
        const double first_wavelength { wavelengths_out(i_conv)
                                        - n_vals_half * step };
        const auto first_idx { lookupIdx(first_wavelength) };
        data_out(i_conv) =
          kernel.dot(data_in(Eigen::seqN(first_idx, kernel.size())));
    }
    return data_out;
}

} // namespace tango
