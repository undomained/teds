// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for constructing and convolving with an ISRF

#pragma once

#include "eigen.h"

namespace tango {

class ISRF
{
private:
    // x and y-values of the ISRF data as read from file or from
    // generalized Gaussian parameters.
    Eigen::ArrayXd wavelength_diffs {};
    ArrayXXd isrf_data {};
    // Wavelength grid of unconvolved data
    Eigen::ArrayXd wavelengths_in {};
    // Wavelength grid of convolved data
    Eigen::ArrayXd wavelengths_out {};
    // Number of kernel values to be precomputed
    int n_vals_half {};
    int n_vals {};
    // Step size of wavelengths_in
    double step {};
    // Difference of first and last value of wavelength_in
    double wavelength_in_range {};
    // Precomputed kernel values
    ArrayXXd kernels {};
    // Precomputed kernel derivative values (limited to just one kernel)
    Eigen::VectorXd kernel_der {};

    // Find index in wavelengths_in corresponding to a wavelength
    [[nodiscard]] auto lookupIdx(const double lambda) const -> int;

    // Core of the convolution algorithm. The kernel can be the normal
    // kernel or the ISRF derivative.
    [[nodiscard]] auto convolve(
      const Eigen::Ref<const Eigen::VectorXd> kernel,
      const Eigen::Ref<const Eigen::VectorXd> data_in) const -> Eigen::ArrayXd;

public:
    ISRF() = default;

    // Parameters
    // ----------
    // wave_cutoff
    //     reduce the computational cost by considering the limited
    //     extent of the ISRF. Assuming fixed-step input and output
    //     wavelength grids.
    ISRF(const Eigen::ArrayXd& wavelength_diffs,
         const ArrayXXd& isrf_data,
         const Eigen::ArrayXd& wavelengths_in,
         const Eigen::ArrayXd& wavelengths_out,
         const double wave_cutoff = 0.58);

    // Interpolate input data onto the wavelength range that is
    // relevant for convolution. This is recallable part of the
    // constructor in case we need to regenerate the ISRF because of
    // spectral shift.
    auto regenerate(const double lambda0) -> void;

    // Read ISRF data from file. First ALT bin and the number of ALT
    // bins are only relevant for reading the heterogeneous ISRF.
    auto fromFile(const std::string& filename,
                  const Eigen::ArrayXd& wavelengths_in,
                  const Eigen::ArrayXd& wavelengths_out,
                  const size_t alt_beg = 0,
                  const size_t n_alt = 0) -> void;

    // Construct ISRF from generalized Gaussian parameters
    auto fromGauss(const Eigen::ArrayXd& wavelengths_in,
                   const Eigen::ArrayXd& wavelengths_out,
                   const double fwhm,
                   const double shape) -> void;

    // Convolve the ith kernel with data. For a simple ISRF i_kernel = 0.
    [[nodiscard]] auto convolve(const Eigen::Ref<const Eigen::VectorXd> data_in,
                                const int i_kernel = 0) const -> Eigen::ArrayXd;

    // Convolve with ISRF derivative with respect to lambda0
    [[nodiscard]] auto convolveDer(
      const Eigen::Ref<const Eigen::VectorXd> data_in) const -> Eigen::ArrayXd;
};

} // namespace tango
