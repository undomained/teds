// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for constructing and convolving with an ISRF

#pragma once

#include "eigen.h"

namespace tango {

class ISRF
{
private:
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

    // Find index in wavelengths_in corresponding to a wavelength
    auto lookupIdx(const double lambda) const -> int;

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
         const double wave_cutoff = 0.7);

    // Read ISRF data from file
    auto fromFile(const std::string& filename,
                  const Eigen::ArrayXd& wavelengths_in,
                  const Eigen::ArrayXd& wavelengths_out) -> void;

    // Construct ISRF from generalized Gaussian parameters
    auto fromGauss(const Eigen::ArrayXd& wavelengths_in,
                   const Eigen::ArrayXd& wavelengths_out,
                   const double fwhm,
                   const double shape) -> void;

    // Convolve the ith kernel with data. For a simple ISRF i_kernel = 0.
    auto convolve(const Eigen::Ref<const Eigen::VectorXd> data_in,
                  const int i_kernel = 0) const -> Eigen::ArrayXd;
};

} // namespace tango
