// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for constructing an ISRF and convolving data with the ISRF.

#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace tango {

class ISRF
{
private:
    // Wavelength grid of unconvolved data
    std::vector<double> wavelengths_in {};
    // Wavelength grid of convolved data
    std::vector<double> wavelengths_out {};
    // Number ofkernel values to be precomputed
    int n_vals_half {};
    int n_vals {};
    // Step size of wavelengths_in
    double step {};
    // Difference of first and last value of wavelength_in
    double wavelength_in_range {};
    // Precomputed kernel values
    std::vector<Eigen::VectorXd> kernels {};

    // Find index in wavelengths_in corresponding to a wavelength
    auto lookupIdx(const double lambda) const -> size_t;

public:
    ISRF() = default;

    // Parameters
    // ----------
    // wave_cutoff
    //     reduce the computational cost by considering the limited
    //     extent of the ISRF. Assuming fixed-step input and output
    //     wavelength grids.
    ISRF(const std::vector<double>& wavelength_diffs,
         const std::vector<double>& isrf_data,
         const std::vector<double>& wavelengths_in,
         const std::vector<double>& wavelengths_out,
         const double wave_cutoff = 0.7);

    // Read ISRF data from file
    auto fromFile(const std::string& filename,
                  const std::vector<double>& wavelengths_in,
                  const std::vector<double>& wavelengths_out) -> void;

    // Construct ISRF from generalized Gaussian parameters
    auto fromGauss(const std::vector<double>& wavelengths_in,
                   const std::vector<double>& wavelengths_out,
                   const double fwhm,
                   const double shape) -> void;

    // Convolve ith kernel with data. For a simple ISRF i_kernel = 0
    // always.
    auto convolve(const size_t i_kernel,
                  double* const data_in,
                  double* data_out) const -> void;
};

} // namespace tango
