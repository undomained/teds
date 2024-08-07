// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Functions related to fast Fourier transforms

#pragma once

#include <complex>
#include <vector>

namespace tango {

// Return the optimal size in units of std::complex<double> which to
// store the Fourier transform.
auto getFFTSize(const int n_rows, const int n_cols) -> int;

// Fourier transform of a 2D real array (real-to-complex)
auto fft_r2c(const int n_rows,
             const int n_cols,
             std::vector<double>& image,
             std::vector<std::complex<double>>& image_fft) -> void;

// Inverse Fourier transform of a 2D complex array that stores roughly
// half the elements (complex-to-real)
auto fft_c2r(const int n_rows,
             const int n_cols,
             std::vector<std::complex<double>>& image_fft,
             std::vector<double>& image) -> void;

// Convolve image_in with kernel_fft and store the result in image_out.
//   image_n_rows - number of image rows
//   image_n_cols - number of image columns
//   image_in - image in real space
//   kernel_n_rows - number of rows of the real space kernel
//   kernel_n_cols - number of columns of the real space kernel
//   kernel_fft - Fourier transform of the kernel
//   kernel_image - preallocated buffer for taking the image FFT
//   image_out - convolution result in real space
// Note that while kernel_n_rows and kernel_n_cols are required as
// arguments the real space kernel itself is not used, only its
// Fourier transform.  Importantly, the size of kernel_fft is not
// kernel_n_rows * kernel_n_cols but is determined by getFFTSize.
auto convolve(const int image_n_rows,
              const int image_n_cols,
              const std::vector<double>& image_in,
              const int kernel_n_rows,
              const int kernel_n_cols,
              const std::vector<std::complex<double>>& kernel_fft,
              std::vector<std::complex<double>>& image_fft,
              std::vector<double>& image_out) -> void;

} // namespace tango
