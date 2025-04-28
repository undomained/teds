// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Functions related to fast Fourier transforms

#pragma once

#include "constants.h"
#include "eigen.h"

namespace tango {

// Return the optimal size in units of std::complex<double> which to
// store the Fourier transform.
auto getFFTSize(const int n_rows, const int n_cols) -> int;

// Fourier transform of a 2D real array (real-to-complex)
auto fft_r2c(const ArrayXXd& image, Eigen::ArrayXcd& image_fft) -> void;

// Inverse Fourier transform of a 2D complex array that stores roughly
// half the elements (complex-to-real)
auto fft_c2r(Eigen::ArrayXcd& image_fft, ArrayXXd& image) -> void;

// Convolve image_in with kernel_fft and store the result in image_out.
//   image_in - image in real space
//   kernel_n_rows - number of rows of the real space kernel
//   kernel_n_cols - number of columns of the real space kernel
//   kernel_fft - Fourier transform of the kernel
//   kernel_image - preallocated buffer for taking the image FFT
//   image_out - convolution result in real space
// Note that while kernel_n_rows and kernel_n_cols are required as
// arguments the real space kernel itself is not used, only its
// Fourier transform. Importantly, the size of kernel_fft is not
// kernel_n_rows * kernel_n_cols but is determined by getFFTSize.
auto convolve(const Eigen::Ref<const ArrayXXd> image_in,
              const int kernel_n_rows,
              const int kernel_n_cols,
              const Eigen::ArrayXcd& kernel_fft,
              Eigen::ArrayXcd& image_fft,
              ArrayXXd& image_out) -> void;

// Convolution with multiple kernels, each acting on a different part
// of the image.
auto convolveMulti(
  const Eigen::Array<int, Eigen::Dynamic, box::n, Eigen::RowMajor>& edges,
  const ArrayXXd& weights,
  const Eigen::Ref<const ArrayXXd> image_in,
  const Eigen::VectorXi& kernel_rows,
  const Eigen::VectorXi& kernel_cols,
  const std::vector<Eigen::ArrayXcd>& kernels_fft,
  Eigen::ArrayXcd& image_fft,
  ArrayXXd& image_out) -> void;

} // namespace tango
