// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "fourier.h"

#include <cassert>
#include <pocketfft.h>

namespace tango {

auto getFFTSize(const int n_rows, const int n_cols) -> int
{
    return (n_rows / 2 + 1) * n_cols;
}

auto fft_r2c(const ArrayXXd& image, Eigen::ArrayXcd& image_fft) -> void
{
    const std::vector<std::size_t> shape {
        static_cast<std::size_t>(image.rows()),
        static_cast<std::size_t>(image.cols())
    };
    const std::vector<std::ptrdiff_t> stride_in {
        static_cast<std::ptrdiff_t>(image.cols() * sizeof(double)),
        sizeof(double)
    };
    const std::vector<std::ptrdiff_t> stride_out {
        static_cast<std::ptrdiff_t>(image.cols()
                                    * sizeof(std::complex<double>)),
        sizeof(std::complex<double>)
    };
    const std::vector<std::size_t> axes { 1, 0 };
    pocketfft::r2c(shape,
                   stride_in,
                   stride_out,
                   axes,
                   pocketfft::FORWARD,
                   image.data(),
                   image_fft.data(),
                   1.0);
}

auto fft_c2r(Eigen::ArrayXcd& image_fft, ArrayXXd& image) -> void
{
    const std::vector<std::size_t> shape {
        static_cast<std::size_t>(image.rows()),
        static_cast<std::size_t>(image.cols())
    };
    const std::vector<std::ptrdiff_t> stride_in {
        static_cast<std::ptrdiff_t>(image.cols()
                                    * sizeof(std::complex<double>)),
        sizeof(std::complex<double>)
    };
    const std::vector<std::ptrdiff_t> stride_out {
        static_cast<std::ptrdiff_t>(image.cols() * sizeof(double)),
        sizeof(double)
    };
    const std::vector<std::size_t> axes { 1, 0 };
    pocketfft::c2r(shape,
                   stride_in,
                   stride_out,
                   axes,
                   pocketfft::BACKWARD,
                   image_fft.data(),
                   image.data(),
                   1.0);
}

auto convolve(const Eigen::Ref<const ArrayXXd> image_in,
              const int kernel_n_rows,
              const int kernel_n_cols,
              const Eigen::ArrayXcd& kernel_fft,
              Eigen::ArrayXcd& image_fft,
              ArrayXXd& image_out) -> void
{
    assert(image_in.rows() <= kernel_n_rows
           && "Image dimension must be less than the kernel dimension");
    assert(image_in.cols() <= kernel_n_cols
           && "Image dimension must be less than the kernel dimension");
    // Enlarge image_in to the size of the kernel and pad with zeros
    ArrayXXd image(ArrayXXd::Zero(kernel_n_rows, kernel_n_cols));
    image.topLeftCorner(image_in.rows(), image_in.cols()) = image_in;
    fft_r2c(image, image_fft);
    // Image FFT x kernel FFT
    image_fft *= kernel_fft;
    // Image real
    fft_c2r(image_fft, image);
    const int norm { kernel_n_rows * kernel_n_cols };
    image_out = image.topLeftCorner(image_in.rows(), image_in.cols()) / norm;
}

auto convolveMulti(
  const Eigen::Array<int, Eigen::Dynamic, box::n, Eigen::RowMajor>& edges,
  const ArrayXXd& weights,
  const Eigen::Ref<const ArrayXXd> image_in,
  const Eigen::VectorXi& kernel_rows,
  const Eigen::VectorXi& kernel_cols,
  const std::vector<Eigen::ArrayXcd>& kernels_fft,
  Eigen::ArrayXcd& image_fft,
  ArrayXXd& image_out) -> void
{
    image_out = 0.0;
    for (int i_kernel {}; i_kernel < static_cast<int>(kernels_fft.size());
         ++i_kernel) {
        // Number of rows and column in this subsignal
        const int b { edges(i_kernel, box::b) };
        const int t { edges(i_kernel, box::t) - 1 };
        const int l { edges(i_kernel, box::l) };
        const int r { edges(i_kernel, box::r) - 1 };
        auto weights_sub(weights.row(i_kernel).reshaped<Eigen::RowMajor>(
          image_in.rows(), image_in.cols()));
        const ArrayXXd sub_signal(
          image_in(Eigen::seq(b, t), Eigen::seq(l, r))
          * weights_sub(Eigen::seq(b, t), Eigen::seq(l, r)));
        // Result of convolution using one of the subimages and kernels
        ArrayXXd conv_result(sub_signal.rows(), sub_signal.cols());
        convolve(sub_signal,
                 kernel_rows(i_kernel),
                 kernel_cols(i_kernel),
                 kernels_fft[i_kernel],
                 image_fft,
                 conv_result);
        // The full convolution is a sum over all convolutions
        image_out(Eigen::seq(b, t), Eigen::seq(l, r)) += conv_result;
    }
}

} // namespace tango
