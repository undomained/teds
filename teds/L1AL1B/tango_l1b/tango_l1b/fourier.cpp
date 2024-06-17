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

auto fft_r2c(const int n_rows,
             const int n_cols,
             std::vector<double>& image,
             std::vector<std::complex<double>>& image_fft) -> void
{
    const std::vector<std::size_t> shape { static_cast<std::size_t>(n_rows),
                                           static_cast<std::size_t>(n_cols) };
    const std::vector<std::ptrdiff_t> stride_in {
        static_cast<std::ptrdiff_t>(n_cols * sizeof(double)), sizeof(double)
    };
    const std::vector<std::ptrdiff_t> stride_out {
        static_cast<std::ptrdiff_t>(n_cols * sizeof(std::complex<double>)),
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

auto fft_c2r(const int n_rows,
             const int n_cols,
             std::vector<std::complex<double>>& image_fft,
             std::vector<double>& image) -> void
{
    const std::vector<std::size_t> shape { static_cast<std::size_t>(n_rows),
                                           static_cast<std::size_t>(n_cols) };
    const std::vector<std::ptrdiff_t> stride_in {
        static_cast<std::ptrdiff_t>(n_cols * sizeof(std::complex<double>)),
        sizeof(std::complex<double>)
    };
    const std::vector<std::ptrdiff_t> stride_out {
        static_cast<std::ptrdiff_t>(n_cols * sizeof(double)), sizeof(double)
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

auto convolve(const int image_n_rows,
              const int image_n_cols,
              const std::vector<double>& image_in,
              const int kernel_n_rows,
              const int kernel_n_cols,
              const std::vector<std::complex<double>>& kernel_fft,
              std::vector<std::complex<double>>& image_fft,
              std::vector<double>& image_out) -> void
{
    const int kernel_size { kernel_n_rows * kernel_n_cols };
    assert(image_n_rows <= kernel_n_rows
           && "Image dimension must be less than the kernel dimension");
    assert(image_n_cols <= kernel_n_cols
           && "Image dimension must be less than the kernel dimension");
    // Enlarge image_in to the size of the kernel and pad with zeros
    std::vector<double> image(kernel_size, 0.0);
    for (int i {}; i < image_n_rows; ++i) {
        for (int j {}; j < image_n_cols; ++j) {
            image[i * kernel_n_cols + j] = image_in[i * image_n_cols + j];
        }
    }
    fft_r2c(kernel_n_rows, kernel_n_cols, image, image_fft);
    // Image FFT x kernel FFT
    for (int i {}; i < static_cast<int>(kernel_fft.size()); ++i) {
        image_fft[i] *= kernel_fft[i];
    }
    // Image real
    fft_c2r(kernel_n_rows, kernel_n_cols, image_fft, image);
    // Shrink the resulting image to the original size
    image_out.resize(image_n_rows * image_n_cols);
    const double kernel_size_inv { 1.0 / kernel_size };
    for (int i {}; i < image_n_rows; ++i) {
        for (int j {}; j < image_n_cols; ++j) {
            image_out[i * image_n_cols + j] =
              image[i * kernel_n_cols + j] * kernel_size_inv;
        }
    }
}

} // namespace tango
