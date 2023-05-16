// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "fftw_data.h"

FFTWData::FFTWData(const int kernel_n_rows, const int kernel_n_cols) :
    kernel_size { kernel_n_rows*kernel_n_cols },
    kernel_fft_size { kernel_n_rows*(kernel_n_cols/2+1) },
    image { fftw_alloc_real(kernel_size) },
    image_fft { fftw_alloc_complex(kernel_fft_size) }
{
    const auto &flags { FFTW_ESTIMATE };
    forward_plan = fftw_plan_dft_r2c_2d(kernel_n_rows, kernel_n_cols,
                                        image,
                                        image_fft,
                                        flags);
    reverse_plan = fftw_plan_dft_c2r_2d(kernel_n_rows, kernel_n_cols,
                                        image_fft,
                                        image,
                                        flags);
}

void FFTWData::execute(const Direction direction) const
{
    if (direction == FORWARD) {
        fftw_execute(forward_plan);
    } else {
        fftw_execute(reverse_plan);
    }
}

FFTWData::~FFTWData()
{
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(reverse_plan);
    fftw_free(image);
    fftw_free(image_fft);
}
