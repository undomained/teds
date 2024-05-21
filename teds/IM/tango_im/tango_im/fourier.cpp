#include "header.h"
#include "functions.h"
#include "fourier.h"

#include <cassert>
#include <complex>
#include <fftw3.h>

namespace tango {

void convolve_fft( // {{{
    const int image_n_rows,
    const int image_n_cols,
    const vector<double> &image_in,
    const int kernel_n_rows,
    const int kernel_n_cols,
    const int kernel_fft_size,
    const vector<complex<double> > &kernel_fft,
    vector<double> &image_out
)
{
    const int kernel_size { kernel_n_rows*kernel_n_cols };
    assert(image_n_rows <= kernel_n_rows && "Image dimension must be less than the kernel dimension");
    assert(image_n_cols <= kernel_n_cols && "Image dimension must be less than the kernel dimension");
    vector<double> image_wrk(kernel_size, 0.0);
    // Expand image to the size of the kernel
    for (int i {}; i < image_n_rows; ++i) {
        for (int j {}; j < image_n_cols; ++j) {
            image_wrk[i*kernel_n_cols+j] = image_in[i*image_n_cols+j];
        }
    }
    for (int i {}; i < kernel_n_rows; ++i) {
        for (int j = image_n_cols; j < kernel_n_cols; ++j) {
            image_wrk[i*kernel_n_cols+j] = 0.0;
        }
    }
    for (int i = image_n_rows; i < kernel_n_rows; ++i) {
        for (int j {}; j < image_n_cols; ++j) {
            image_wrk[i*kernel_n_cols+j] = 0.0;
        }
    }

    vector<complex<double> > image_fft(kernel_fft_size, 0.0);

    fftw_plan fft_plan = fftw_plan_dft_r2c_2d (
        kernel_n_rows, kernel_n_cols,
        image_wrk.data(),
        reinterpret_cast<fftw_complex*>(image_fft.data()),
        FFTW_ESTIMATE

    );

    fftw_execute(fft_plan);
    fftw_destroy_plan(fft_plan);

    for (int i {}; i < kernel_fft_size; ++i) image_fft[i] *= kernel_fft[i];// * by elements wise, not as arrays(matrixes)

    fft_plan = fftw_plan_dft_c2r_2d (
        kernel_n_rows, kernel_n_cols,
        reinterpret_cast<fftw_complex*>(image_fft.data()),
        image_wrk.data(),
        FFTW_ESTIMATE
    );

    fftw_execute(fft_plan);
    fftw_destroy_plan(fft_plan);

    image_out.resize(image_n_rows * image_n_cols);
    for (int i {}; i < image_n_rows; ++i) {
        const int row_idx { i*image_n_cols };
        const int kernel_row_idx { i*kernel_n_cols };
        for (int j {}; j < image_n_cols; ++j) {
            image_out[row_idx+j] = image_wrk[kernel_row_idx+j]/kernel_size;
        }
    }

} // }}}

} // namespace tango
