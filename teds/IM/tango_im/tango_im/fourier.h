#pragma once

#include "header.h"
#include <complex>

namespace tango {

void convolve_fft(
    const int image_n_rows,
    const int image_n_cols,
    const vector<double> &image_in,
    const int kernel_n_rows,
    const int kernel_n_cols,
    const int kernel_fft_size,
    const vector<complex<double>> &kernel_fft,
    vector<double> &image_out
);

} // namespace tango
