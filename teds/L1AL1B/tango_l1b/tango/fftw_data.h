// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include <fftw3.h>

class FFTWData
{
private:
    const int kernel_size;
    const int kernel_fft_size;
    fftw_plan forward_plan;
    fftw_plan reverse_plan;

public:
    enum Direction
    {
     FORWARD,
     REVERSE,
    };

    double *image;
    fftw_complex *image_fft;

    FFTWData(const int kernel_n_rows, const int kernel_n_cols);
    void execute(const Direction direction) const;
    ~FFTWData();
};
