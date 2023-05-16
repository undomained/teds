// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef FOURIER_H
#define FOURIER_H

#include "header.h"
#include <complex>

// Fourier transform of an array of complex numbers. The numbers can be
// inside a bigger array.
int complex_fft(
    size_t sz, // Size of the array (must be a power of two).
    double *x, // Pointer to first element of array. Turns into the Fourier transform.
    bool inverse=false, // Flag for performing inverse Fourier transform.
    ptrdiff_t elshift=2, // Index shift to next complex element.
    ptrdiff_t imagshift=1 // Index shift from real to imaginary part of one number.
);

// Fourier transform of real array. First two elements are real parts of
// F_0 and F_N/2. Followed by real and imaginary parts of F_1, F_2, until
// F_N/2-1. This uses the symmetry of mirror-complex conjugate about F_N/2
// (or about F_N=F_0, note periodic boundary conditions).
int real_fft(
    size_t sz, // Size of the array (must be a power of two).
    double *x, // Pointer to first element of array. Turns into the Fourier transform.
    bool inverse=false, // Flag for performing inverse Fourier transform.
    ptrdiff_t elshift=1 // Index shift to next element.
);

// Perform real-number two-dimensional Fourier transform for real numbers. The
// corresponding Fourier transform will be.
// R 0 0        R 0 X/2        R 0 1     I 0 1     R 0 2    ...    I 0 X/2-1
// R Y/2 0      R Y/2 X/2      R 1 1     I 1 1     R 1 2    ...    I 1 X/2-1
// R 1 0        R 1 X/2        R 2 1     I 2 1     R 2 2    ...    I 2 X/2-1
// I 1 0        I 1 X/2        R 3 1     I 3 1     R 3 2    ...    I 3 X/2-1
// R 2 0        R 2 X/2        R 4 1     I 4 1     R 4 2    ...    I 4 X/2-1
// ............................................................................
// I Y/2-1 0    I Y/2-1 X/2    R Y-1 1   I Y-1 1   R Y-1 2  ...    I Y-1 X/2-1
// Here, the horizontal direction is the quick dimension. This means that most
// of the complex pairs are adjacent, but the first two columns have vertically
// adjacent complex pairs with pointer difference sz_quick. The upper left corner
// only have real parts. If you are really interested in the Fourier transform,
// F x y = F* X-x Y-y. We only use this matrix representation to know how to
// multiply. Subsequently, we will inverse Fourier transform as quickly as possible.
int real_2d_fft(
    size_t sz_slow, // Size of the slow dimension (must be a power of two).
    size_t sz_quick, // Size of the quick dimension (must be a power of two).
    double *x, // The data. Will be turned into a representation of the Fourier transform.
    bool inverse=false, // Flag for performing inverse Fourier transform.
    size_t nslow_nonzero=0 // Number of slow-dimension indices that are not zero. Use zero to revert to full (sz_slow).
);

// Perform a proper multiplication of two equally sized Fourier transforms.
// It assumes the 2D Fourier transforms were created with real_2d_fft or has
// the same conventions on the orientations of the complex numbers.
int multiply_2d(
    size_t sz_slow, // Size of the slow dimension (must be a power of two).
    size_t sz_quick, // Size of the quick dimension (must be a power of two).
    double *fact, // Fourier-transformed data to multiply with.
    double *res // Result. Starts with one Fourier transformed data set. Ends like the product.
);

// Perform a proper division of two equally sized Fourier transforms.
// It assumes the 2D Fourier transforms were created with real_2d_fft or has
// the same conventions on the orientations of the complex numbers.
int divide_2d(
    size_t sz_slow, // Size of the slow dimension (must be a power of two).
    size_t sz_quick, // Size of the quick dimension (must be a power of two).
    double *denom, // Fourier-transformed data to divide by.
    double *res // Result. Starts with one Fourier transformed data set. Ends like the product.
);

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

#endif
