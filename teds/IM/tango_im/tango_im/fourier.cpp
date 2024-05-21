#include "header.h"
#include "functions.h"
#include "fourier.h"

#include <cassert>
#include <complex>
#include <fftw3.h>

int complex_fft( // {{{
    size_t sz, // Number of complex elements. The size of your array is twice this number. Must be a power of two.
    double *x, // Pointer to first element of array. Turns into the Fourier transform.
    bool inverse, // Flag for performing inverse Fourier transform.
    ptrdiff_t elshift, // Index shift to next complex element.
    ptrdiff_t imagshift // Index shift from real to imaginary part of one number.
)
{

    // Assert that sz is a power of two.
    size_t attempt = 1;
    while (attempt != sz) {
        attempt *= 2;
        if (attempt > sz) return 1; // Number of elements should be a power of two.
    }
    // Step one. Do the re-ordering.
    // Element i should be switched with element i counted with the bits reversed.
    size_t idx_rev = 0; // Reverse index.
    for (size_t idx=0 ; idx<sz ; idx++) { // Forward index.
        // We start with both indices at zero. That is fine, because all bits are zero.
        // Indices idx and idx_rev should be flipped, the pairs will occur twice, one
        // with idx and idx_rev reversed. We only want to flip once.
        if (idx_rev > idx) {
            // Flip real parts.
            double flip = x[idx*elshift];
            x[idx*elshift] = x[idx_rev*elshift];
            x[idx_rev*elshift] = flip;
            // Flip imaginary parts.
            flip = x[idx*elshift+imagshift];
            x[idx*elshift+imagshift] = x[idx_rev*elshift+imagshift];
            x[idx_rev*elshift+imagshift] = flip;
        }
        // Increase the reverse index.
        // Like counting, turn on the bit. If that bit is already turned on,
        // turn it off and continue to the next bit. The only thing is that
        // the bits are reversed.
        size_t bit = sz/2; // Most significant bit.
        while (idx_rev >= bit and bit != 0) {
            // Bit is already on. Turn off and go to next (thus previous) bit.
            idx_rev -= bit;
            bit /= 2;
        }
        // Now turn on this bit.
        idx_rev += bit;
        // If all bits are turned off, we escaped the while loop with bit != 0.
        // That only happens at the end after which idx_rev no longer does anything.
    }

    // Step two. Do the actual work.
    size_t subsize = 1; // Size of E or O.
    size_t comsize = 2; // Size of the combination.
    while (subsize < sz) {
        // The twiddle factor is e^(+/- 2*pi*i*k / comsize)
        // As we iterate over k, we will multiply with e^(+/- 2*pi*i/comsize) every iteration.
        // And we start with one.
        double t_r = 1.0;
        double t_i = 0.0;
        double ang = 360.0*DEGREES/comsize;
        double t_mult_r = cos(ang); // Real part of e^(+/- 2*pi*i/comsize).
        double t_mult_i = sin(ang); // Imaginary part of e^(+ 2*pi*i/comsize), for inverse Fourier transform..
        if (!inverse) t_mult_i *= -1.0; // Turn imagshift part into the one of e^(- 2*pi*i/comsize).
        for (size_t k=0 ; k<subsize ; k++) {
            // This is done for each of the sz/comsize pairs of E and O at the same time.
            // Each of those are a Fourier list that start comsize elements apart.
            // We consider element k now, so we start at element k and iterate by jumping
            // comsize elements so the same k of the next combination.
            for (size_t ek=k ; ek<sz ; ek+=comsize) {
                // Element ek is E_k.
                size_t ok = ek + subsize; // That is element O_k.
                // In the combined list, at element ek, we will have F_k.
                // And at element ok, we will have F_(k+subsize).
                // Thus, at element ek, we have to write ek + t*ok
                // and at element ok, we have to write ek - t*ok.
                double tok_r = t_r * x[ok*elshift] - t_i * x[ok*elshift+imagshift];
                double tok_i = t_i * x[ok*elshift] + t_r * x[ok*elshift+imagshift];
                // As we still need ek, we first overwrite ok.
                x[ok*elshift] = x[ek*elshift] - tok_r;
                x[ok*elshift+imagshift] = x[ek*elshift+imagshift] - tok_i;
                x[ek*elshift] = x[ek*elshift] + tok_r;
                x[ek*elshift+imagshift] = x[ek*elshift+imagshift] + tok_i;
            }
            // Update twiddle factor.
            double temp_r = t_r * t_mult_r - t_i * t_mult_i; // To preserve original t_r one more line.
            t_i = t_i * t_mult_r + t_r * t_mult_i;
            t_r = temp_r;
        }
        // Now, we have twice as large lists.
        subsize *= 2;
        comsize *= 2;
    }
    // Divide by sz for inverse Fourier transform.
    if (inverse) {
        for (size_t k=0 ; k<sz ; k++) {
            x[k*elshift] /= sz;
            x[k*elshift+imagshift] /= sz;
        }
    }

    return 0;

} // }}}

int real_fft( // {{{
    size_t sz, // Size of the array (must be a power of two).
    double *x, // Pointer to first element of array. Turns into the Fourier transform.
    bool inverse, // Flag for performing inverse Fourier transform.
    ptrdiff_t elshift // Index shift to next element.
)
{

    // Nothing to be done for N=1.
    if (sz == 1) return 0;

    // Any other odd sz is illegal, because the underlying complex Fourier transform
    // cannot be executed with a non-integral number of complex numbers.
    if (sz % 2 == 1) return 1; // Number of elements (%zu) should be even.
    // Stricly speaking, sz need not be a power of two. In the implemented complex
    // FFT, it needs to be, but in principle, it is not needed.

    // This is the role of N, because of half-sized complex Fourier transform will be done.
    size_t n = sz/2;
    // Deviation for the two branches. The forward simulation must first do the complex
    // Fourier transform and then the correction calculus. For the inverse, first
    // correction calculus is inverted, followed by an inverse complex Fourier transform.
    int cosfac;
    if (inverse) {
        cosfac = -1; // Turns twiddle factor T into -T*, by negating cosine part.
    } else {
        cosfac = 1; // Preserves twiddle factor T as it is.
        // Complex Fourier transform is done first.
        handle_nonlethal(complex_fft(n,x,inverse,2*elshift,elshift));
    }

    // Set up the twiddle factor gain.
    double ang = 180.0*DEGREES/n; // Angle for twiddle factor.
    double t_mult_r = cos(ang);
    double t_mult_i = -sin(ang);

    // Initialize twiddle factor for k=1.
    double t_r = t_mult_r;
    double t_i = t_mult_i;

    // Apply the correction calculus. Let Z be the complex Fourier transform of the
    // array and F is the desired Fourier transform, then,
    // F_k = (Z_k + Z*_n-k)/2 -iT * (Z_k - Z*_n-k)/2
    // F_n-k = (Z*_k + Z_n-k)/2 -iT* * (Z*_k - Z_n-k)/2
    // Or for the inverse procedure, transform T into -T*, which is done by the
    // cosine factor.
    // For k=N/2, we just have to conjugate, not matter what. Normally, there is
    // a k=N/2, because N is a power of two. For sz=2 not, or when an alternative FFT
    // for complex numbers is implemented, this may be not the case.
    for (size_t k=1 ; 2*k<n ; k++) {

        size_t n_k = n-k;

        // Combinations with plus or minus signs for real and imaginary parts.
        double plus_r = 0.5*(x[2*k*elshift] + x[2*n_k*elshift]);
        double minus_r = plus_r - x[2*n_k*elshift];
        double plus_i = 0.5*(x[2*k*elshift+elshift] + x[2*n_k*elshift+elshift]);
        double minus_i = plus_i - x[2*n_k*elshift+elshift];

        // Write down the matrix and multiply each t_r with the cosine factor, which is
        // one for a normal transform, but minus one for an inverse transform.
        x[2*k*elshift] = plus_r + cosfac * t_r * plus_i + t_i * minus_r;
        x[2*k*elshift+elshift] = minus_i - cosfac * t_r * minus_r + t_i * plus_i;
        x[2*n_k*elshift] = plus_r - cosfac * t_r * plus_i - t_i * minus_r;
        x[2*n_k*elshift+elshift] = -minus_i - cosfac * t_r * minus_r + t_i * plus_i;

        // Update twiddle factor for a one higher k.
        double temp_r = t_r * t_mult_r - t_i * t_mult_i; // To preserve original t_r one more line.
        t_i = t_i * t_mult_r + t_r * t_mult_i;
        t_r = temp_r;
    }

    // The k = N/2 element becomes her complex conjugate.
    if (n%2 == 0) x[n*elshift+elshift] *= -1; // Index is 2*k*elshift+elshift.

    // At last, the last to F_0 and F_N factors. Note the factor 0.5 for the inverse transform.

    if (inverse) {
        double temp_x0 = 0.5 * (x[0] + x[elshift]); // To preserve original x[0] one more line.
        x[elshift] = 0.5 * (x[0] - x[elshift]);
        x[0] = temp_x0;
        // And not forget the complex inverse Fourier transform.
        handle_nonlethal(complex_fft(n,x,inverse,2*elshift,elshift));
    } else {
        double temp_x0 = x[0] + x[elshift]; // To preserve original x[0] one more line.
        x[elshift] = x[0] - x[elshift];
        x[0] = temp_x0;
    }

    return 0;

} // }}}

int real_2d_fft( // {{{
    size_t sz_slow, // Size of the slow dimension (must be a power of two).
    size_t sz_quick, // Size of the quick dimension (must be a power of two).
    double *x, // The data. Will be turned into a representation of the Fourier transform.
    bool inverse, // Flag for performing inverse Fourier transform.
    size_t nslow_nonzero // Number of slow-dimension indices that are not zero. Use zero to revert to full (sz_slow).
)
{
    // Assert that it is not one-dimensional. I do not know if we want to support
    // this in the future. Even a non-power-of-two FFT is more probable.
    if (sz_slow == 1 || sz_quick == 1) // One of the dimensions is only one sized. That is not two-dimensional. Maybe revert to a 1D Fourier transform.

    // Default nslow_nonzero is sz_slow.
    if (nslow_nonzero == 0) nslow_nonzero = sz_slow;

    if (inverse) {
        // In slow dimnension: Left-most two columns are real, the rest is complex.
        double *arr = x;
        for (size_t iquick=0 ; iquick<2 ; iquick++) {
            handle_nonlethal(real_fft(sz_slow,arr,inverse,sz_quick));
            arr++;
        }
        for (size_t iquick=2 ; iquick<sz_quick ; iquick+=2) {
            handle_nonlethal(complex_fft(sz_slow,arr,inverse,sz_quick,1));
            arr += 2;
        }
        // In quick dimension, all are real.
        arr = x;
        for (size_t islow=0 ; islow<nslow_nonzero ; islow++) {
            handle_nonlethal(real_fft(sz_quick,arr,inverse,1));
            arr += sz_quick; // Progress pointer.
        }
    } else {
        // In quick dimension, all are real.
        double *arr = x;
        for (size_t islow=0 ; islow<nslow_nonzero ; islow++) {
            handle_nonlethal(real_fft(sz_quick,arr,inverse,1));
            arr += sz_quick; // Progress pointer.
        }
        // Keep the all-zero arrays zero.
        // In slow dimnension: Left-most two columns are real, the rest is complex.
        arr = x;
        for (size_t iquick=0 ; iquick<2 ; iquick++) {
            handle_nonlethal(real_fft(sz_slow,arr,inverse,sz_quick));
            arr++;
        }
        for (size_t iquick=2 ; iquick<sz_quick ; iquick+=2) {
            handle_nonlethal(complex_fft(sz_slow,arr,inverse,sz_quick,1));
            arr += 2;
        }
    }

    return 0;

} // }}}

int multiply_2d( // {{{
    size_t sz_slow, // Size of the slow dimension (must be a power of two).
    size_t sz_quick, // Size of the quick dimension (must be a power of two).
    double *fact, // Fourier-transformed data to multiply with.
    double *res // Result. Starts with one Fourier transformed data set. Ends like the product.
)
{

    // Real multiplication for the upper left corner.
    for (size_t islow=0 ; islow<2 ; islow++) {
        for (size_t iquick=0 ; iquick<2 ; iquick++) {
            res[islow*sz_quick+iquick] *= fact[islow*sz_quick+iquick];
        }
    }
    // Perform complex multiplications for the complex numbers oriented in the
    // slow-dimension direction. That is only for the first two indices in the
    // quick dimension.
    for (size_t islow=2 ; islow<sz_slow ; islow+=2) {
        for (size_t iquick=0 ; iquick<2 ; iquick++) {
            double temp_r = fact[islow*sz_quick+iquick] * res[islow*sz_quick+iquick] - fact[(islow+1)*sz_quick+iquick] * res[(islow+1)*sz_quick+iquick]; // To have the original real part survive one more line.
            res[(islow+1)*sz_quick+iquick] = fact[islow*sz_quick+iquick] * res[(islow+1)*sz_quick+iquick] + fact[(islow+1)*sz_quick+iquick] * res[islow*sz_quick+iquick];
            res[islow*sz_quick+iquick] = temp_r;
        }
    }
    // Perform complex multiplication for complex numbers oriented in the quick
    // dimension for the whole matrix where the quick index is at least two.
    for (size_t islow=0 ; islow<sz_slow ; islow++) {
        for (size_t iquick=2 ; iquick<sz_quick ; iquick+=2) {
            double temp_r = fact[islow*sz_quick+iquick] * res[islow*sz_quick+iquick] - fact[islow*sz_quick+iquick+1] * res[islow*sz_quick+iquick+1]; // To have the original real part survive one more line.
            res[islow*sz_quick+iquick+1] = fact[islow*sz_quick+iquick] * res[islow*sz_quick+iquick+1] + fact[islow*sz_quick+iquick+1] * res[islow*sz_quick+iquick];
            res[islow*sz_quick+iquick] = temp_r;
        }
    }

    return 0;

} // }}}

int divide_2d( // {{{
    size_t sz_slow, // Size of the slow dimension (must be a power of two).
    size_t sz_quick, // Size of the quick dimension (must be a power of two).
    double *denom, // Fourier-transformed data to divide by.
    double *res // Result. Starts with one Fourier transformed data set. Ends like the product.
)
{

    // Real multiplication for the upper left corner.
    for (size_t islow=0 ; islow<2 ; islow++) {
        for (size_t iquick=0 ; iquick<2 ; iquick++) {
            res[islow*sz_quick+iquick] /= denom[islow*sz_quick+iquick];
        }
    }
    // Perform complex multiplications for the complex numbers oriented in the
    // slow-dimension direction. That is only for the first two indices in the
    // quick dimension.
    for (size_t islow=2 ; islow<sz_slow ; islow+=2) {
        for (size_t iquick=0 ; iquick<2 ; iquick++) {
            double mag2 = pow(denom[islow*sz_quick+iquick],2.0) + pow(denom[(islow+1)*sz_quick+iquick],2.0);
            double temp_r = (denom[islow*sz_quick+iquick] * res[islow*sz_quick+iquick] + denom[(islow+1)*sz_quick+iquick] * res[(islow+1)*sz_quick+iquick]) / mag2; // To have the original real part survive one more line.
            res[(islow+1)*sz_quick+iquick] = (denom[islow*sz_quick+iquick] * res[(islow+1)*sz_quick+iquick] - denom[(islow+1)*sz_quick+iquick] * res[islow*sz_quick+iquick]) / mag2;
            res[islow*sz_quick+iquick] = temp_r;
        }
    }
    // Perform complex multiplication for complex numbers oriented in the quick
    // dimension for the whole matrix where the quick index is at least two.
    for (size_t islow=0 ; islow<sz_slow ; islow++) {
        for (size_t iquick=2 ; iquick<sz_quick ; iquick+=2) {
            double mag2 = pow(denom[islow*sz_quick+iquick],2.0) + pow(denom[islow*sz_quick+iquick+1],2.0);
            double temp_r = (denom[islow*sz_quick+iquick] * res[islow*sz_quick+iquick] + denom[islow*sz_quick+iquick+1] * res[islow*sz_quick+iquick+1]) / mag2; // To have the original real part survive one more line.
            res[islow*sz_quick+iquick+1] = (denom[islow*sz_quick+iquick] * res[islow*sz_quick+iquick+1] - denom[islow*sz_quick+iquick+1] * res[islow*sz_quick+iquick]) / mag2;
            res[islow*sz_quick+iquick] = temp_r;
        }
    }

    return 0;

} // }}}

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

