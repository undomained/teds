// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing the calibration key data

#pragma once

#include "linear_spline.h"

#include <complex>

namespace tango {

class CKD
{
private:
    // Generate integer pixel indices and weights from the fractional
    // row indices. The argument is the width of a spectrum on the
    // detector in units of pixel rows.
    auto genPixelIndices(const double spectrum_width) -> void;

public:
    // Construct the CKD from a NetCDF file. The second argument
    // defines the number of pixels comprising one spectral element.
    CKD(const std::string& filename, const double spectrum_width = 1.0);

    // Number of detector pixels in the spatial direction
    int n_detector_rows {};
    // Number of detector pixels in the spectral direction
    int n_detector_cols {};
    // Number of pixels in one unbinned image. Shorthand for
    // (n_detector_rows x n_detector_cols)
    int npix {};
    // Number of L1B spectra (spatial samples across track)
    int n_act {};
    // Number of samples per L1B (spectral samples)
    int n_lbl {};
    int n_isrf_samples;

    // Bad pixel mask
    std::vector<bool> pixel_mask {};

    // Dark
    struct
    {
        // This and other processing steps are enabled if the relevant
        // section is found in the CKD file.
        bool enabled { false };
        // Dark offset (independent of integration time)
        std::vector<double> offset {};
        // Dark current per second of integration time
        std::vector<double> current {};
    } dark;

    // Noise
    struct
    {
        bool enabled { false };
        // Conversion gain (signal dependent noise term)
        std::vector<double> g {};
        // Square or read noise (signal independent noise term)
        std::vector<double> n2 {};
    } noise;

    // Nonlinearity
    struct
    {
        bool enabled { false };
        // Nonlinearity is stored in a single pixel-independent spline
        LinearSpline spline {};
    } nonlin;

    // PRNU
    struct
    {
        bool enabled { false };
        // Photoresponse non-uniformity
        std::vector<double> prnu {};
    } prnu;

    // Stray light
    struct
    {
        bool enabled { false };
        // Number of stray light kernels
        int n_kernels {};
        // Number of rows of each kernel
        std::vector<int> kernel_rows {};
        // Number of columns of each kernel
        std::vector<int> kernel_cols {};
        // The FFT array size of each kernel
        std::vector<int> kernel_fft_sizes {};
        // Fourier transforms of the kernels
        std::vector<std::vector<std::complex<double>>> kernels_fft {};
        // Total internal scattering factor
        std::vector<double> eta {};
        // Kernel weights
        std::vector<std::vector<double>> weights {};
        // Boundaries of subimages that must be extracted for the
        // convolutions. The order of coefficients is 'bottom', 'top',
        // 'left', 'right'.
        std::vector<int> edges {};
    } stray;

    // Swath
    struct
    {
        bool enabled { false };
        // Row indices at which to extract each spectrum
        // Dimensions: (n_act, n_detector_cols)
        std::vector<std::vector<double>> row_indices {};
        // Number of pixels from which one spectral element is constructed
        int n_indices {};
        // Pixels from which one spectral element is constructed
        std::vector<std::vector<int>> pix_indices {};
        // Weights of pixels comprising a given spectral element
        std::vector<std::vector<double>> weights {};
        // Line of sight vectors
        // Dimensions: (n_act, 3)
        std::vector<double> los {};
    } swath;
    auto genFovIndices(const double spectrum_width) -> void;

    // Spectral
    struct
    {
        bool enabled { false };
        // Wavelengths assigned to each detector column of each L1B spectrum.
        // Dimensions: (n_act, n_detector_cols).
        std::vector<std::vector<double>> wavelength {};
    } wave;

    // Radiometric
    struct
    {
        bool enabled { false };
        // Radiometric calibration constants
        // Dimensions: (n_act, n_detector_cols).
        std::vector<std::vector<double>> rad {};
        std::vector<std::vector<std::vector<double>>> isrf {};
        std::vector<std::vector<std::vector<double>>> isrf_wl {};
    } rad;

    ~CKD() = default;
};

} // namespace tango
