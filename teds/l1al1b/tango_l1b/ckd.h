// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing the calibration key data

#pragma once

#include "linear_spline.h"

#include <complex>

namespace tango {

class BinningTable;

class CKD
{
public:
    CKD() = default;
    // Construct the CKD from a NetCDF file. The second argument
    // defines the number of pixels comprising one spectral element.
    CKD(const std::string& filename);
    // Bin all CKDs and update the binned detector dimensions
    auto bin(const BinningTable& binning_table) -> void;

    // Number of detector pixels in the spatial direction
    int n_detector_rows {};
    // Number of detector pixels in the spectral direction
    int n_detector_cols {};
    // Number of pixels in one unbinned image. Shorthand for
    // (n_detector_rows x n_detector_cols)
    int npix {};
    // Binned detector dimensions
    int n_detector_rows_binned {};
    int n_detector_cols_binned {};
    int npix_binned {};
    // Number of L1B spectra (spatial samples across track)
    int n_act {};
    // Number of L1B wavelengths
    int n_wavelengths {};

    // Bad pixel mask
    std::vector<bool> pixel_mask {};

    // Dark
    struct
    {
        // Dark offset (independent of integration time)
        std::vector<double> offset {};
        // Dark current per second of integration time
        std::vector<double> current {};
    } dark;

    // Noise
    struct
    {
        // Conversion gain (signal dependent noise term)
        std::vector<double> g {};
        // Square of read noise (signal independent noise term)
        std::vector<double> n2 {};
    } noise;

    // Nonlinearity
    struct
    {
        // Nonlinearity is stored in a single pixel-independent spline
        LinearSpline spline {};
    } nonlin;

    // PRNU
    struct
    {
        // Photoresponse non-uniformity
        std::vector<double> prnu {};
    } prnu;

    // Stray light
    struct
    {
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
        // Across track angles
        std::vector<double> act_angles {};
        // ACT angle of each detector pixel
        std::vector<double> act_map {};
        // Wavelength of each detector pixel
        std::vector<double> wavelength_map {};
        // Row index of each L1B spectral element
        std::vector<double> row_map {};
        // Column index of each L1B spectral element
        std::vector<double> col_map {};
        // Line of sight vectors
        // Dimensions: (n_act, 3)
        std::vector<double> los {};
    } swath;

    // Spectral
    struct
    {
        // Wavelengths assigned to each detector column of each L1B spectrum.
        // Dimensions: (n_act, n_wavelengths).
        std::vector<double> wavelengths {};
    } wave;

    // Radiometric
    struct
    {
        // Radiometric calibration constants
        // Dimensions: (n_act, n_detector_cols).
        std::vector<std::vector<double>> rad {};
    } rad;

    ~CKD() = default;
};

} // namespace tango
