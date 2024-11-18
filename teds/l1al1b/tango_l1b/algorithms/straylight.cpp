// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "../datasets.h"
#include "straylight.h"
#include "../cubic_spline.h"

namespace tango {

//Straylight::Straylight() {
//}

//Straylight::~Straylight() {
//}

std::string Straylight::getName() const {
    return std::string("Straylight");
}

bool Straylight::isInteger(double N)
{
    int X = N;
    double rem = N - X;
    if (rem > 0){
        return false;
    }
    return true;
}

bool Straylight::algoCheckInput(const CKD& ckd, L1& l1) {

    int crows = ckd.stray.coarse_rows;
    int ccols = ckd.stray.coarse_cols;
    // Check if binsizes are integers
    double bin_row = static_cast<double>(ckd.n_detector_rows) / crows;
    double bin_col = static_cast<double>(ckd.n_detector_cols) / ccols;
    if ((isInteger(bin_row)) && (isInteger(bin_col)))
    {
        return true;
    }
    spdlog::info("CKD Straylight bins not of integer size");
    return true;
}

//void Straylight::unloadData() {
//    spdlog::info("Straylight unload fct still to be filled in");
//}

static auto fillHoles(const std::vector<bool>& pixel_mask,
                      std::vector<double>& image) -> void
{
    const int i_last { static_cast<int>(image.size() - 1) };
    for (int i {}; i < static_cast<int>(image.size()); ++i) {
        if (pixel_mask[i]) {
            // Unless we are at either end of the image array or one
            // or both of the neighboring pixels are bad, take the
            // average of the neighboring values.
            int n_good {};
            image[i] = 0.0;
            if (i > 0 && !pixel_mask[i - 1]) {
                image[i] += image[i - 1];
                ++n_good;
            }
            if (i < i_last && !pixel_mask[i + 1]) {
                image[i] += image[i + 1];
                ++n_good;
            }
            image[i] /= std::max(1, n_good);
        }
    }
}


void Straylight::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");

    fillHoles(l1.pixel_mask, l1.image);

    // First bin image to coarse resolution of straylight ckd
    int cpix = ckd.stray.coarse_rows * ckd.stray.coarse_cols;
    std::vector<double> image_binned(cpix);
    binImage(ckd, l1.image, image_binned);

    // Allocate the image FFT once to save time
    std::vector<std::complex<double>> image_fft(
      *std::max_element(ckd.stray.kernel_fft_sizes.begin(), ckd.stray.kernel_fft_sizes.end()));
    
    
    std::vector<double> conv_result(cpix); // Result of whole convolution
    std::vector<double> image_result(cpix, 0.0);
    if (getModelType() == "L1B")
    { // L1 implementation (van Cittert Convolution Subtraction)
        
        std::vector<double> image_ideal {image_binned}; // Start with binned image

        // Van Cittert algorithm
        for (int i_vc = 0; i_vc < ckd.stray.n_van_cittert; i_vc++){            
            
            // Clear previous convolution results
            std::fill(conv_result.begin(), conv_result.end(), 0.0);

            // Perform kernel convolutions and update conv_result
            convolveKernels(ckd, image_ideal, image_fft, conv_result);

            // Update the estimate of the "ideal" image
            for (int i = 0; i < cpix; ++i) 
            {
                image_ideal[i] = (image_binned[i] - conv_result[i]) / (1 - ckd.stray.eta[i]);
            }
        }
        image_result = std::move(image_ideal);

    } else if (getModelType() == "IM")
    {   // IM implementation (single Convolution Addition)
        convolveKernels(ckd, image_binned, image_fft, conv_result);
        for (int i = 0; i < cpix; ++i) 
        {
            image_result[i] = (1 - ckd.stray.eta[i]) * image_binned[i] + conv_result[i];
        }
    }
    
    std::vector<double> image_ideal_unbinned(ckd.npix);
    unbinImage(ckd, image_result, image_ideal_unbinned);

    l1.image = std::move(image_ideal_unbinned);
}

void Straylight::convolveKernels(
    const CKD& ckd, 
    const std::vector<double>& image_binned, 
    std::vector<std::complex<double>>& image_fft,
    std::vector<double>& conv_result)
{
    int binned_rows = ckd.stray.coarse_rows;
    int binned_cols = ckd.stray.coarse_cols;
    int cpix = binned_rows * binned_cols;


    for (int i_kernel = 0; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
        int edge_b = ckd.stray.edges[i_kernel * 4];
        int edge_t = ckd.stray.edges[i_kernel * 4 + 1];
        int edge_l = ckd.stray.edges[i_kernel * 4 + 2];
        int edge_r = ckd.stray.edges[i_kernel * 4 + 3];
        int image_n_rows = edge_t - edge_b;
        int image_n_cols = edge_r - edge_l;

        std::vector<double> sub_image(image_n_rows * image_n_cols);
        for (int i = 0; i < image_n_rows; ++i) {
            for (int j = 0; j < image_n_cols; ++j) {
                int idx = (i + edge_b) * binned_cols + j + edge_l;
                sub_image[i * image_n_cols + j] = image_binned[idx] * ckd.stray.weights[i_kernel][idx];
            }
        }

        std::vector<double> conv_result_sub;
        convolve(image_n_rows, image_n_cols, sub_image,
                    ckd.stray.kernel_rows[i_kernel],
                    ckd.stray.kernel_cols[i_kernel],
                    ckd.stray.kernels_fft[i_kernel],
                    image_fft,
                    conv_result_sub);

        // Total result is sum of contributions of all kernels
        for (int i = 0; i < image_n_rows; ++i) 
        {
            int row_ix = i + edge_b;
            for (int j = 0; j < image_n_cols; ++j) 
            {
                int col_ix = j + edge_l;
                conv_result[row_ix * binned_cols + col_ix] += conv_result_sub[i * image_n_cols + j];
            }
        }
    }
}


void Straylight::binImage(const CKD& ckd, const std::vector<double>& image, std::vector<double>& image_binned) 
{
    int binned_rows = ckd.stray.coarse_rows;
    int binned_cols = ckd.stray.coarse_cols;
    int binsize_row = static_cast<double>(ckd.n_detector_rows) / binned_rows;
    int binsize_col = static_cast<double>(ckd.n_detector_cols) / binned_cols;

    for (int br = 0; br < binned_rows; br++) {
        for (int bc = 0; bc < binned_cols; bc++) {
            double bin_sum = 0.0;
            for (int r = 0; r < binsize_row; r++) {
                for (int c = 0; c < binsize_col; c++) {
                    int row_index = br * binsize_row + r;
                    int col_index = bc * binsize_col + c;
                    int index = row_index * ckd.n_detector_cols + col_index;
                    bin_sum += image[index];
                }
            }
            image_binned[br * binned_cols + bc] = bin_sum;
        }
    }
}


void Straylight::unbinImage(
    const CKD& ckd, const std::vector<double>& image, 
    std::vector<double>& image_unbinned) {
    
    int binned_rows = ckd.stray.coarse_rows;
    int binned_cols = ckd.stray.coarse_cols;
    int binsize_row = static_cast<double>(ckd.n_detector_rows) / binned_rows;
    int binsize_col = static_cast<double>(ckd.n_detector_cols) / binned_cols;

    std::vector<double> binned_col_indices(binned_cols);
    for (int i_col = 0; i_col < binned_cols; i_col++) {
        binned_col_indices[i_col] = i_col * binsize_col + binsize_col / 2;
    }

    std::vector<double> binned_row_indices(binned_rows);
    std::vector<double> image_buf(binned_rows * ckd.n_detector_cols);

    for (int i_row = 0; i_row < binned_rows; i_row++) {
        binned_row_indices[i_row] = i_row * binsize_row + binsize_row / 2;
        std::vector<double> thisrow(binned_cols);
        for (int i_col = 0; i_col < binned_cols; i_col++) {
            thisrow[i_col] = image[i_row * binned_cols + i_col];
        }
        const CubicSpline spline {binned_col_indices, thisrow};
        for (int i_col = 0; i_col < ckd.n_detector_cols; i_col++) {
            image_buf[i_row * ckd.n_detector_cols + i_col] = spline.eval(i_col);
        }
    }

    for (int i_col = 0; i_col < ckd.n_detector_cols; i_col++) {
        std::vector<double> thiscol(binned_rows);
        for (int i_row = 0; i_row < binned_rows; i_row++) {
            thiscol[i_row] = image_buf[i_row * ckd.n_detector_cols + i_col];
        }
        const CubicSpline spline {binned_row_indices, thiscol};
        for (int i_row = 0; i_row < ckd.n_detector_rows; i_row++) {
            image_unbinned[i_row * ckd.n_detector_cols + i_col] = spline.eval(i_row);
        }
    }

    // Normalize
    double bin_area = binsize_row * binsize_col;
    for (int i = 0; i < ckd.npix; ++i) {
        image_unbinned[i] /= bin_area;
    }
}


} // namespace tango
