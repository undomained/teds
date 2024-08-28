// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "straylight.h"
#include <cstdint>
#include <spdlog/spdlog.h>
#include <ranges>

namespace tango {

//Straylight::Straylight() {
//}

//Straylight::~Straylight() {
//}

std::string Straylight::getName() const {
    return std::string("Straylight");
}

bool Straylight::algoCheckInput(const CKD& ckd, L1& l1) {
    spdlog::info("Straylight algoCheckInput fct still to be filled in");

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


void Straylight::algoExecute(const CKD& ckd, L1& l1) {
    fillHoles(l1.pixel_mask, l1.image);
    // Unbin the image using the stray light binning table
    std::vector<double> image_unbin(ckd.npix);
    //if (binned_detector_image) {
    //    binning_table.unbin(l1.image, image_unbin);
    //} else {
    //    image_unbin = std::move(l1.image);
    //}
    // "Ideal" image, i.e. the one without stray light
    std::vector<double> image_ideal { image_unbin };
    // Part of the unbinned detector image
    std::vector<double> sub_image {};
    // Result of convoling one kernel with a subimage
    std::vector<double> conv_result_sub {};
    // Allocate the image FFT once to save time
    std::vector<std::complex<double>> image_fft(
      *std::max_element(ckd.stray.kernel_fft_sizes.begin(), ckd.stray.kernel_fft_sizes.end()));
    // Result of taking the whole convolution
    std::vector<double> conv_result(ckd.npix);
    // Van Cittert algorithm
    int n_van_cittert = 2;
    spdlog::warn("Straylight: Move n_van_cittert to ckd");
    spdlog::warn("Straylight: binned images not implemented yet");
    for (int i_vc {}; i_vc < n_van_cittert; ++i_vc) {
        std::fill(conv_result.begin(), conv_result.end(), 0.0);
        for (int i_kernel {}; i_kernel < ckd.stray.n_kernels; ++i_kernel) {
            // Number of rows and column in this subimage
            const int image_n_rows {
                ckd.stray.edges[i_kernel * box::n + box::t]
                - ckd.stray.edges[i_kernel * box::n + box::b]
            };
            const int image_n_cols {
                ckd.stray.edges[i_kernel * box::n + box::r]
                - ckd.stray.edges[i_kernel * box::n + box::l]
            };
            // Each iteration starts with image_ideal as the best
            // current estimate.
            sub_image.resize(image_n_rows * image_n_cols);
            for (int i {}; i < image_n_rows; ++i) {
                for (int j {}; j < image_n_cols; ++j) {
                    const int idx {
                        (i + ckd.stray.edges[i_kernel * box::n + box::b])
                          * ckd.n_detector_cols
                        + j + ckd.stray.edges[i_kernel * box::n + box::l]
                    };
                    sub_image[i * image_n_cols + j] =
                      image_ideal[idx] * ckd.stray.weights[i_kernel][idx];
                }
            }
            // Result of taking a convolution using one of the
            // subimages and kernels
            convolve(image_n_rows,
                     image_n_cols,
                     sub_image,
                     ckd.stray.kernel_rows[i_kernel],
                     ckd.stray.kernel_cols[i_kernel],
                     ckd.stray.kernels_fft[i_kernel],
                     image_fft,
                     conv_result_sub);
            // The full convolution is a sum over all convolutions
            for (int i {}; i < image_n_rows; ++i) {
                for (int j {}; j < image_n_cols; ++j) {
                    conv_result
                      [(i + ckd.stray.edges[i_kernel * box::n + box::b])
                         * ckd.n_detector_cols
                       + j + ckd.stray.edges[i_kernel * box::n + box::l]] +=
                      conv_result_sub[i * image_n_cols + j];
                }
            }
        }
        for (int i {}; i < ckd.npix; ++i) {
            image_ideal[i] =
              (image_unbin[i] - conv_result[i]) / (1 - ckd.stray.eta[i]);
        }
    }
    l1.image = std::move(image_ideal);
    //if (binned_detector_image) {
    //    binning_table.bin(image_ideal, l1.image);
    //} else {
    //    l1.image = std::move(image_ideal);
    //}
    
}



} // namespace tango
