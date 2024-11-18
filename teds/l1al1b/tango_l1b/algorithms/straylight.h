// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "base_algo.h"
#include "../fourier.h"

namespace tango {

class CKD;
class L1;

/// Dummy Correction class

class Straylight : public BaseAlgo {
public:

    /// Constructor.
    Straylight() = default;

    /// Destructor.
    ~Straylight() = default;

    /// Return the name of the class.
    std::string getName() const override;

    bool isInteger(double N);

    /// Retrieve the required datasets
    bool algoCheckInput(const CKD& ckd, L1& l1) override;

//    /// Set all loaded data to null.
//    virtusl void unloadData();

    /// Perform the algorithm
    void algoExecute(L1& l1, const Dataset& input_data) override;

    // Bin image to straylight kernel dimensions
    void binImage(const CKD& ckd, const std::vector<double>& image, std::vector<double>& image_binned);

    // Interpolate image to original dimensions
    void unbinImage(const CKD& ckd, const std::vector<double>& image, std::vector<double>& image_unbinned);

    // Convole kernels with image
    void convolveKernels(
        const CKD& ckd, 
        const std::vector<double>& image_binned, 
        std::vector<std::complex<double>>& image_fft, 
        std::vector<double>& conv_result);
    
};
} // namespace tango
