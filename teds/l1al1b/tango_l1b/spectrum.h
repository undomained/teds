// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for extracting one spectrum from a detector image and
// radiometrically calibrating it.

#pragma once

#include <vector>

namespace tango {

class CKD;

class Spectrum
{
public:
    std::vector<double> signal {};
    std::vector<double> stdev {};
    std::vector<bool> mask {};
    Spectrum() = default;
    // Extract spectra and the corresponding noise values
    auto extract(const CKD& ckd,
                 const std::vector<double>& image,
                 const std::vector<double>& image_stdev,
                 const std::vector<bool>& pixel_mask,
                 const int i_act) -> void;
    // Radiometrically calibrate
    auto calibrate(const CKD& ckd, const double exposure_time, const int i_act)
      -> void;
};

} // namespace tango
