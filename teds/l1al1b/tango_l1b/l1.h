// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Struct to hold partially or fully calibrated level 1 data. The data
// level can range from L1A to L1B depending on the calibration
// level. The relevant data is held in these member variables:
//   image_i32 - detector image (int32) for level L1A
//   image - detector image (float64) for levels [raw, stray] (inclusive)
//   stdev - detector noise for levels [raw, stray]
//   spectra - detector spectra for levels [swath, l1b]
//   wavelength - corresponding wavelengths [swath, l1b]

#pragma once

#include "constants.h"

#include <memory>

namespace tango {

struct L1
{
    // Science data
    std::vector<double> image {};
    std::vector<int> image_i32 {};
    std::vector<double> stdev {};
    std::vector<double> spectra {};
    std::vector<double> spectra_stdev {};
    std::vector<double> spectra_mask {};
    // Wavelength grid currently associated with the spectra (either
    // line-by-line or that of the CKD). It is the same for all
    // images.
    std::shared_ptr<std::vector<std::vector<double>>> wavelengths {};

    // Geolocation data
    struct
    {
        std::vector<float> lat {};
        std::vector<float> lon {};
        std::vector<float> height {};
        std::vector<float> vza {}; // Viewing zenith angle
        std::vector<float> vaa {}; // Viewing azimuth angle
        std::vector<float> sza {}; // Solar zenith angle
        std::vector<float> saa {}; // Solar azimuth angle
    } geo;

    // Detector image meta data
    double image_time {};
    uint8_t binning_table_id {};
    uint16_t nr_coadditions {};
    double exposure_time {};

    // Current calibration level
    ProcLevel level { ProcLevel::l1b };
};

} // namespace tango
