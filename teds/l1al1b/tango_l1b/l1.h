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

// #include "binning_table.h"

#pragma once
#include "constants.h"
#include "spectrum.h"
#include <memory>



namespace tango {

struct L1 {
    // Science data
    std::vector<bool> pixel_mask;
    std::vector<double> image;
    std::vector<double> noise2; // squared to indicate variance
    std::vector<int> image_i32;
    std::vector<double> stdev; // stdev is a statistical property, using noise2

    std::vector<Spectrum> spectra; // for carbon

    std::vector<std::vector<double>> observation_sig; // spectra.signal alternative for nitro
    std::vector<std::vector<double>> observation_std; // spectra.stdev alternative for nitro

    // shared pointer between all instances of l1 because they are all the same for every along track
    std::shared_ptr<std::vector<std::vector<double>>> lbl_wavelength; // wavelength of incoming line-by-line spectra
    std::shared_ptr<std::vector<std::vector<double>>> wavelength;

    // Geolocation data
    struct {
        std::vector<float> lat;
        std::vector<float> lon;
        std::vector<float> height;
        std::vector<float> vza; // Viewing zenith angle
        std::vector<float> vaa; // Viewing azimuth angle
        std::vector<float> sza; // Solar zenith angle
        std::vector<float> saa; // Solar azimuth angle
    } geo;

    // Detector image meta data
    double image_time;
    uint8_t binning_table_id;
    uint16_t nr_coadditions;
    double exposure_time;

    // Current calibration level
    ProcLevel level;

    // Constructor
    L1();
};


} // namespace tango
