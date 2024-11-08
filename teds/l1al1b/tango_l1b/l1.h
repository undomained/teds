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
    std::vector<double> stdev; // this is used for noise, is not the standard deviation
    // TODO: Should we make wavelenght 1D as well? For improved perfomance and consistency
    std::shared_ptr<std::vector<std::vector<double>>> wavelength; //shared_ptr, same wl for all alt
    std::string units;

    // Detector image meta data
    double image_time;
    uint8_t binning_table_id;
    uint16_t nr_coadditions;
    double exposure_time;
    int i_alt;

    // Carbon Instrument 
    std::vector<Spectrum> spectra; 
    std::vector<int> image_i32;

    // Nitro instrument, 
    // TODO: Should we make these 1D as well? For improved perfomance and consistency
    std::vector<std::vector<double>> observation_sig; // spectra.signal alternative for nitro
    std::vector<std::vector<double>> observation_std; // spectra.stdev alternative for nitro
    std::shared_ptr<std::vector<std::vector<double>>> observation_wl; // wavelength of incoming line-by-line spectra
    std::shared_ptr<std::vector<std::vector<double>>> wavelength_binned; //for wavelength map in case of binning

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

    // Current calibration level
    ProcLevel level;

    // Constructor
    L1();
};


} // namespace tango
