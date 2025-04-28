// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Struct to hold partially or fully calibrated level 1 data. The data
// level can range from L1A to L1B depending on the calibration
// level. All detector image or spectra, depending on the data level,
// are held in a single L1 instance.

#pragma once

#include "eigen.h"
#include "geometry.h"

#include <Eigen/Geometry>
#include <array>
#include <memory>

namespace tango {

struct L1
{
    // Number of along-track positions. This is the cropped number
    // after alt_beg and alt_end have been applied. Number of
    // across-track positions is stored in the CKD and is not a user
    // parameter.
    int n_alt {};

    // Science data
    // Detector signal for all along track (ALT) positions. Units
    // depend on the calibration level. Dimensions are ALT x bins
    // where the number of bins per detector image depends on
    // binning. When using Eigen slicing, signal.row(i) returns the
    // ith detector image, not a detector row.
    ArrayXXd signal {};
    // Detector signal noise (stdev) levels
    ArrayXXd noise {};
    // Partially or fully calibrated spectra extracted from detector
    // images. Dimensions are (ACT x ALT) x wavelengths. When slicing,
    // spectra.row(i) returns all the spectrum corresponding to one
    // ALT and ACT positions.
    ArrayXXd spectra {};
    // Spectra noise values
    ArrayXXd spectra_noise {};
    // Wavelength grid associated with spectra. Can change between
    // calibration steps. For instance, the SGM data has the
    // line-by-line grid but otherwise the grid is read from the CKD.
    Eigen::ArrayXd wavelengths {};

    // Navigation data
    // Orbit positions
    std::vector<double> orb_pos0 {};
    ArrayXNd<dims::vec> orb_pos {};
    // Attitude quaternions
    std::vector<Eigen::Quaterniond> att_quat {};

    // Geolocation data
    Geometry geo {};

    // Detector image meta data
    Eigen::ArrayXd time {};
    std::vector<uint32_t> tai_seconds {};
    std::vector<double> tai_subsec {};
    uint8_t binning_table_id { 1 };
    uint16_t nr_coadditions { 1 };
    double exposure_time {};

    // Current calibration level
    ProcLevel level { ProcLevel::l1b };
};

} // namespace tango
