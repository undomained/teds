// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Struct to hold partially or fully calibrated level 1 data. The data
// level can range from L1A to L1B depending on the calibration

#include "l1.h"

namespace tango {

L1::L1():
    pixel_mask(),
    image(),
    image_i32(),
    stdev(),
    geo(),

    image_time(0.0),
    binning_table_id(0),
    nr_coadditions(0),
    exposure_time(0.0)
    {
}

} // namespace tango
