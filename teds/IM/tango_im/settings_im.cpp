// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "settings_im.h"

#include <tango_l1b/io.h>

namespace tango {

auto SettingsIM::scanKeys() -> void
{
    scan(instrument);
    scan(processing_version);
    scan(cal_level);
    scan(image_start);
    scan(image_end);

    scan(detector.binning_table_id);
    scan(detector.exposure_time);
    scan(detector.nr_coadditions);

    scan(isrf.enabled);
    scan(isrf.fwhm_gauss);

    scan(dark.enabled);

    scan(noise.enabled);
    scan(noise.seed);

    scan(prnu.enabled);

    scan(stray.enabled);

    scan(rad.enabled);

    scan(io.ckd);
    scan(io.binning_table);
    scan(io.l1a);
    scan(io.sgm);
}

auto SettingsIM::checkParameters() -> void
{
    // These files must always be present
    checkPresenceOfFile(io.ckd, true);
    checkPresenceOfFile(io.binning_table, true);
    checkPresenceOfFile(io.sgm, true);

    checkFileWritable(io.l1a);
}

} // namespace tango
