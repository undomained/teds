// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "settings_l1b.h"

#include "io.h"

namespace tango {

auto SettingsL1B::scanKeys() -> void
{
    scan(instrument);
    scan(processing_version);
    scan(cal_level);
    scan(reverse_wavelength);
    scan(image_start);
    scan(image_end);
    scan(unbinning);
    scan(proctable.file);
    scan(proctable.algo_list);

    scan(dark.enabled);

    scan(noise.enabled);

    scan(nonlin.enabled);

    scan(prnu.enabled);

    scan(stray.enabled);
    scan(stray.van_cittert_steps);

    scan(swath.enabled);
    scan(swath.spectrum_width);

    scan(rad.enabled);

    scan(io.ckd);
    scan(io.binning_table);
    scan(io.l1a);
    scan(io.geometry);
    scan(io.l1b);
}

auto SettingsL1B::checkParameters() -> void
{
    // These files must always be present
    checkPresenceOfFile(io.ckd, true);
    checkPresenceOfFile(io.binning_table, true);
    checkPresenceOfFile(io.l1a, true);
    checkPresenceOfFile(io.geometry, true);

    checkFileWritable(io.l1b);
}

} // namespace tango
