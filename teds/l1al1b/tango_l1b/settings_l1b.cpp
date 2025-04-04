// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "settings_l1b.h"

#include "io.h"

namespace tango {

auto SettingsL1B::scanKeys() -> void
{
    scan(processing_version);
    scan(cal_level);
    scan(alt_beg);
    scan(alt_end);
    scan(bin_spectra);
    scan(compress);

    scan(dark.enabled);

    scan(noise.enabled);

    scan(nonlin.enabled);

    scan(prnu.enabled);

    scan(stray.van_cittert_steps);

    scan(swath.b_spline_order);
    scan(swath.geolocation);

    scan(rad.enabled);

    scan(io_files.ckd);
    scan(io_files.binning_table);
    scan(io_files.l1a);
    scan(io_files.l1b);
    scan(io_files.dem);
    scan(io_files.geometry);
}

auto SettingsL1B::checkParameters() -> void
{
    // These files must always be present
    checkPresenceOfFile(io_files.ckd, true);
    checkPresenceOfFile(io_files.binning_table, true);
    checkPresenceOfFile(io_files.l1a, true);
    checkPresenceOfFile(io_files.geometry, !swath.geolocation);

    checkFileWritable(io_files.l1b);
}

} // namespace tango
