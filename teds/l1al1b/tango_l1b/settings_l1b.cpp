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
    scan(alt_beg);
    scan(alt_end);
    scan(bin_spectra);

    scan(dark.enabled);

    scan(noise.enabled);

    scan(nonlin.enabled);

    scan(prnu.enabled);

    scan(stray.van_cittert_steps);

    scan(swath.b_spline_order);
    scan(swath.exact_drawing);
    scan(swath.geolocation);

    scan(rad.enabled);

    scan(io.ckd);
    scan(io.binning_table);
    scan(io.l1a);
    scan(io.l1b);
    scan(io.dem);
    scan(io.geometry);
}

auto SettingsL1B::checkParameters() -> void
{
    // These files must always be present
    checkPresenceOfFile(io.ckd, true);
    checkPresenceOfFile(io.binning_table, true);
    checkPresenceOfFile(io.l1a, true);
    checkPresenceOfFile(io.geometry, !swath.geolocation);

    checkFileWritable(io.l1b);

    if (swath.exact_drawing) {
        noise.artificial_scaling =
          1 / std::sqrt(static_cast<double>(bin_spectra));
        bin_spectra = 1;
    }
}

} // namespace tango
