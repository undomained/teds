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
    scan(alt_beg);
    scan(alt_end);

    scan(detector.binning_table_id);
    scan(detector.exposure_time);
    scan(detector.nr_coadditions);

    scan(optimal_coadd.enabled);
    scan(optimal_coadd.FMC);
    scan(optimal_coadd.full_well);
    scan(optimal_coadd.t_dead);
    scan(optimal_coadd.f_sat);

    scan(isrf.enabled);
    scan(isrf.fwhm_gauss);
    scan(isrf.shape);
    scan(isrf.in_memory);

    scan(dark.enabled);

    scan(noise.enabled);
    scan(noise.seed);

    scan(nonlin.enabled);

    scan(prnu.enabled);

    scan(stray.enabled);

    scan(swath.b_spline_order);
    scan(swath.exact_drawing);

    scan(rad.enabled);

    scan(io.ckd);
    scan(io.binning_table);
    scan(io.l1a);
    scan(io.sgm);
    scan(io.geometry);
    scan(io.navigation);
}

auto SettingsIM::checkParameters() -> void
{
    // These files must always be present
    checkPresenceOfFile(io.ckd, true);
    checkPresenceOfFile(io.binning_table, true);
    checkPresenceOfFile(io.sgm, true);
    checkPresenceOfFile(io.geometry, true);

    checkFileWritable(io.l1a);
}

} // namespace tango
