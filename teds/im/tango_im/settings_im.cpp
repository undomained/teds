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

    scan(io_files.ckd);
    scan(io_files.binning_table);
    scan(io_files.l1a);
    scan(io_files.sgm);
    scan(io_files.geometry);
    scan(io_files.navigation);
}

auto SettingsIM::checkParameters() -> void
{
    // These files must always be present
    checkPresenceOfFile(io_files.ckd, true);
    checkPresenceOfFile(io_files.binning_table, true);
    checkPresenceOfFile(io_files.sgm, true);
    checkPresenceOfFile(io_files.geometry, true);

    checkFileWritable(io_files.l1a);

    if (swath.exact_drawing) {
        noise.artificial_scaling =
          1.0 / std::sqrt(static_cast<double>(detector.binning_table_id));
        detector.binning_table_id = 1;
    }
}

} // namespace tango
