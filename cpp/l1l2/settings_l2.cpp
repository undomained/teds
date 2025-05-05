// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "settings_l2.h"

#include <common/io.h>

namespace tango {

auto SettingsL2::scanKeys() -> void
{
    scan(processing_version);
    scan(alt_beg);
    scan(alt_end);
    scan(compress);

    scan(retrieval.max_iter);
    scan(retrieval.chi2_lim);
    scan(retrieval.n_albedos);
    scan(retrieval.gases);
    scan(retrieval.initial_concentrations);

    scan(atmosphere.n_layers);
    scan(atmosphere.layer_thickness);
    scan(atmosphere.surface_pressure);

    scan(spec_settings.wave_start);
    scan(spec_settings.wave_end);
    scan(spec_settings.wave_extend);
    scan(spec_settings.dwave);

    scan(isrf.tabulated);
    scan(isrf.fwhm);
    scan(isrf.shape);

    scan(io_files.isrf);
    scan(io_files.l1b);
    scan(io_files.atmosphere);
    scan(io_files.afgl);
    scan(io_files.sun_reference);
    scan(io_files.dump_xsec);
    scan(io_files.l2);
    scan(io_files.l2_diag);
}

auto SettingsL2::checkParameters() -> void
{
    // For consistency, gas names should always be uppercase
    for (std::string& gas : retrieval.gases) {
        gas = upper(gas);
    }
}

} // namespace tango
