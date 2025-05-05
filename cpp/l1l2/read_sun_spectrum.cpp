// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "read_sun_spectrum.h"

#include <common/constants.h>
#include <common/linear_spline.h>

#include <netcdf>

namespace tango {

auto readSunSpectrum(const std::string& filename,
                     const Eigen::ArrayXd& wave_lbl) -> Eigen::ArrayXd
{
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    const auto n_wavelengths { nc.getDim("wavelength").getSize() };
    Eigen::ArrayXd sun_wavelengths(n_wavelengths);
    Eigen::ArrayXd sun_spectrum(n_wavelengths);
    nc.getVar("Vacuum Wavelength").getVar(sun_wavelengths.data());
    nc.getVar("SSI").getVar(sun_spectrum.data());
    sun_spectrum *= sun_wavelengths * 1e-9 / (atm::hplanck * atm::clight);
    const LinearSpline spline { sun_wavelengths, sun_spectrum };
    return spline.eval(wave_lbl);
}

} // namespace tango {
