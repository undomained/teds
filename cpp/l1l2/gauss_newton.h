// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Level 2 retrieval based on the Gauss-Newton algorithm.

#pragma once

#include <Eigen/Dense>

namespace tango {

class Atmosphere;
class ISRF;
class L2;
class OpticAbsProp;
class RefProfiles;
class SettingsL2;

// Parameters
// ----------
// i_gp
//     Target ground point (ALT and ACT position)
// settings
//     Retrieval parameters such as number of iterations and chi2
//     convergence criterion
// ref_profiles
//     Model gas columns
// atm
//     Retrieved gas columns. Starts with a copy of the model atmosphere.
// wave_lbl
//     Line-by-line (LBL) wavelength grid
// sun_lbl
//     Solar irradiance spectrum on the LBL grid
// spectrum
//     L1B radiance spectrum corresponding to ground point i_gp
// S_meas
//     Corresponding variance vector
// mu0
//     cos(SZA)
// muv
//     cos(VZA)
// ISRF
//     ISRF class instance for convolutions
// optics
//     Absorption cross-sections and optical thicknesses
// sun
//     ISRF-convolved solar irradiance spectrum
// l2
//     Level 2 product, incrementally filled with each call to this
//     function
auto gaussNewton(const int i_gp,
                 const SettingsL2& settings,
                 const RefProfiles& ref_profiles,
                 const Atmosphere& atm,
                 const Eigen::ArrayXd& wave_lbl,
                 const Eigen::ArrayXd& sun_lbl,
                 const Eigen::Ref<const Eigen::ArrayXd> spectrum,
                 const Eigen::Ref<const Eigen::ArrayXd> S_meas,
                 const double mu0,
                 const double muv,
                 const ISRF& isrf,
                 OpticAbsProp& optics,
                 const Eigen::ArrayXd& sun,
                 L2& l2) -> void;

} // namespace tango
