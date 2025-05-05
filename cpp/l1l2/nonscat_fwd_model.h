// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Non-scattering forward model

#pragma once

#include <common/eigen.h>
#include <map>
#include <string>
#include <vector>

namespace tango {

class Atmosphere;
class ISRF;
class OpticAbsProp;
class Surface;

auto nonscatFwdModel(const Eigen::VectorXd& state_vector,
                     const std::vector<std::string>& gas_names,
                     const int n_albedo,
                     const ISRF& isrf,
                     const Eigen::ArrayXd& sun_lbl,
                     const Atmosphere& atm,
                     const Surface& surface,
                     const double mu0,
                     const double muv,
                     OpticAbsProp& optics,
                     Eigen::ArrayXd& rad_lbl,
                     Eigen::ArrayXd& dev_tau_lbl,
                     Eigen::ArrayXd& dev_alb_lbl,
                     Eigen::VectorXd& rad,
                     Eigen::MatrixXd& jacobian) -> void;

// Derivative with respect to a scaling of the layer optical depth
auto speciesDerivLayers(const std::string& gas,
                        const ISRF& isrf,
                        ArrayXXd& tau_alt,
                        Eigen::ArrayXd& dev_tau_lbl,
                        Eigen::MatrixXd& deriv_layers) -> void;

} // namespace tango
