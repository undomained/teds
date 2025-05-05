// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "nonscat_fwd_model.h"

#include "optic_abs_prop.h"
#include "surface.h"

#include <common/isrf.h>
#include <numbers>

namespace tango {

static auto transmission(const Eigen::ArrayXd& sun_lbl,
                         OpticAbsProp& optics,
                         const Eigen::ArrayXd& albedo,
                         const double mu0,
                         const double muv,
                         Eigen::ArrayXd& rad,
                         Eigen::ArrayXd& dev_tau,
                         Eigen::ArrayXd& dev_alb) -> void
{
    const double mueff { std::abs(1.0 / mu0) + std::abs(1.0 / muv) };
    const double fact { mu0 / std::numbers::pi };
    const Eigen::ArrayXd exp_tot = Eigen::exp(-optics.tau_tot * mueff);
    dev_alb = fact * exp_tot * sun_lbl;
    rad = albedo * dev_alb;
    dev_tau = -mueff * rad;
}

auto nonscatFwdModel(const Eigen::VectorXd& state_vector,
                     const std::vector<std::string>& gas_names,
                     const int n_albedo,
                     const int n_shift,
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
                     Eigen::MatrixXd& jacobian) -> void
{
    const int n_gases { static_cast<int>(gas_names.size()) };
    optics.setOptDepth(gas_names, state_vector);

    transmission(sun_lbl,
                 optics,
                 surface.alb,
                 mu0,
                 muv,
                 rad_lbl,
                 dev_tau_lbl,
                 dev_alb_lbl);

    // Radiance
    rad = isrf.convolve(rad_lbl);

    // Derivative with respect to a scaling of the total optical depth
    for (int i {}; i < n_gases; ++i) {
        const std::string& name { gas_names[i] };
        jacobian.col(i) = isrf.convolve(optics.tau.at(name) * dev_tau_lbl);
    }

    // Albedo coefficients
    jacobian.col(n_gases) = isrf.convolve(dev_alb_lbl);
    for (int i_deriv { 1 }; i_deriv < n_albedo; ++i_deriv) {
        jacobian.col(n_gases + i_deriv) =
          isrf.convolve(dev_alb_lbl * surface.waves.col(i_deriv - 1));
    }

    // Spectral Shift
    if (n_shift > 0) {
        jacobian.col(n_gases + n_albedo) = isrf.convolveDer(rad_lbl);
    }
}

auto speciesDerivLayers(const std::string& gas,
                        const ISRF& isrf,
                        ArrayXXd& tau_alt,
                        Eigen::ArrayXd& dev_tau_lbl,
                        Eigen::MatrixXd& deriv_layers) -> void
{
    for (long int k {}; k < deriv_layers.cols(); ++k) {
        deriv_layers.col(k) = isrf.convolve(tau_alt.col(k) * dev_tau_lbl);
    }
}

} // namespace tango
