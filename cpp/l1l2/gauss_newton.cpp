// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "gauss_newton.h"

#include "atmosphere.h"
#include "nonscat_fwd_model.h"
#include "optic_abs_prop.h"
#include "ref_profiles.h"
#include "settings_l2.h"
#include "surface.h"

#include <common/isrf.h>
#include <common/l2.h>
#include <numeric>

namespace tango {

auto gaussNewton(const int i_gp,
                 const SettingsL2& settings,
                 const RefProfiles& ref_profiles,
                 const Atmosphere& atm_orig,
                 const Eigen::ArrayXd& wave_lbl,
                 const Eigen::ArrayXd& sun_lbl,
                 const Eigen::Ref<const Eigen::ArrayXd> spectrum,
                 const Eigen::Ref<const Eigen::ArrayXd> S_meas,
                 const double mu0,
                 const double muv,
                 const ISRF& isrf,
                 OpticAbsProp& optics,
                 const Eigen::ArrayXd& sun,
                 L2& l2) -> void
{
    // A copy of the standard atmosphere is use for storing the
    // retrieved gas columns.
    Atmosphere atm { atm_orig };

    // Gases to be retrieved
    std::vector<std::string> gas_names {};
    for (auto it { ref_profiles.gases.begin() }; it != ref_profiles.gases.end();
         ++it) {
        gas_names.push_back(it->first);
    }
    const int n_gases { static_cast<int>(gas_names.size()) };
    const int n_albedos { settings.retrieval.n_albedos };

    // State vector, to be populated with gas concentrations and
    // albedo coefficients.
    Eigen::VectorXd state_vector(n_gases + n_albedos);

    // Start by adding gas scaling initial guesses to the state vector
    for (int i {}; i < n_gases; ++i) {
        state_vector(i) =
          ref_profiles.initial.at(gas_names[i])
          / (ref_profiles.gas_sums.at(gas_names[i]) / atm.air_sum);
    }

    // Next add albedo coefficient initial guesses. The first entries
    // of the state vector correspond to gases.
    state_vector(Eigen::seqN(n_gases, n_albedos)).array() = 0.0;
    // Derive first guess albedo from maximum reflectance
    Eigen::Index idx_max {};
    const double max_rad { spectrum.maxCoeff(&idx_max) };
    state_vector(n_gases) = max_rad / sun(idx_max) * std::numbers::pi / mu0;

    // Initialize Gauss-Newton work variables
    const int n_dof { static_cast<int>(spectrum.size() - state_vector.size()) };
    double chi2_prev { 1e10 };
    Surface surface { n_albedos, wave_lbl };
    Eigen::VectorXd rad(spectrum.size());
    Eigen::MatrixXd jacobian(spectrum.size(), state_vector.size());
    Eigen::ArrayXd rad_lbl {};
    Eigen::ArrayXd dev_tau_lbl {};
    Eigen::ArrayXd dev_alb_lbl {};
    Eigen::ArrayXd precisions {};
    Eigen::MatrixXd gains {};

    for (int iterations {}; iterations < settings.retrieval.max_iter;
         ++iterations) {
        for (int i_gas {}; i_gas < n_gases; ++i_gas) {
            atm.gases.at(gas_names[i_gas]) =
              state_vector[i_gas] * ref_profiles.gases.at(gas_names[i_gas]);
        }

        surface.getAlbedoPoly(state_vector(Eigen::seqN(n_gases, n_albedos)));

        nonscatFwdModel(state_vector,
                        gas_names,
                        n_albedos,
                        isrf,
                        sun_lbl,
                        atm,
                        surface,
                        mu0,
                        muv,
                        optics,
                        rad_lbl,
                        dev_tau_lbl,
                        dev_alb_lbl,
                        rad,
                        jacobian);

        const Eigen::VectorXd y_diff = spectrum.matrix() - rad;

        // Calculate least square solution
        const Eigen::MatrixXd Syinv = (1.0 / S_meas).matrix().asDiagonal();
        const Eigen::MatrixXd JTJ = jacobian.transpose() * (Syinv * jacobian);
        const Eigen::MatrixXd Sx =
          JTJ.llt().solve(Eigen::MatrixXd::Identity(JTJ.rows(), JTJ.cols()));
        gains = Sx * (jacobian.transpose() * Syinv);
        const Eigen::VectorXd state_delta = gains * y_diff;

        state_vector += state_delta;

        l2.chi2(i_gp) = (y_diff.array().square() / S_meas).sum() / n_dof;
        l2.iterations(i_gp) = iterations + 1;

        const bool converged { chi2_prev - l2.chi2(i_gp)
                               < settings.retrieval.chi2_lim };
        if (converged || iterations == settings.retrieval.max_iter - 1) {
            precisions = Sx.diagonal().array().sqrt();
        }
        if (converged) {
            l2.converged(i_gp) = true;
            break;
        }
        chi2_prev = l2.chi2(i_gp);
    }

    // Update gas concentrations
    for (int i_gas {}; i_gas < n_gases; ++i_gas) {
        atm.gases.at(gas_names[i_gas]) =
          state_vector[i_gas] * ref_profiles.gases.at(gas_names[i_gas]);
    }

    // Column mixing ratios, precisions, and albedo values
    for (int i_gas {}; i_gas < n_gases; ++i_gas) {
        const std::string& gas { gas_names[i_gas] };
        l2.mixing_ratios.at(gas)(i_gp) = atm.gases.at(gas).sum() / atm.air_sum;
        l2.precisions.at(gas)(i_gp) =
          l2.mixing_ratios.at(gas)(i_gp) * precisions[i_gas];
        l2.gains.at(gas).row(i_gp) = gains.row(i_gas)
                                     * l2.mixing_ratios.at(gas)(i_gp)
                                     * (gas == "CH4" ? 1e9 : 1e6);
    }
    l2.albedo0(i_gp) = state_vector(n_gases);

    // Column averaging kernels
    optics.setOptDepthLayers(atm.gases, gas_names);
    Eigen::MatrixXd deriv_gas_layers(spectrum.size(), atm.z_lay.size());
    for (int i_gas {}; i_gas < n_gases; ++i_gas) {
        const std::string& gas { gas_names[i_gas] };
        speciesDerivLayers(
          gas, isrf, optics.tau_alt.at(gas), dev_tau_lbl, deriv_gas_layers);
        l2.col_avg_kernels.at(gas).row(i_gp) =
          (gains.row(i_gas) * deriv_gas_layers).array()
          * ref_profiles.gas_sums.at(gas)
          / ref_profiles.gases.at(gas).transpose();
    }
}

} // namespace tango
