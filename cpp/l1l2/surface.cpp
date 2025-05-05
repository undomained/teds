// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "surface.h"

#include <algorithm>
#include <numeric>

namespace tango {

Surface::Surface(const int n_albedos, const Eigen::ArrayXd& wave_lbl)
{
    waves.resize(wave_lbl.size(), n_albedos);
    waves.col(0) = (wave_lbl - wave_lbl.mean())
                   / (wave_lbl.maxCoeff() - wave_lbl.minCoeff());
    for (int i { 1 }; i < n_albedos; ++i) {
        waves.col(i) = waves.col(i - 1) * waves.col(0);
    }
    alb.resize(wave_lbl.size());
}

auto Surface::getAlbedoPoly(const Eigen::Ref<const Eigen::ArrayXd> albedo_coeff)
  -> void
{
    alb = albedo_coeff(0);
    for (long int i { 1 }; i < albedo_coeff.size(); ++i) {
        alb += albedo_coeff(i) * waves.col(i - 1);
    }
}

} // namespace tango
