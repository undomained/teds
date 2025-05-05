// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include <Eigen/Dense>

namespace tango {

struct Surface
{
    Surface(const int n_albedos, const Eigen::ArrayXd& wave_lbl);

    Eigen::ArrayXXd waves {};
    Eigen::ArrayXd alb {};

    auto getAlbedoPoly(const Eigen::Ref<const Eigen::ArrayXd> albedo_coeff)
      -> void;
};

} // namespace tango
