// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class representing a model atmosphere
#pragma once

#include <Eigen/Dense>
#include <map>
#include <string>

namespace tango {

struct Atmosphere
{
    Atmosphere(const int n_layers,
               const double layer_thickness,
               const double surface_pressure,
               const std::string& afgl_filename);

    // Atmospheric layers
    Eigen::ArrayXd z_lay {};
    // Atmospheric layer boundaries
    Eigen::ArrayXd z_lev {};

    // Air and trace gas columns
    Eigen::ArrayXd air {};
    double air_sum {};
    std::map<std::string, Eigen::ArrayXd> gases {};
};

} // namespace tango
