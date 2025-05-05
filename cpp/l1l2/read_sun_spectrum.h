// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Read sun spectrum TSIS-1 HSRS

#pragma once

#include <Eigen/Dense>
#include <string>

namespace tango {

auto readSunSpectrum(const std::string& filename,
                     const Eigen::ArrayXd& wave_lbl) -> Eigen::ArrayXd;

} // namespace tango {
