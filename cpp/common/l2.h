// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Level 2 product definition

#pragma once

#include "eigen.h"

#include <map>
#include <string>
#include <vector>

namespace tango {

class Geometry;

struct L2
{
    Eigen::ArrayXd chi2 {};
    ArrayXb converged {};
    Eigen::ArrayXi iterations {};
    Eigen::ArrayXd albedo0 {};
    std::map<std::string, Eigen::ArrayXd> mixing_ratios {};
    std::map<std::string, Eigen::ArrayXd> precisions {};
    std::map<std::string, ArrayXXd> gains {};
    std::map<std::string, ArrayXXd> col_avg_kernels {};
    std::map<std::string, ArrayXXd> proxys {};
    std::map<std::string, Eigen::ArrayXd> proxy_precisions {};
    Eigen::ArrayXd spec_shift {};
    Eigen::ArrayXd spec_squeeze {};

    L2(const int n_alt,
       const int n_act,
       const int n_wave,
       const int n_lay,
       const std::vector<std::string>& gas_names);
};

} // namespace tango
