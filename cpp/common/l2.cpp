// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "l2.h"

#include "geometry.h"
#include "time.h"

#include <netcdf>
#include <yaml-cpp/yaml.h>

namespace tango {

L2::L2(const int n_alt,
       const int n_act,
       const int n_wave,
       const int n_lay,
       const std::vector<std::string>& gas_names)
{
    // Number of ground points
    const int n_gp { n_alt * n_act };
    chi2.resize(n_gp);
    converged.resize(n_gp);
    iterations.resize(n_gp);
    albedo0.resize(n_gp);
    for (const auto& gas_name : gas_names) {
        mixing_ratios[gas_name] = Eigen::ArrayXd(n_gp);
        precisions[gas_name] = Eigen::ArrayXd(n_gp);
        gains[gas_name] = ArrayXXd(n_gp, n_wave);
        col_avg_kernels[gas_name] = ArrayXXd(n_gp, n_lay);
    }
    for (const auto& gas_name : { "CO2", "CH4" }) {
        proxys[gas_name] = Eigen::ArrayXd(n_gp);
        proxy_precisions[gas_name] = Eigen::ArrayXd(n_gp);
    }
    spec_shift.resize(n_gp);
    spec_squeeze.resize(n_gp);
}

} // namespace tango
