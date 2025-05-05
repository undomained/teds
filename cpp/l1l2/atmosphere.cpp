// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "atmosphere.h"

#include <common/constants.h>
#include <common/linear_spline.h>

#include <fstream>
#include <sstream>

namespace tango {

// Read standard atmosphere data into a 2D array
static auto arraysFromText(const std::string& filename) -> Eigen::ArrayXXd
{
    std::ifstream in { filename };
    std::string line {};
    std::vector<double> buf {};
    int n_cols {};
    while (std::getline(in, line)) {
        if (line[0] == '!') {
            continue;
        }
        std::stringstream ss { line };
        double value {};
        while (ss >> value) {
            buf.push_back(value);
        }
        if (n_cols == 0) {
            n_cols = buf.size();
        }
    }
    auto raw_it { buf.cbegin() };
    Eigen::ArrayXXd data(buf.size() / n_cols, n_cols);
    for (auto& x : data.transpose().reshaped()) {
        x = *raw_it;
        ++raw_it;
    }
    return data;
}

Atmosphere::Atmosphere(const int n_layers,
                       const double layer_thickness,
                       const double surface_pressure,
                       const std::string& afgl_filename)
{
    z_lay = Eigen::ArrayXd::LinSpaced(n_layers, n_layers - 0.5, 0.5)
            * layer_thickness;
    z_lev =
      Eigen::ArrayXd::LinSpaced(n_layers + 1, n_layers, 0.0) * layer_thickness;

    Eigen::ArrayXXd afgl_data = arraysFromText(afgl_filename);

    // Height [km] -> [m]
    Eigen::ArrayXd zalt_in = afgl_data.col(0) * 1e3;
    // Pressure [Pa]
    Eigen::ArrayXd press_in = afgl_data.col(1) * 1e2;
    // Temperature [K]
    Eigen::ArrayXd temp_in = afgl_data.col(2);
    // Air number density [#/cm^3]
    Eigen::ArrayXd air_in = afgl_data.col(3);
    // O3 number density -> mole fraction [-]
    Eigen::ArrayXd o3_in = afgl_data.col(4) / air_in;
    // O2 number density -> mole fraction [-]
    Eigen::ArrayXd o2_in = afgl_data.col(5) / air_in;
    // H2O number density -> mole fraction [-]
    Eigen::ArrayXd h2o_in = afgl_data.col(6) / air_in;
    // CO2 number density -> mole fraction [-]
    Eigen::ArrayXd co2_in = afgl_data.col(7) / air_in;
    // NO2 number density -> mole fraction [-]
    Eigen::ArrayXd no2_in = afgl_data.col(8) / air_in;

    // Truncate or extrapolate the AFGL profile depending on surface
    // pressure.
    if (surface_pressure - press_in(Eigen::last) > 1e-10) {
        throw std::runtime_error { "not implemented yet" };
    } else if (press_in(Eigen::last) - surface_pressure > 1e-10) {
        throw std::runtime_error { "not implemented yet" };
    }

    // Interpolate quantities to output layers

    LinearSpline spline { zalt_in, press_in };
    Eigen::ArrayXd plev = spline.eval(z_lev);
    Eigen::ArrayXd play = spline.eval(z_lay);

    constexpr double sp { atm::na / (atm::mdryair * atm::g0) };
    Eigen::ArrayXd vc_air = sp * plev;
    air = vc_air(Eigen::seqN(1, vc_air.size() - 1))
          - vc_air(Eigen::seqN(0, vc_air.size() - 1));
    air(0) = vc_air(0);
    air_sum = air.sum();

    spline = { zalt_in, o3_in };
    gases["O3"] = spline.eval(z_lay) * air;

    spline = { zalt_in, h2o_in };
    gases["H2O"] = spline.eval(z_lay) * air;

    spline = { zalt_in, co2_in };
    gases["CO2"] = spline.eval(z_lay) * air;

    spline = { zalt_in, no2_in };
    gases["NO2"] = spline.eval(z_lay) * air;

    gases["O2"] = Eigen::ArrayXd::Constant(air.size(), atm::xo2) * air;

    gases["CH4"] = Eigen::ArrayXd::Constant(air.size(), atm::xch4) * air;
}

} // namespace tango
