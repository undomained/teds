// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "optic_abs_prop.h"

#include "atmosphere.h"

#include <netcdf>

namespace tango {

OpticAbsProp::OpticAbsProp(const std::string& xsec_filename)
{
    const netCDF::NcFile nc { xsec_filename, netCDF::NcFile::read };
    const auto n_layers { nc.getDim("layer").getSize() };
    const auto n_wavelengths { nc.getDim("wavelength").getSize() };
    for (auto& [name, group] : nc.getGroups()) {
        // Cross sections are given in cm^2, atmospheric densities
        // in m^2
        constexpr double conv_f { 1e-4 };
        xsec[name].resize(n_wavelengths, n_layers);
        group.getVar("xsec").getVar(xsec[name].data());
        xsec[name] *= conv_f;
        tau_alt[name].resize(n_wavelengths, n_layers);
        tau_ref[name].resize(n_wavelengths);
        tau[name].resize(n_wavelengths);
    }
    tau_tot.resize(n_wavelengths);
}

auto OpticAbsProp::checkWavelengthGrid(const Eigen::ArrayXd& wave) const -> void
{
    if (wave.size() != tau_tot.rows()) {
        throw std::runtime_error("absorption cross-sections have been "
                                 "generated for a different wavelength grid");
    }
}

auto OpticAbsProp::setOptDepthRef(const Atmosphere& atm,
                                  const std::vector<std::string>& species_names)
  -> void
{
    for (int i_spec {}; i_spec < static_cast<int>(species_names.size());
         ++i_spec) {
        const std::string& name { species_names[i_spec] };
        tau_ref.at(name) =
          (xsec.at(name).rowwise() * atm.gases.at(name).transpose())
            .rowwise()
            .sum();
    }
}

auto OpticAbsProp::setOptDepthLayers(
  const std::map<std::string, Eigen::ArrayXd>& gases,
  const std::vector<std::string>& species_names) -> void
{
    for (int i_spec {}; i_spec < static_cast<int>(species_names.size());
         ++i_spec) {
        const std::string& name { species_names[i_spec] };
        tau_alt.at(name) = xsec.at(name).rowwise() * gases.at(name).transpose();
    }
}

auto OpticAbsProp::setOptDepth(const std::vector<std::string>& species_names,
                               const Eigen::VectorXd& state_vector) -> void
{
    for (int i_spec {}; i_spec < static_cast<int>(species_names.size());
         ++i_spec) {
        const std::string& name { species_names[i_spec] };
        tau.at(name) = state_vector(i_spec) * tau_ref.at(name);
        if (i_spec == 0) {
            tau_tot = tau.at(name);
        } else {
            tau_tot += tau.at(name);
        }
    }
}

} // namespace tango {
