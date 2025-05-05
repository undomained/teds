// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Optical absorption properties

#pragma once

#include <common/eigen.h>
#include <map>
#include <string>
#include <vector>

namespace tango {

class Atmosphere;

struct OpticAbsProp
{
    OpticAbsProp(const std::string& xsec_filename);

    // Absorption cross-sections, cm^2
    std::map<std::string, ArrayXXd> xsec {};
    // Optical thickness per species (trace gas) per atmosphere layer
    std::map<std::string, ArrayXXd> tau_alt {};
    // Optical thickness per species
    std::map<std::string, Eigen::ArrayXd> tau {};
    // Reference optical thickness profile per species
    std::map<std::string, Eigen::ArrayXd> tau_ref {};
    // Total optical thickness
    Eigen::ArrayXd tau_tot {};

    // Check if the current wavelength grid has the same length as
    // that of the cross-sections. If then the cross-section file
    // needs to be regenerated using the Python version of the L2
    // processor. There is currently no HAPI-like API for C++.
    auto checkWavelengthGrid(const Eigen::ArrayXd& wave) const -> void;

    // Compute reference optical thicknesses from cross-sections and
    // reference gas columns
    auto setOptDepthRef(const Atmosphere& atm,
                        const std::vector<std::string>& species_names) -> void;

    // Compute optical thicknesses from cross-sections and gas columns
    // per layer
    auto setOptDepthLayers(const std::map<std::string, Eigen::ArrayXd>& gases,
                           const std::vector<std::string>& species_names)
      -> void;

    // Update optical thickness of each species from the Gauss-Newton
    // state vector and the reference optical thickness.
    auto setOptDepth(const std::vector<std::string>& species_names,
                     const Eigen::VectorXd& state_vector) -> void;
};

} // namespace tango {
