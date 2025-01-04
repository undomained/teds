// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Digital elevation model (DEM) for computing the surface height for
// target latitudes and longitudes.

#pragma once

#include <string>
#include <vector>

namespace tango {

class DEM
{
private:
    // DEM file
    std::string filename {};
    // DEM grid
    std::vector<double> latitudes {};
    std::vector<double> longitudes {};
    std::vector<int16_t> heights {};
    // Whether the lat/lon and heights have been loaded into memory
    bool in_memory { false };

public:
    DEM(const std::string& filename) : filename { filename } {}
    // Load a portion of the DEM file around a target grid point into
    // memory. Automatically called on first call to getHeight.
    auto loadToMemory(const double lat, const double lon) -> void;
    // Return the height of a target grid point
    [[nodiscard]] auto getHeight(const double lat, const double lon) -> int16_t;
    ~DEM() = default;
};

} // namespace tango
