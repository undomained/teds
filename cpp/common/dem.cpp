// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dem.h"

#include "constants.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <netcdf>
#include <numbers>
#include <numeric>

namespace tango {

// Find index of the first element larger than a value in a sorted array
static auto findIdx(const std::vector<double>& values,
                    const double value) -> size_t
{
    return static_cast<size_t>(std::distance(
      values.begin(), std::ranges::find_if(values, [value](const double x) {
          return x > value;
      })));
}

auto DEM::loadToMemory(const double lat, const double lon) -> void
{
    if (in_memory) {
        return;
    }
    // Read in the DEM grid
    const netCDF::NcFile nc { filename, netCDF::NcFile::read };
    const auto n_lat { nc.getDim("lat").getSize() };
    const auto n_lon { nc.getDim("lon").getSize() };
    latitudes.resize(n_lat);
    longitudes.resize(n_lon);
    nc.getVar("lat").getVar(latitudes.data());
    nc.getVar("lon").getVar(longitudes.data());
    for (double& val : latitudes) {
        val *= math::deg_to_rad;
    }
    for (double& val : longitudes) {
        val *= math::deg_to_rad;
    }
    // Construct a rectangle around the target lat/lon point
    constexpr double lat_margin { 2.0 * math::deg_to_rad };
    const double lon_margin { 2.0 * math::deg_to_rad / std::cos(lat) };
    const double lat_beg { lat - lat_margin };
    const double lat_end { lat + lat_margin };
    const double lon_beg { lon - lon_margin };
    const double lon_end { lon + lon_margin };
    // Convert the rectangle corners from lat/lon to indices of DEM
    // file lat/lon arrays.
    const size_t idx_lat_beg { findIdx(latitudes, lat_beg) };
    const size_t idx_lat_end { findIdx(latitudes, lat_end) };
    const size_t idx_lon_beg { findIdx(longitudes, lon_beg) };
    const size_t idx_lon_end { findIdx(longitudes, lon_end) };
    assert(idx_lat_beg < latitudes.size() && idx_lat_end < latitudes.size()
           && idx_lon_beg < longitudes.size()
           && idx_lon_end < longitudes.size());
    // Retrieve all heights for the target area
    const std::vector<size_t> starts { idx_lat_beg, idx_lon_beg };
    const std::vector<size_t> counts { idx_lat_end - idx_lat_beg,
                                       idx_lon_end - idx_lon_beg };
    heights.resize(counts[0] * counts[1]);
    nc.getVar("height").getVar(starts, counts, heights.data());
    // Set water height to zero
    for (int16_t& height : heights) {
        if (height < 0) {
            height = 0;
        }
    }
    // Trim the lat/lon grid
    latitudes = { latitudes.begin() + idx_lat_beg,
                  latitudes.begin() + idx_lat_beg + counts[0] };
    longitudes = { longitudes.begin() + idx_lon_beg,
                   longitudes.begin() + idx_lon_beg + counts[1] };
    in_memory = true;
}

[[nodiscard]] auto DEM::getHeight(const double lat, const double lon) -> int16_t
{
    if (!in_memory) {
        // This should be called only once for the entire run
#pragma omp critical
        loadToMemory(lat, lon);
    }
    size_t idx_lat { findIdx(latitudes, lat) };
    size_t idx_lon { findIdx(longitudes, lon) };
    return heights[idx_lat * longitudes.size() + idx_lon];
}

} // namespace tango
