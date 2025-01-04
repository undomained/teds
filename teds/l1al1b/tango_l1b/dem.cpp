// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "dem.h"

#include "constants.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <hdf5.h>
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

// Read part or all of an HDF5 dataset
template <typename T>
static auto readHDF5Dset(const hid_t& file,
                         const std::string& dataset_name,
                         const std::vector<size_t> starts,
                         const std::vector<size_t> counts,
                         std::vector<T>& dest) -> void
{
    const hid_t dataset { H5Dopen(file, dataset_name.c_str(), H5P_DEFAULT) };
    const hid_t dataspace { H5Dget_space(dataset) };
    const auto n_dim { H5Sget_simple_extent_ndims(dataspace) };
    std::vector<hsize_t> h_starts(n_dim);
    std::vector<hsize_t> h_counts(n_dim);
    if (starts.empty()) {
        std::ranges::fill(h_starts, 0);
        H5Sget_simple_extent_dims(dataspace, h_counts.data(), NULL);
    } else {
        std::copy(starts.cbegin(), starts.cend(), h_starts.begin());
        std::copy(counts.cbegin(), counts.cend(), h_counts.begin());
    }
    dest.resize(std::accumulate(
      h_counts.cbegin(), h_counts.cend(), 1, std::multiplies<hsize_t>()));
    H5Sselect_hyperslab(
      dataspace, H5S_SELECT_SET, h_starts.data(), NULL, h_counts.data(), NULL);
    const auto memspace { H5Screate_simple(n_dim, h_counts.data(), NULL) };
    H5Dread(dataset,
            H5Dget_type(dataset),
            memspace,
            dataspace,
            H5P_DEFAULT,
            dest.data());
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}

auto DEM::loadToMemory(const double lat, const double lon) -> void
{
    if (in_memory) {
        return;
    }
    // Read in the DEM grid
    const hid_t file { H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT) };
    readHDF5Dset(file, "/lat", {}, {}, latitudes);
    readHDF5Dset(file, "/lon", {}, {}, longitudes);
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
    const hsize_t idx_lat_beg { static_cast<hsize_t>(
      findIdx(latitudes, lat_beg)) };
    const hsize_t idx_lat_end { static_cast<hsize_t>(
      findIdx(latitudes, lat_end)) };
    const hsize_t idx_lon_beg { static_cast<hsize_t>(
      findIdx(longitudes, lon_beg)) };
    const hsize_t idx_lon_end { static_cast<hsize_t>(
      findIdx(longitudes, lon_end)) };
    assert(idx_lat_beg < latitudes.size() && idx_lat_end < latitudes.size()
           && idx_lon_beg < longitudes.size()
           && idx_lon_end < longitudes.size());
    // Retrieve all heights for the target area
    const std::vector<size_t> starts { idx_lat_beg, idx_lon_beg };
    const std::vector<size_t> counts { idx_lat_end - idx_lat_beg,
                                       idx_lon_end - idx_lon_beg };
    readHDF5Dset(file, "/height", starts, counts, heights);
    H5Fclose(file);
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
