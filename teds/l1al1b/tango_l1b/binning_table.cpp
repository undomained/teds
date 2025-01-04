// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "binning_table.h"

#include <netcdf>
#include <numeric>
#include <sstream>

namespace tango {

BinningTable::BinningTable(const int detector_n_rows,
                           const int detector_n_cols,
                           const std::string& binning_table,
                           const int table_id)
{
    if (table_id < 2) {
        // Default binning table is the same size as the unbinned detector
        bin_indices.resize(detector_n_rows * detector_n_cols);
        std::iota(bin_indices.begin(), bin_indices.end(), 0);
        count_table.resize(detector_n_rows * detector_n_cols, 1);
        return;
    }
    const netCDF::NcFile nc { binning_table, netCDF::NcFile::read };
    std::stringstream group_id {};
    group_id << table_id;
    const std::string group_name { "Table_" + group_id.str() };
    const auto n_rows { nc.getDim("row").getSize() };
    const auto n_cols { nc.getDim("column").getSize() };
    const auto n_bins { nc.getGroup(group_name).getDim("bins").getSize() };
    bin_indices.resize(n_rows * n_cols);
    nc.getGroup(group_name).getVar("binning_table").getVar(bin_indices.data());
    count_table.resize(n_bins);
    nc.getGroup(group_name).getVar("count_table").getVar(count_table.data());
}

auto BinningTable::bin(const std::vector<double>& data,
                       std::vector<double>& data_binned) const -> void
{
    std::ranges::fill(data_binned, 0.0);
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        data_binned[bin_indices[i]] += data[i];
    }
    for (int i {}; i < static_cast<int>(data_binned.size()); ++i) {
        data_binned[i] /= count_table[i];
    }
}

auto BinningTable::binUnscaled(std::vector<double>& data) const -> void
{
    const size_t n_full { bin_indices.size() };
    const size_t n_binned { count_table.size() };
    const size_t n_alt { data.size() / n_full };
    std::vector<double> data_binned(n_alt * count_table.size(), 0.0);
#pragma omp parallel for
    for (size_t i_alt = 0; i_alt < n_alt; ++i_alt) {
        for (size_t i {}; i < n_full; ++i) {
            data_binned[i_alt * n_binned + bin_indices[i]] +=
              data[i_alt * n_full + i];
        }
    }
    data = std::move(data_binned);
}

auto BinningTable::bin(std::vector<double>& data) const -> void
{
    std::vector<double> data_binned(count_table.size());
    bin(data, data_binned);
    data = std::move(data_binned);
}

auto BinningTable::bin(std::vector<bool>& data) const -> void
{
    std::vector<bool> data_binned(count_table.size(), false);
    for (int i {}; i < static_cast<int>(data.size()); ++i) {
        data_binned[bin_indices[i]] = data_binned[bin_indices[i]] || data[i];
    }
    data = std::move(data_binned);
}

} // namespace tango
