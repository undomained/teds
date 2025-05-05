// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "binning_table.h"

#include <netcdf>
#include <numeric>

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
        count_table =
          Eigen::ArrayXd::Constant(detector_n_rows * detector_n_cols, 1.0);
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

auto BinningTable::bin(const Eigen::ArrayXd& data) const -> Eigen::ArrayXd
{
    Eigen::ArrayXd data_binned = Eigen::ArrayXd::Zero(count_table.size());
    for (long int i {}; i < data.size(); ++i) {
        data_binned(bin_indices[i]) += data(i);
    }
    data_binned /= count_table;
    return data_binned;
}

auto BinningTable::bin(const ArrayXXd& data) const -> Eigen::ArrayXd
{
    Eigen::ArrayXd data_binned = Eigen::ArrayXd::Zero(count_table.size());
    for (long int i {}; i < data.size(); ++i) {
        data_binned(bin_indices[i]) += data.data()[i];
    }
    data_binned /= count_table;
    return data_binned;
}

auto BinningTable::binMulti(const ArrayXXd& data,
                            const bool scale_by_bin_size) const -> ArrayXXd
{
    ArrayXXd data_binned = ArrayXXd::Zero(data.rows(), count_table.size());
    for (long int i_row {}; i_row < data.rows(); ++i_row) {
        for (long int i {}; i < data.cols(); ++i) {
            data_binned(i_row, bin_indices[i]) += data(i_row, i);
        }
        if (scale_by_bin_size) {
            data_binned.row(i_row) /= count_table;
        }
    }
    return data_binned;
}

auto BinningTable::bin(const ArrayXb& data) const -> ArrayXb
{
    ArrayXb data_binned = ArrayXb::Constant(count_table.size(), false);
    for (long int i {}; i < data.size(); ++i) {
        data_binned[bin_indices[i]] = data_binned[bin_indices[i]] || data[i];
    }
    return data_binned;
}

auto BinningTable::unbin(const Eigen::Ref<const Eigen::ArrayXd> data) const
  -> Eigen::ArrayXd
{
    Eigen::ArrayXd data_unbinned = Eigen::ArrayXd::Zero(bin_indices.size());
    for (size_t i {}; i < bin_indices.size(); ++i) {
        data_unbinned(i) = data(bin_indices[i]);
    }
    return data_unbinned;
}

} // namespace tango
