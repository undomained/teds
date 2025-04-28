// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing the binning table and functions for
// binning/unbinning different types of data.

#pragma once

#include "eigen.h"

namespace tango {

class BinningTable
{
private:
    // bin_indices[i] is the index of the ith element of the unbinned
    // data. The data could represent, e.g., an image or part of the
    // CKD.
    std::vector<uint32_t> bin_indices {};

public:
    // count_table[i] is the superpixel size of the ith element of
    // binned data
    Eigen::ArrayXd count_table {};

    // Construct a binning table from the binning table name given by
    // user input and the binning table ID. The latter is stored in an
    // L1A file and corresponds to a group name in the binning table
    // file. For example, ID 5 means use the binning table from group
    // Table_5 of the binning table file. If not using a trivial
    // binning table (ID > 1) then detector rows and columns are
    // ignored.
    BinningTable(const int detector_n_rows,
                 const int detector_n_cols,
                 const std::string& binning_table,
                 const int table_id);
    // Size of binned data
    [[nodiscard]] auto nBins() const -> int
    {
        return static_cast<int>(count_table.size());
    }
    // Return a binned index
    [[nodiscard]] auto binIndex(const int idx) const -> int
    {
        return static_cast<int>(bin_indices[idx]);
    }
    // Return a bin (superpixel) size
    [[nodiscard]] auto binSize(const int idx) const -> int
    {
        return static_cast<int>(count_table(idx));
    }
    // Bin an array
    auto bin(const Eigen::ArrayXd& data) const -> Eigen::ArrayXd;
    auto bin(const ArrayXXd& data) const -> Eigen::ArrayXd;
    // Bin a boolean array
    auto bin(const ArrayXb& data) const -> ArrayXb;
    // Bin multiple images in a 2D Eigen array. One row corresponds to
    // one image. If scale_by_bin_size then binned array is divided by
    // the count table (false for normal L1A product).
    auto binMulti(const ArrayXXd& data,
                  const bool scale_by_bin_size) const -> ArrayXXd;

    // Unbin data using the binning table. Output array is the same
    // size as the binning table (bin_indices.size()).
    auto unbin(const Eigen::Ref<const Eigen::ArrayXd> data) const
      -> Eigen::ArrayXd;
};

} // namespace tango
