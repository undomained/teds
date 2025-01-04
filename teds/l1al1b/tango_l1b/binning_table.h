// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing the binning table and functions for
// binning/unbinning different types of data.

#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace tango {

class BinningTable
{
private:
    // bin_indices[i] is the index of the ith element of the unbinned
    // data. The data could represent, e.g., an image or part of the
    // CKD.
    std::vector<uint32_t> bin_indices {};
    // count_table[i] is the superpixel size of the ith element of
    // binned data
    std::vector<uint16_t> count_table {};

public:
    // Construct a binning table from the binning table name given by
    // user input and the binning table ID. The latter is stored in an
    // L1A file and corresponds to a group name in the binning table
    // file. For example, ID 5 means use the binning table from group
    // Table_5 of the binning table file.
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
        return static_cast<int>(count_table[idx]);
    }
    // Bin an array and save result in data_binned
    auto bin(const std::vector<double>& data,
             std::vector<double>& data_binned) const -> void;
    // Bin an array and save the result in the same array
    auto bin(std::vector<double>& data) const -> void;
    // Bin a boolean array
    auto bin(std::vector<bool>& data) const -> void;
    // Like bin but don't multiply with the count table
    auto binUnscaled(std::vector<double>& data) const -> void;
    // Unbin data using the binning table. Output array is the same
    // size as the binning table (bin_indices.size()).
    auto unbin(const std::ranges::range auto& data,
               std::vector<double>& data_unbinned) const -> void
    {
        std::ranges::fill(data_unbinned, 0.0);
        for (int i {}; i < static_cast<int>(bin_indices.size()); ++i) {
            data_unbinned[i] = data[bin_indices[i]];
        }
    }
};

} // namespace tango
