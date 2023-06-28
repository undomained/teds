// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include <vector>

class PackedMatrix
{
public:
    // See if these can be encapsulated
    int n_rows {};
    int n_cols {};
    std::vector<int> row_locations {};
    std::vector<int> global_col_indices {};
    std::vector<int> row_sizes {};
    std::vector<double> data {};

    auto init(const int n_rows, const int n_cols) -> void;

    // Y =  X A^T where X is self
    auto multRightTransposed(const std::vector<double>& A,
                             std::vector<double>& Y) const -> void;

    // Y =  X A^T where X is self
    auto multRightTransposed(const PackedMatrix& A,
                             std::vector<double>& Y) const -> void;

    // Y = A X^T where X is self
    auto multLeftTransposed(const std::vector<double>& A,
                            std::vector<double>& Y) const -> void;

    // Y^T = A X^T where X is self
    auto multLeftTransposedT(const std::vector<double>& A,
                             std::vector<double>& Y) const -> void;

    auto transpose(const PackedMatrix& X) -> void;

    auto fullMatrix(std::vector<double>& full_matrix) const -> void;
};
