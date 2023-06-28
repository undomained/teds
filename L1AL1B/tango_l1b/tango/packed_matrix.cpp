// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "packed_matrix.h"

#include <cassert>

constexpr int i_one { 1 };

extern "C" {
    auto ddot_(const int*,
               const double*,
               const int*,
               const double*,
               const int*) -> double;
}

auto PackedMatrix::init(const int n_rows, const int n_cols) -> void
{
    this->n_rows = n_rows;
    this->n_cols = n_cols;
    row_locations.resize(n_rows);
    global_col_indices.resize(n_rows, -1);
    row_sizes.resize(n_rows, 0);
}

auto PackedMatrix::multRightTransposed(const std::vector<double>& A,
                                       std::vector<double>& Y) const -> void
{
    const int A_cols { n_cols };
    const int A_rows { static_cast<int>(A.size()) / A_cols };
    for (int i {}; i < n_rows; ++i) {
        const int row_idx_start { row_locations[i] };
        const int row_idx { i * A_rows };
        for (int j {}; j < A_rows; ++j) {
            const int col_idx_start { j * A_cols + global_col_indices[i] };
            double sum {};
            for (int k {}; k < row_sizes[i]; ++k) {
                sum += data[row_idx_start + k]
                       * A[col_idx_start + k];
            }
            Y[row_idx + j] = sum;
        }
        // for (int j {}; j < A_rows; ++j) {
        //     double sum {};
        //     for (int k {}; k < row_sizes[i]; ++k) {
        //         sum += data[row_locations[i] + k]
        //                * A[j * A_cols + global_col_indices[i] + k];
        //     }
        //     Y[i * A_rows + j] = sum;
        //     // for (int k {}; k < row_sizes[i]; ++k) {
        //     //     Y[i * A_rows + j] +=
        //     //       data[row_locations[i] + k]
        //     //       * A[j * A_cols + global_col_indices[i] + k];
        //     // }
        //     // Y[i * A_rows + j] = ddot_(&row_sizes[i],
        //     //                           &data[row_locations[i]],
        //     //                           &i_one,
        //     //                           &A[j * A_cols + global_col_indices[i]],
        //     //                           &i_one);
        // }
    }
}

auto PackedMatrix::multRightTransposed(const PackedMatrix& A,
                                       std::vector<double>& Y) const -> void
{
    for (int i {}; i < n_rows; ++i) {
        for (int j {}; j < A.n_rows; ++j) {
            const int i_beg { global_col_indices[i] };
            const int i_end { i_beg + row_sizes[i] };
            const int j_beg { A.global_col_indices[j] };
            const int j_end { j_beg + A.row_sizes[j] };
            if (i_beg <= j_beg) {
                const int i_start_shift { j_beg - i_beg };
                const int ni { row_sizes[i] - i_start_shift };
                const int nj { A.row_sizes[j] };
                const int n_el { std::min(ni, nj) };
                for (int k {}; k < n_el; ++k) {
                    const int idx_i { row_locations[i] + i_start_shift + k };
                    const int idx_j { A.row_locations[j] + k };
                    Y[i * A.n_rows + j] += data[idx_i] * A.data[idx_j];
                }
            }
            if (i_beg > j_beg) {
                const int j_start_shift { i_beg - j_beg };
                const int nj { A.row_sizes[j] - j_start_shift };
                const int ni { row_sizes[i] };
                const int n_el { std::min(ni, nj) };
                for (int k {}; k < n_el; ++k) {
                    const int idx_j { A.row_locations[j] + j_start_shift + k };
                    const int idx_i { row_locations[i] + k };
                    Y[i * A.n_rows + j] += data[idx_i] * A.data[idx_j];
                }
            }
        }
    }
}

auto PackedMatrix::multLeftTransposed(const std::vector<double>& A,
                                      std::vector<double>& Y) const -> void
{
    const int A_cols { n_cols };
    const int A_rows { static_cast<int>(A.size()) / A_cols };
    assert(static_cast<int>(Y.size()) >= A_rows * n_rows);
    for (int i {}; i < A_rows; ++i) {
        for (int j {}; j < n_rows; ++j) {
            Y[i * n_rows + j] = ddot_(&row_sizes[j],
                                      &A[i * n_cols + global_col_indices[j]],
                                      &i_one,
                                      &data[row_locations[j]],
                                      &i_one);
        }
    }
}

auto PackedMatrix::multLeftTransposedT(const std::vector<double>& A,
                                       std::vector<double>& Y) const -> void
{
    const int A_cols { n_cols };
    const int A_rows { static_cast<int>(A.size()) / A_cols };
    assert(static_cast<int>(Y.size()) >= A_rows * n_rows);
    for (int j {}; j < n_rows; ++j) {
        for (int i {}; i < A_rows; ++i) {
            Y[j * A_rows + i] = ddot_(&row_sizes[j],
                                      &A[i * n_cols + global_col_indices[j]],
                                      &i_one,
                                      &data[row_locations[j]],
                                      &i_one);
        }
    }
}

auto PackedMatrix::transpose(const PackedMatrix& X) -> void
{
    data.resize(X.data.size());
    int data_idx {};
    for (int i_row {}; i_row < n_rows; ++i_row) {
        row_locations[i_row] = data_idx;
        for (int m {}; m < n_cols; ++m) {
            const int col_idx_first { X.global_col_indices[m] };
            const int col_idx_last { col_idx_first + X.row_sizes[m] };
            if (i_row >= col_idx_first && i_row < col_idx_last) {
                if (global_col_indices[i_row] == -1) {
                    global_col_indices[i_row] = m;
                }
                ++row_sizes[i_row];
                data[data_idx++] = X.data[X.row_locations[m] + i_row - col_idx_first];
            }
        }
    }
}

auto PackedMatrix::fullMatrix(std::vector<double>& full_matrix) const -> void
{
    full_matrix.resize(n_rows * n_cols, 0.0);
    for (int i {}; i < n_rows; ++i) {
        for (int k {}; k < row_sizes[i]; ++k) {
            full_matrix[i * n_cols + global_col_indices[i] + k] =
              data[row_locations[i] + k];
        }
    }
}
