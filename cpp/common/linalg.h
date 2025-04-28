// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Linear algebra operations

#pragma once

#include "eigen.h"

#include <Eigen/Sparse>

namespace tango {

// Convert a sparse matrix to a banded matrix. The matrix is assumed
// to symmetric and only the lower half is stored.
auto sparseToBanded(const Eigen::SparseMatrix<double>& sparse)
  -> Eigen::MatrixXd;

// LDLT decomposition
auto ldlt_decompose(Eigen::MatrixXd& L, Eigen::VectorXd& D) -> void;

// Solve linear system Ax=b for x where A=LDLT has been computed using
// the previous function. Returns x.
auto ldlt_solve(const Eigen::MatrixXd& L,
                const Eigen::VectorXd& D,
                const Eigen::MatrixXd& B) -> Eigen::MatrixXd;

} // namespace tango
