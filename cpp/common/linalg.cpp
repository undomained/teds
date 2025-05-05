// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "linalg.h"

namespace tango {

auto sparseToBanded(const Eigen::SparseMatrix<double>& mat_sparse)
  -> Eigen::MatrixXd
{
    // Number of non-zero elements in the first column determines the
    // number of bands including the diagonal.
    const int n_bands { mat_sparse.outerIndexPtr()[1]
                        - mat_sparse.outerIndexPtr()[0] };
    Eigen::MatrixXd mat_banded(n_bands, mat_sparse.cols());
    for (int j {}; j < mat_sparse.outerSize(); ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat_sparse, j); it;
             ++it) {
            if (it.row() >= j) {
                mat_banded(it.row() - j, j) = it.value();
            }
        }
    }
    return mat_banded;
}

auto ldlt_decompose(Eigen::MatrixXd& L, Eigen::VectorXd& D) -> void
{
    const int n_diag { static_cast<int>(L.cols()) };
    const int n_rows { static_cast<int>(L.rows()) };
    D.resize(n_diag);
    Eigen::VectorXd v(n_diag);
    for (int j {}; j < n_diag; ++j) {
        const int row_start { std::max(0, j - n_rows + 1) };
        for (int i { row_start }; i < j; ++i) {
            v(i) = L(j - i, i) * D(i);
        }
        D(j) = L(0, j);
        for (int k { row_start }; k < j; ++k) {
            D(j) -= L(j - k, k) * v(k);
        }
        for (int i { j }; i < j + n_rows; ++i) {
            const int row_start { std::max(0, i - n_rows + 1) };
            for (int k { row_start }; k < j; ++k) {
                L(i - j, j) -= L(i - k, k) * v(k);
            }
            L(i - j, j) /= D(j);
        }
    }
}

static auto banded_forward_substitution(const Eigen::MatrixXd& A,
                                        const Eigen::VectorXd& b)
  -> Eigen::VectorXd
{
    const int n_rows { static_cast<int>(A.rows()) };
    Eigen::VectorXd y(b.size());
    for (int i {}; i < A.cols(); ++i) {
        const int row_start { std::max(0, i - n_rows + 1) };
        y(i) = b(i);
        for (int j { row_start }; j < i; ++j) {
            y(i) -= A(i - j, j) * y(j);
        }
        y(i) /= A(0, i);
    }
    return y;
}

static auto banded_back_substitution(const Eigen::MatrixXd& A,
                                     const Eigen::VectorXd& b)
  -> Eigen::VectorXd
{
    const int n_diag { static_cast<int>(A.cols()) };
    const int n_rows { static_cast<int>(A.rows()) };
    Eigen::VectorXd y(b.size());
    // Backward substitution (solving A^T*x=y)
    for (int i { n_diag - 1 }; i >= 0; --i) {
        const int max_row { i + n_rows - 1 };
        y(i) = b(i);
        for (int j { std::min(n_diag - 1, max_row) }; j > i; --j) {
            // No need to explicitly use A^T here
            y(i) -= A(j - i, i) * y(j);
        }
        // In the case of LDLT algorithm the diagonal elements are
        // zero. Otherwise use y(i) /= A(0, i) here.
    }
    return y;
}

auto ldlt_solve(const Eigen::MatrixXd& L,
                const Eigen::VectorXd& D,
                const Eigen::MatrixXd& B) -> Eigen::MatrixXd
{
    Eigen::MatrixXd X(B.rows(), B.cols());
    for (int i_col {}; i_col < B.cols(); ++i_col) {
        Eigen::VectorXd b = B.col(i_col);
        // Solve Lz = b
        Eigen::VectorXd z = banded_forward_substitution(L, b);
        // Solve Dy = z
        Eigen::VectorXd y(z.array() / D.array());
        // Solve L^Tx = y
        Eigen::VectorXd x = banded_back_substitution(L, y);
        X.col(i_col) = x;
    }
    return X;
}

} // namespace tango
