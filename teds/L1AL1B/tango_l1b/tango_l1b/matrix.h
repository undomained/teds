// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef MATRIX_H
#define MATRIX_H

#include "header.h"

// Forward declaration.
class Logger;

enum squareoption_t {
    OPT_NONE = 0,
    OPT_DIAG = 1,
    OPT_FULL = 2
};

// Matrix multiplication module.
// This is a class with just public static routines.
// No instances should ever be created (constructor is private).
class Matrix {

    private:
    Matrix() {} // Private constructor.

    public:

    // Regular matrix multiplications.
    // l stands for a one-dimensional array.
    // s stands for a square (triangle-defined) symmetric matrix.
    // r stands for a rectangular or square matrix without symmetry.
    static void matmul_sl_sym(
        size_t sz,
        squareoption_t opt,
        double *sq,
        double *vec,
        double *res // Also sz.
    );
    static void matmul_rl_fold_quick(
        size_t nslow,
        size_t nquick,
        double *mat,
        double *vec,
        double *res // Slow.
    );
    static void matmul_rl_fold_slow(
        size_t nslow,
        size_t nquick,
        double *mat,
        double *vec,
        double *res // Quick.
    );
    static void matmul_rs_fold_quick_sym(
        size_t nslow,
        size_t nquick,
        double *mat,
        squareoption_t opt_quick, // Only DIAG or FULL. For NONE, this routine should not be called.
        double *sq_quick,
        double *res
    );
    static void matmul_rs_fold_slow_sym(
        size_t nslow,
        size_t nquick,
        double *mat,
        squareoption_t opt_slow, // Only DIAG or FULL. For NONE, this routine should not be called.
        double *sq_slow,
        double *res
    );
    static void matmul_rr_fold_quick_sym(
        size_t nslow,
        size_t nquick,
        double *mat1,
        double *mat2,
        double *res // Slow by slow, symmetric triangle definition.
    );
    static void matmul_rr_fold_slow_sym(
        size_t nslow,
        size_t nquick,
        double *mat1,
        double *mat2,
        double *res // Quick by quick, symmetric triangle definition.
    );
    static void matmul_rr_fold_quick_asym(
        size_t nslow_slowside,
        size_t nslow_quickside,
        size_t nquick,
        double *mat_slowside, // Surviving dimension will be the slow one in res.
        double *mat_quickside, // Surviving dimension will be the quick one in res.
        double *res // Slow by slow, the slow side of mat_quickside is here quick.
    );
    static void matmul_rr_fold_slow_asym(
        size_t nslow,
        size_t nquick_slowside,
        size_t nquick_quickside,
        double *mat_slowside, // Surviving dimension will be the slow one in res.
        double *mat_quickside, // Surviving dimension will be the quick one in res.
        double *res // Quick by quick, the quick side of mat_slowside is here slow.
    );

    // Triple or quadruple matrix multiplications.
    // Here, we use metaphors from linear inversion, because that is the most
    // common practise.
    // k is a matrix like a Jacobian.
    // The quick dimensions folds left and the slow dimension folds right.
    // t is a matrix like the transpose of a Jacobian, defined just as the Jacobian.
    // The quick dimensions folds right and the slow dimension folds left.
    // Any time k and t are in the same operator, they are assumed to be the same
    // matrix besides transposition.
    // s is a matrix like an error covariance matrix, a symmetric matrix defined
    // in the triangular way.
    // g is a matrix like the gain matrix, the same shape as t, but this is a
    // different matrix.
    // Triple matrix multiplications.
    // k is a matrix like a Jacobian.
    // The quick dimensions folds left and the slow dimension folds right.
    // t is a matrix like the transpose of a Jacobian, defined just as the Jacobian.
    // The quick dimensions folds right and the slow dimension folds left.
    // Any time k and t are in the same operator, they are assumed to be the same
    // matrix besides transposition.
    // s is a matrix like an error covariance matrix, a symmetric matrix defined
    // in the triangular way.
    // g is a matrix like the gain matrix, the same shape as t, but this is a
    // different matrix.
    // y is a vector, like a measurement vector, it is one-dimensional.

    // tsk is K^T * S_y * K
    // Used for well-posed matrix-to-be-inverted.
    static void matmul_tsk(
        size_t nslow,
        size_t nquick,
        double *mat,
        squareoption_t opt_quick,
        double *sq_quick, // Matrix S, Diagonals D or nothing (NULL).
        double *res
    );

    // sts = S_x * K^T * S_y
    // Used for well-posed if slow S is inverted matrix and quick S is error matrix.
    // Used for ill-posed if slow S is error matrix and quick S is inverted matrix.
    static void matmul_sts(
        size_t nslow,
        size_t nquick,
        double *mat,
        squareoption_t opt_slow,
        double *sq_slow,
        squareoption_t opt_quick,
        double *sq_quick,
        double *res
    );

    // kst = K * S_x * K^T
    // Used for ill-posed matrix-to-be-inverted.
    static void matmul_kst(
        size_t nslow,
        size_t nquick,
        double *mat,
        squareoption_t opt_slow,
        double *sq_slow,
        double *res
    );
    // stsy = S_x * K^T * S_y * y
    // Take S_y^{-1} as square y matrix and inverted matrix as square x matrix
    // and you have a retrieval without saving intermediate products besides
    // the inverted matrix.
    static void matmul_stsy(
        size_t nslow,
        size_t nquick,
        double *mat,
        squareoption_t opt_slow,
        double *sq_slow,
        squareoption_t opt_quick,
        double *sq_quick,
        double *vec_quick,
        double *res
    );

    static void printmatrix_sym(
        Logger *writer,
        size_t sz,
        double *mat
    );
    static void matmul_gk(
        size_t nslow,
        size_t nquick,
        double *mat,
        double *gain, // Reverse slow and quick. They are from the Jacobian perspective.
        double *res
    );

};

#endif
