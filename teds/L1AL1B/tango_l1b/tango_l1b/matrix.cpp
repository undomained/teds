// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "matrix.h"

void Matrix::matmul_sl_sym( // {{{
    size_t sz,
    squareoption_t opt,
    double *sq,
    double *vec,
    double *res // Also sz.
)
{
    if (opt == OPT_FULL) {
        double *row = sq;
        for (size_t ihigh=0 ; ihigh<sz ; ihigh++) {
            row += ihigh;
            res[ihigh] = 0.0;
            for (size_t ilow=0 ; ilow<=ihigh ; ilow++) {
                res[ihigh] += row[ilow] * vec[ilow];
                if (ilow != ihigh) res[ilow] += row[ilow] * vec[ihigh];
            }
        }
    } else { // OPT_DIAG, otherwise, do not call this routine.
        for (size_t i=0 ; i<sz ; i++) res[i] = sq[i] * vec[i];
    }
} // }}}
void Matrix::matmul_rl_fold_quick( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    double *vec,
    double *res // Slow.
)
{
    double *mat_row = mat;
    for (size_t islow=0 ; islow<nslow ; islow++) {
        res[islow] = 0.0;
        for (size_t iquick=0 ; iquick<nquick ; iquick++) res[islow] += mat_row[iquick] * vec[iquick];
        mat_row += nquick;
    }
} // }}}
void Matrix::matmul_rl_fold_slow( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    double *vec,
    double *res // Quick.
)
{
    for (size_t iquick=0 ; iquick<nquick ; iquick++) res[iquick] = 0.0;
    double *mat_row = mat;
    for (size_t islow=0 ; islow<nslow ; islow++) {
        for (size_t iquick=0 ; iquick<nquick ; iquick++) res[iquick] += mat_row[iquick] * vec[islow];
        mat_row += nquick;
    }
} // }}}

void Matrix::matmul_rs_fold_quick_sym( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    squareoption_t opt_quick, // Only DIAG or FULL. For NONE, this routine should not be called.
    double *sq_quick,
    double *res
)
{
    double *row_res = res;
    double *row_mat = mat;
    for (size_t islow=0 ; islow<nslow ; islow++) {
        if (opt_quick == OPT_FULL) {
            double *row_quick = sq_quick;
            for (size_t iquick_high=0 ; iquick_high<nquick ; iquick_high++) {
                row_quick += iquick_high;
                row_res[iquick_high] = 0.0;
                for (size_t iquick_low=0 ; iquick_low<=iquick_high ; iquick_low++) {
                    row_res[iquick_high] += row_quick[iquick_low] * row_mat[iquick_low];
                    if (iquick_low != iquick_high) row_res[iquick_low] += row_quick[iquick_low] * row_mat[iquick_high];
                }
            }
        } else { // opt_quick is OPT_DIAG.
            for (size_t iquick=0 ; iquick<nquick ; iquick++) row_res[iquick] = row_mat[iquick] * sq_quick[iquick];
        }
        row_mat += nquick;
        row_res += nquick;
    }
} // }}}
void Matrix::matmul_rs_fold_slow_sym( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    squareoption_t opt_slow, // Only DIAG or FULL. For NONE, this routine should not be called.
    double *sq_slow,
    double *res
)
{
    if (opt_slow == OPT_FULL) {
        double *row_slow = sq_slow;
        double *row_res_high = res;
        double *row_mat_high = mat;
        for (size_t islow_high=0 ; islow_high<nslow ; islow_high++) {
            row_slow += islow_high;
            for (size_t iquick=0 ; iquick<nquick ; iquick++) row_res_high[iquick] = 0.0;
            double *row_res_low = res;
            double *row_mat_low = mat;
            for (size_t islow_low=0 ; islow_low<=islow_high ; islow_low++) {
                for (size_t iquick=0 ; iquick<nquick ; iquick++) {
                    row_res_high[iquick] += row_slow[islow_low] * row_mat_low[iquick];
                    if (islow_low != islow_high) row_res_low[iquick] += row_slow[islow_low] * row_mat_high[iquick];
                }
                row_res_low += nquick;
                row_mat_low += nquick;
            }
            row_res_high += nquick;
            row_mat_high += nquick;
        }
    } else { // opt_slow is OPT_DIAG.
        double *row_res = res;
        double *row_mat = mat;
        for (size_t islow=0 ; islow<nslow ; islow++) {
            for (size_t iquick=0 ; iquick<nquick ; iquick++) row_res[iquick] = row_mat[iquick] * sq_slow[islow];
            row_mat += nquick;
            row_res += nquick;
        }
    }
} // }}}
void Matrix::matmul_rr_fold_quick_sym( // {{{
    size_t nslow,
    size_t nquick,
    double *mat1,
    double *mat2,
    double *res // Slow by slow, symmetric triangle definition.
)
{
    double *row_res = res;
    double *row1 = mat1;
    for (size_t islow_high=0 ; islow_high<nslow ; islow_high++) {
        row_res += islow_high;
        double *row2 = mat2;
        for (size_t islow_low=0 ; islow_low<nslow ; islow_low++) {
            row_res[islow_low] = 0.0;
            for (size_t iquick=0 ; iquick<nquick ; iquick++) row_res[islow_low] += row1[iquick] * row2[iquick];
            row2 += nquick;
        }
        row1 += nquick;
    }
} // }}}
void Matrix::matmul_rr_fold_slow_sym( // {{{
    size_t nslow,
    size_t nquick,
    double *mat1,
    double *mat2,
    double *res // Quick by quick, symmetric triangle definition.
)
{
    double *row1 = mat1;
    double *row2 = mat2;
    for (size_t ires=0 ; ires<nquick*(nquick+1)/2 ; ires++) res[ires] = 0.0;
    for (size_t islow=0 ; islow<nslow ; islow++) {
        double *row_res = res;
        for (size_t iquick_high=0 ; iquick_high<nquick ; iquick_high++) {
            row_res += iquick_high;
            for (size_t iquick_low=0 ; iquick_low<=iquick_high ; iquick_low++) {
                row_res[iquick_low] += row1[iquick_high] * row2[iquick_low];
            }
        }
        row1 += nquick;
        row2 += nquick;
    }
} // }}}

void Matrix::matmul_rr_fold_quick_asym( // {{{
    size_t nslow_slowside,
    size_t nslow_quickside,
    size_t nquick,
    double *mat_slowside, // Surviving dimension will be the slow one in res.
    double *mat_quickside, // Surviving dimension will be the quick one in res.
    double *res // Slow by slow, the slow side of mat_quickside is here quick.
)
{
    double *row_res = res;
    double *row_slowside = mat_slowside;
    for (size_t islow_slowside=0 ; islow_slowside<nslow_slowside ; islow_slowside++) {
        double *row_quickside = mat_quickside;
        for (size_t islow_quickside=0 ; islow_quickside<nslow_quickside ; islow_quickside++) {
            row_res[islow_quickside] = 0.0;
            for (size_t iquick=0 ; iquick<nquick ; iquick++) row_res[islow_quickside] += row_slowside[iquick] * row_quickside[iquick];
            row_quickside += nquick;
        }
        row_res += nslow_quickside;
        row_slowside += nquick;
    }
} // }}}
void Matrix::matmul_rr_fold_slow_asym( // {{{
    size_t nslow,
    size_t nquick_slowside,
    size_t nquick_quickside,
    double *mat_slowside, // Surviving dimension will be the slow one in res.
    double *mat_quickside, // Surviving dimension will be the quick one in res.
    double *res // Quick by quick, the quick side of mat_slowside is here slow.
)
{
    for (size_t ires=0 ; ires<nquick_slowside*nquick_quickside ; ires++) res[ires] = 0.0;
    double *row_slowside = mat_slowside;
    double *row_quickside = mat_quickside;
    for (size_t islow=0 ; islow<nslow ; islow++) {
        double *row_res = res;
        for (size_t iquick_slowside=0 ; iquick_slowside<nquick_slowside ; iquick_slowside++) {
            for (size_t iquick_quickside=0 ; iquick_quickside<nquick_quickside ; iquick_quickside++) row_res[iquick_quickside] += row_slowside[iquick_slowside] * row_quickside[iquick_quickside];
            row_res += nquick_quickside;
        }
        row_slowside += nquick_slowside;
        row_quickside += nquick_quickside;
    }
} // }}}

void Matrix::matmul_tsk( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    squareoption_t opt_quick,
    double *sq_quick,
    double *res // Slow by slow.
)
{
    // First calculate s*k.
    vector<double> sk_new;
    double *sk;
    if (opt_quick == OPT_NONE) {
        sk = mat;
    } else {
        sk_new.resize(nslow*nquick);
        sk = sk_new.data();
        matmul_rs_fold_quick_sym(nslow,nquick,mat,opt_quick,sq_quick,sk);
    }
    matmul_rr_fold_quick_sym(nslow,nquick,mat,sk,res);
} // }}}
void Matrix::matmul_sts( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    squareoption_t opt_slow,
    double *sq_slow,
    squareoption_t opt_quick,
    double *sq_quick,
    double *res
)
{
    vector<double> ts_new;
    double *ts;
    if (opt_quick == OPT_NONE) {
        if (opt_slow == OPT_NONE) memcpy(res,mat,nslow*nquick*sizeof(double)); // res = mat.
        else ts = mat;
    } else {
        if (opt_slow == OPT_NONE) ts = res;
        else {
            ts_new.resize(nslow*nquick);
            ts = ts_new.data();
        }
        matmul_rs_fold_quick_sym(nslow,nquick,mat,opt_quick,sq_quick,ts);
    }
    if (opt_slow != OPT_NONE) matmul_rs_fold_slow_sym(nslow,nquick,ts,opt_slow,sq_slow,res);
} // }}}
void Matrix::matmul_kst( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    squareoption_t opt_slow,
    double *sq_slow,
    double *res // Quick by quick.
)
{
    // First calculate s*t.
    vector<double> st_new;
    double *st;
    if (opt_slow == OPT_NONE) {
        st = mat;
    } else {
        st_new.resize(nslow*nquick);
        st = st_new.data();
        matmul_rs_fold_slow_sym(nslow,nquick,mat,opt_slow,sq_slow,st);
    }
    matmul_rr_fold_slow_sym(nslow,nquick,mat,st,res);
} // }}}
// Matrix multiplication skipping intermediate products. For instance,
// to do one calculation and not want to put the entire gain matrix on the
// memory.
// From inverted matrix to final retrieval.
void Matrix::matmul_stsy( // {{{
    size_t nslow,
    size_t nquick,
    double *mat,
    squareoption_t opt_slow,
    double *sq_slow,
    squareoption_t opt_quick,
    double *sq_quick,
    double *vec,
    double *res // Slow.
)
{
    // Calculate s*y.
    vector<double> sy_new;
    double *sy;
    if (opt_quick == OPT_NONE) {
        sy = vec;
    } else {
        sy_new.resize(nquick);
        sy = sy_new.data();
        matmul_sl_sym(nquick,opt_quick,sq_quick,vec,sy);
    }
    // Calculate t*sy.
    vector<double> tsy_new;
    double *tsy;
    if (opt_slow == OPT_NONE) {
        tsy = res;
    } else {
        tsy_new.resize(nslow);
        tsy = tsy_new.data();
    }
    matmul_rl_fold_quick(nslow,nquick,mat,sy,tsy);
    // Calculate s*tsy
    if (opt_slow != OPT_NONE) matmul_sl_sym(nslow,opt_slow,sq_slow,tsy,res);
} // }}}

void Matrix::printmatrix_sym( // {{{
    Logger *writer,
    size_t sz,
    double *mat
)
{
    vector<double> expanded(sz*sz);
    double *row = mat;
    for (size_t ihigh=0 ; ihigh<sz ; ihigh++) {
        row += ihigh;
        for (size_t ilow=0 ; ilow<=ihigh ; ilow++) {
            expanded[ihigh*sz+ilow] = row[ilow];
            if (ilow != ihigh) expanded[ilow*sz+ihigh] = row[ilow];
        }
    }
    double *row_exp = expanded.data();
    for (size_t islow=0 ; islow<sz ; islow++) {
        string line = "";
        for (size_t iquick=0 ; iquick<sz ; iquick++) line += format("%16.6e",row_exp[iquick]);
        writer->writelog(log_debug,"%s",line.c_str());
        row_exp += sz;
    }
} // }}}

