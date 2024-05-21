// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "array.h"

string now_timestring( // {{{
)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int millisec = lrint(tv.tv_usec/1000.0); // Round to nearest millisecond.
    // If rounding up to nearest millisecond rounds up to the next second.
    if (millisec >= 1000) {
        millisec -= 1000;
        tv.tv_sec++;
    }
    struct tm *tm_info = localtime(&tv.tv_sec);
    const size_t sz = 20; // With \0.
    char buffer[sz];
    strftime(buffer,sz,"%FT%T",tm_info);
    return format("%s.%03d",buffer,millisec);
} // }}}

// Julian date. That is the day number from January 1st 4713 BC. That
// is the start of the Julian calander. Note that the year 4713 BC should
// be referred to as year -4712, because there is no year 0. Just in case
// we ever want to launch a satellite before the year 1.
int juliandate( // {{{
    int day,
    int month,
    int year
)
{
    // This code comes from the internet.
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;
    if (year > 1582 || (year == 1582 && month > 10) || (year == 1582 && month == 10 && day >= 15)) return day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
    else return day + (153 * m + 2) / 5 + 365 * y + y / 4 - 32083;
} // }}}

// For some reason, the standard library <format> is not available, so we copy and paste the
// implementation from the internet.
string format( // {{{
    const string& fmt,
    ...
)
{
    va_list args;
    va_start (args, fmt);
    size_t len = vsnprintf(NULL, 0, fmt.c_str(), args);
    va_end (args);
    vector<char> vec(len + 1);
    va_start (args, fmt);
    vsnprintf(&vec[0], len + 1, fmt.c_str(), args);
    va_end (args);
    return &vec[0];
} // }}}

// Solve a matrix equation using Cholevsky decomposition.
int chol_solve( // {{{
    size_t nstate, // Size of the matrix (square).
    double *mat // Matrix (will be turned into her inverse. It is represented as just one triangle, because the matrix must be symmetric for this algorithm.
)
{

    // This routine originated from Fortran, but we now represent the matrix as
    // just the independent triangle. For instance, if the dimensions is 5, the
    // array mat only has 15 elements, not 25.
    double *row_state = mat;
    for (size_t istate=0 ; istate < nstate ; istate++) {
        row_state += istate; // Set the pointer to the beginning of your row.
        double c = row_state[istate];
        for (size_t ilow=0 ; ilow < istate ; ilow++) {
            c -= pow(row_state[ilow],2.0);
        }
        if (c < 0.0) return 1;
        row_state[istate] = sqrt(c);

        double *row_high = row_state;
        for (size_t ihigh=istate+1 ; ihigh<nstate ; ihigh++) {
            row_high += ihigh; // Beginning of row ihigh.
            c = row_high[istate];
            for (size_t ilow=0 ; ilow < istate ; ilow++) {
                c -= row_high[ilow]*row_state[ilow];
            }
            row_high[istate] = c/row_state[istate];
        }
    }

    // Promote the matrix to the inverse or the original.
    double *row_low = mat;
    for (size_t ilow=0 ; ilow < nstate ; ilow++) {
        row_low += ilow;
        row_low[ilow] = 1.0/row_low[ilow];
        row_state = row_low;
        for (size_t istate=ilow+1 ; istate < nstate ; istate++) {
            row_state += istate;
            row_state[ilow] *= -row_low[ilow];
        }
        row_state = row_low;
        for (size_t istate=ilow+1 ; istate < nstate ; istate++) {
            row_state += istate;
            row_state[ilow] /= row_state[istate];
            double *row_high = row_state;
            for (size_t ihigh=istate+1 ; ihigh < nstate ; ihigh++) {
                row_high += ihigh;
                row_high[ilow] -= row_high[istate]*row_state[ilow];
            }
        }
    }

    row_low = mat;
    for (size_t ilow=0 ; ilow < nstate ; ilow++) {
        row_low += ilow;
        row_state = row_low;
        for (size_t istate=ilow ; istate < nstate ; istate++) {
            double c = 0.0;
            double *row_high = row_state;
            for (size_t ihigh=istate ; ihigh < nstate ; ihigh++) {
                c += row_high[ilow]*row_high[istate];
                row_high += ihigh+1;
            }
            row_state[ilow] = c;
            row_state += istate+1;
        }
    }

    return 0;

} // }}}

void median( // {{{
    size_t sz,
    const double *arr,
    double &res
)
{
    vector<double> draft(sz);
    for (size_t i=0 ; i<sz ; i++) draft[i] = arr[i];
    Array<double>::select(sz,(sz-1)/2,sz/2,draft.data(),res);
} // }}}
void fancymedian( // {{{
    size_t sz,
    const double *arr,
    double &res
)
{
    median(sz,arr,res);
    vector<double> del(sz);
    vector<double> absdel(sz);
    for (size_t i=0 ; i<sz ; i++) {
        del[i] = arr[i]-res;
        absdel[i] = abs(del[i]);
    }
    double medabsdel;
    median(sz,absdel.data(),medabsdel);
    if (medabsdel == 0) return;
    vector<double> w(sz);
    double wmag = 0;
    for (size_t i=0 ; i<sz ; i++) {
        w[i] = pow(max(0.0,1.0-pow(del[i]/(6.0*medabsdel),2.0)),2.0);
        wmag += w[i];
    }
    for (size_t i=0 ; i<sz ; i++) {
        res += w[i]*del[i]/wmag;
    }
} // }}}

// Linear interpolation (not using B-splines, because it is no rocket science).
// The independent values can be ascending or descending per array.
void linear_interpol( // {{{
    size_t sz1, // Size of the original array.
    size_t sz2, // Size of the array of interpolated values.
    double *z1, // Independent variable of original array.
    double *z2, // Independent variable of target array.
    double *f1, // Dependent variable of original array.
    double *f2 // Dependent variable of target array (output).
)
{

    // This routine performs linear extrapolation if the output abscissa is out of range.
    // You cannot interpolate or extrapolate if the input is only one large.
    if (sz1 == 1) {
        // Output is equal to input independent of abscissa.
        for (size_t k2=0 ; k2<sz2 ; k2++) {
            f2[k2] = f1[0];
        }
        return;
    }

    // This function is literally translated from Fortran.

    // z1,f1 and z2,f2 may go from top to bottom, x1,y1 and x2,y2
    // go from bottom to top.

    // z1,f1 and z2,f2 must be either ascending or descending, not
    // a combination of both.

    // Defines ascending arrays x1 and y1. These are either copies of
    // z1 and f1, or their reverses, depending on the order of z1.
    vector<double> x1(sz1);
    vector<double> y1(sz1);
    if (z1[0] > z1[sz1-1]) {
        for (size_t k1=0 ; k1<sz1 ; k1++) {
            x1[k1] = z1[sz1-1-k1];
            y1[k1] = f1[sz1-1-k1];
        }
    } else {
        for (size_t k1=0 ; k1<sz1 ; k1++) {
            x1[k1] = z1[k1];
            y1[k1] = f1[k1];
        }
    }

    vector<double> y2(sz2);
    double a;
    double b;

    // The same for x2 with respect of z2, only there is no y2 or f2 yet.
    // The array x2 will be ascending. Later on, the interpolated array
    // may be reversed, depending on whether x2 is a reverse of z2 or not.
    vector<double> x2(sz2);
    if (z2[0] > z2[sz2-1]) {
        for (size_t k2=0 ; k2<sz2 ; k2++) {
            x2[k2] = z2[sz2-1-k2];
        }
    } else {
        for (size_t k2=0 ; k2<sz2 ; k2++) {
            x2[k2] = z2[k2];
        }
    }

    // Interpolation indices. Initialize at lower index at 1, which means
    // that the value should be between 1 and 2. These values may increase
    // later on. Even if the value on which is interpolated is outside
    // domain x1, it will be interpreted as if it is between the last two
    // indices at that side.
    size_t m1 = 0;

    for (size_t k2=0 ; k2<sz2 ; k2++) {
        // Update interpolation index m2 if necessary, but not further than
        // sz1. As x2 is increasing, the interpolation indices are expected
        // to increase steadily. At least, they never need to decrease.
        size_t m2;
        for (m2=m1+1 ; m2<sz1 ; m2++) { // So that sz1-1 is the end result if x1 is very large.
            if (x2[k2] <= x1[m2] || m2 == sz1-1) break; // The iterator lives outside the loop.
        }
        // Lower index is always one lower than the upper index.
        m1 = m2 - 1;

        if (y1[m1] == NC_FILL_DOUBLE || y1[m2] == NC_FILL_DOUBLE) {
            y2[k2] = NC_FILL_DOUBLE;
        } else {

            // Linear relationship.
            a = (y1[m1]-y1[m2])/(x1[m1]-x1[m2]);
            b = y1[m1]-a*x1[m1];

            y2[k2] = a*x2[k2] + b;

        }
    }

    // Reverse the result if the target domain was reversed at the start.
    if (z2[0] > z2[sz2-1]) {
        for (size_t k2=0 ; k2<sz2 ; k2++) {
            f2[k2] = y2[sz2-1-k2];
        }
    } else {
        for (size_t k2=0 ; k2<sz2 ; k2++) {
            f2[k2] = y2[k2];
        }
    }
} // }}}

void fill_holes( // {{{
    size_t sz,
    vector<bool> &mask,
    double *arr
)
{
    for (size_t i=0 ; i<sz ; i++) {
        if (!mask[i]) continue;
        for (size_t movement=1 ; ; movement++) {
            bool left;
            bool right;
            if (movement > i) left = false;
            else left = !mask[i-movement];
            if (movement > sz-1-i) right = false;
            else right = !mask[i+movement];
            if (left) {
                if (right) arr[i] = 0.5 * (arr[i-movement] + arr[i+movement]);
                else arr[i] = arr[i-movement];
            } else {
                if (right) arr[i] = arr[i+movement];
                else continue;
            }
            break;
        }
    }
} // }}}

