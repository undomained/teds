// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

// This h-file is always included.

// Include everything of which we want to assume that it is included.
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <memory>
#include <cstdarg>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <bitset>
#include <string.h>
#include <vector>
#include <map>
#include <math.h>
#include <stdio.h>
#include <chrono>
#include <netcdf>

namespace tango {

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

// Error handling macros. {{{

// These use macros instead of regular routines, because
// in error handling, we want to invoke a return instruction once in a while.

// This macro is only called in other macros in this module. The act argument is
// a string that says something like "Error raised" or "Error progressed".
// It writes where in which function in which file an error is raised or progressed.
#define error_message(act) \
    writelog(log_error,"%s on line %d in routine %s in file %s.",act,__LINE__,__func__,__FILE__);

// Raises an error. The arguments work the same as in a formatted print instruction.
#define raise_error(fmt,...) \
    writelog(log_error,fmt,##__VA_ARGS__); \
    error_message("Error raised"); \
    return 1;

// Conditional error raise instruction. The main reason this is a separate macro is
// that a common mistake would be to say "if (condition) raise_error(message);",
// without curly brackets. This is wrong, because the raise_error decomposes to
// multiple instruction, so curly brackets are needed. Using this macro prevents
// the user from making that mistake.
#define check_error(condition,fmt,...) \
    if (condition) { \
        raise_error(fmt,##__VA_ARGS__); \
    }

// Handles an error. This means that if the instruction goes wrong (return a nonzero),
// the calling routine also dies and reports that it progresses the error.
#define handle(funct) \
    { \
        int errorhandling_evaluation = funct; \
        if (errorhandling_evaluation != 0) { \
            error_message("Error progressed"); \
            return errorhandling_evaluation; \
        } \
    }

#define handle_nonlethal(funct) \
    { \
        int errorhandling_evaluation = funct; \
        if (errorhandling_evaluation != 0) return errorhandling_evaluation; \
    }

// }}}

// Fixed dimensions.
const size_t DIM_POL = 1;
const size_t DIM_VEC = 3;
const size_t DIM_QUAT = 4;

// One degree.
const double DEGREES = acos(0.0) / 90.0;
// Pi. That is 180 degrees. Used when 2*pi is more intuitive to say than
// 360*degrees.
const double PI = 180.0*DEGREES;

// Processing ladder.
enum level_t {
    LEVEL_DIMCAL = 0,
    LEVEL_DARKCAL = 1,
    LEVEL_NOISECAL = 2,
    LEVEL_NONLINCAL = 3,
    LEVEL_PRNUCAL = 4,
    LEVEL_STRAYCAL = 5,
    LEVEL_FOVCAL = 6,
    LEVEL_SWATHCAL = 7,
    LEVEL_WAVECAL = 8,
    LEVEL_RADCAL = 9,
    LEVEL_L1B = 11,
    LEVEL_L1C = 12,
    LEVEL_FILLVALUE = 13 // Used as non-initialized level.
};

// Number of levels.
const size_t nlevel = 13;

// L1X steps. Maybe switch to const int?
enum l1x_t {
    L1X_L1A = -1,
    L1X_RAW = 0,
    L1X_DARK = 1,
    L1X_NOISE = 2,
    L1X_NONLIN = 3,
    L1X_PRNU = 4,
    L1X_UNBIN = 5,
    L1X_STRAY = 6,
    L1X_REBIN = 7,
    L1X_FOV = 8,
    L1X_RAD = 9,
    L1X_TRUTH = 10 // Only for geo-location-only runs for simulated orbits.
};

// TODO: Razzia?
const vector<string> l1x_steps = {
    "raw",
    "dark",
    "noise",
    "nonlin",
    "prnu",
    "unbin",
    "stray",
    "rebin",
    "fov",
    "rad"
};

// Number of L1X.
const size_t nl1x = 10;

// TODO: enum type as well?
const int L1A_FLIGHT = 0;
const int L1A_CALIBRATION = 1;
const int L1A_FRAMES = 2;

// Namespace for whenever we deal with rectangles and need to identify
// their bottom, top, left, and right sides
namespace box {

constexpr int b { 0 }; // bottom
constexpr int t { 1 }; // top
constexpr int l { 2 }; // left
constexpr int r { 3 }; // right
constexpr int n { 4 }; // number of sides a box has

} // namespace box

} // namespace tango
