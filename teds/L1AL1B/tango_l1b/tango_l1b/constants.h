// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

namespace tango {

// Dimensions of mathematical spaces
namespace dims {

// Number of vector components
constexpr int vec { 3 };

} // namespace dims

// Fill values to denote a missing or default value
namespace fill {

constexpr int i { -32767 };
constexpr float f { static_cast<float>(i) };
constexpr double d { static_cast<double>(i) };

} // namespace fill

// Process ladder, an ordered list of possible states of
// data. Possible data levels run from L1A to L1B.
enum class ProcLevel
{
    l1a,
    raw, // Like L1A but data stored as floating point numbers
    dark,
    noise,
    nonlin,
    prnu,
    stray,
    swath,
    wave,
    rad,
    l1b,
    n_levels,
};

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
