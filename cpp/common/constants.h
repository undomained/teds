// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

namespace tango {

// Dimensions of mathematical spaces
namespace dims {

// Number of vector components
constexpr int vec { 3 };
// Number of quaternion components
constexpr int quat { 4 };

} // namespace dims

// Fill values to denote a missing or default value
namespace fill {

constexpr int i { -32767 };
constexpr float f { static_cast<float>(i) };
constexpr double d { static_cast<double>(i) };

} // namespace fill

namespace math {

// Multiply with this factor to convert from degrees to radians
constexpr double deg_to_rad { 0.017453292519943295 };

} // namespace math

// Geolocation related
namespace earth {

constexpr double a { 6378137.0 };
constexpr double b { 6356752.314245179499 };
constexpr double f_inv { 298.257223563 };
constexpr double f { 3.3528106647474804385e-3 };
constexpr double e2 { 6.694379990141316435e-3 };
constexpr double e { 8.1819190842621490908e-2 };

} // namespace earth

// Atmosphere related
namespace atm {

// Avogadro number [#/mol]
constexpr double na { 6.02214076e23 };
// Molar mass dry air [kg/mol]
constexpr double mdryair { 28.9647e-3 };
// Standard gravity at ground [m/s2]
constexpr double g0 { 9.80665 };
// Molecular oxygen mole fraction [-]
constexpr double xo2 { 0.2095 };
// std CH4 volume mixing ratio 1800 ppb
constexpr double xch4 { 1.8e-6 };
// std CO2 volume mixing ratio 410 ppm
constexpr double xco2 { 410e-6 };
// Planck's constant [J/s]
constexpr double hplanck { 6.62607015e-34 };
// Speed of light [m/s]
constexpr double clight { 2.99792458e8 };

} // namespace atm

// Process ladder, an ordered list of possible states of
// data. Possible data levels run from L1A to L1B to L2.
enum class ProcLevel
{
    l1a,
    raw, // Like L1A but data stored as floating point numbers
    dark_offset,
    noise,
    dark_current,
    nonlin,
    prnu,
    stray,
    swath,
    l1b, // Radiometrically calibrated
    sgm, // Unconvolved line-by-line spectra from the scene generation module
    l2,
    n_levels,
};

// Unbinning modes for detector images when calibrating. If none then
// the detector image is not unbinned (except for stray light) and
// instead the CKD is binned and applied to binned images.
enum class Unbin
{
    none,
    nearest, // Nearest neighbor
    linear,  // Linear spline
    cubic,   // Cubic spline
    n_types,
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
