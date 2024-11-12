// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// List of processes applied to data by the instrument model. Each
// process undoes or "uncalibrates" a step that is performed by the
// L1A-L1B processor, gradually bringing the data level from L1B to
// L1A or anywhere in between.

// Each function takes the CKD and zero or more settings as arguments
// and operates on a L1 object. The level of the L1 object can be L1B,
// L1A, or anything in between. It is the responsibility of each
// function to set the new level accordingly, although not every
// process has to change the level (see darkCurrent).

// Many of these function take the "enabled" argument, allowing the
// user to switch off a given process, but not all of them. For
// instance, drawOnDetector can not be turned off if it's in between
// starting and ending data level (cannot perform a detector step if
// the spectra are not drawn on the detector).

#pragma once

#include <Eigen/Sparse>
#include <vector>

namespace tango {

class BinningTable;
class CKD;
class L1;
class LinearSpline;

// Convolve each radiance spectrum with the ISRF. This significantly
// reduces data dimensions.
auto applyISRF(const CKD& ckd,
               const bool enabled,
               const Eigen::SparseMatrix<double>& isrf,
               L1& l1) -> void;

// Undo radiometric calibration
auto radiometric(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Map spectra to the detector. After this spectra will be deallocated
// and we work with the detector image.
auto mapToDetector(const CKD& ckd, const int b_spline_order, L1& l1) -> void;

// Add stray light to the image
auto strayLight(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Multiply detector image with photoresponse nonuniformity (PRNU) factors
auto prnu(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Convert from the ideal to measured signal using the nonlinearity CKD
auto nonlinearity(const CKD& ckd,
                  const bool enabled,
                  const LinearSpline& nonlin_spline,
                  L1& l1) -> void;

// Add dark current to the image
auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Add noise to the image
auto noise(const CKD& ckd, const bool enabled, const int seed, L1& l1) -> void;

// Add dark offset to the image
auto darkOffset(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Convert detector image to integer format. After this the detector
// image in floating point format is allocated and with work with the
// short integer format. No "enabled" argument necessary because this
// is the final process (set cal_level to raw or higher to skip this).
auto digitalToAnalog(const BinningTable& binning_table, L1& l1) -> void;

} // namespace tango
