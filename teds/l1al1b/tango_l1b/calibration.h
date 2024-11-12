// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// List of calibration steps applied to data by the L1A-L1B
// processor. Each function should set the new data level accordingly.

#pragma once

#include "constants.h"

namespace tango {

class BinningTable;
class CKD;
class L1;

// If the data level is L1A then the detector images are yet to be
// divided by the bin sizes of the binning table;
auto binningTable(const CKD& ckd,
                  const BinningTable& binning_table,
                  L1& l1) -> void;

// Remove dark offset
auto darkOffset(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Get an initial noise estimate
auto noise(const CKD& ckd,
           const bool enabled,
           const BinningTable& binning_table,
           L1& l1) -> void;

// Remove dark current. This is split from dark offset because noise
// calibration needs to happen in betwee.
auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Correct for nonlinearity
auto nonlinearity(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Correct for photoresponse non-uniformity and quantum efficiency
auto prnu(const CKD& ckd, const bool enabled, L1& l1) -> void;

// Smooth over bad values. This is necessary for algorithms such as
// stray light correction which use all pixels.
auto removeBadValues(const CKD& ckd, L1& l1) -> void;

// Correct for stray light
auto strayLight(const CKD& ckd,
                const BinningTable& binning_table,
                const int n_van_cittert,
                L1& l1) -> void;

// Extract a set of spectra from one detector image
auto mapFromDetector(const CKD& ckd,
                     const BinningTable& binning_table,
                     const int b_spline_order,
                     L1& l1) -> void;

// Mapping from the detector yields spectra on the intermediate
// wavelengths grid. Spectra need to be interpolated onto the CKD
// wavelength grids. This function is separate from mapFromDetector
// because the latter might not always run depending on the input data
// level.
auto changeWavelengthGrid(const CKD& ckd, L1& l1) -> void;

// Radiometrically calibrate
auto radiometric(const CKD& ckd, const bool enabled, L1& l1) -> void;

} // namespace tango
