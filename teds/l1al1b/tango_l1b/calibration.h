// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// List of calibration steps applied to data by the L1A-L1B
// processor. Each function should set the new data level accordingly.

#pragma once

#include "constants.h"

#include <cstdint>

namespace tango {

class BinningTable;
class CKD;
class L1;

// If the data level is L1A then the detector images are yet to be
// divided by the coadding factors and bin sizes of the binning
// table. Also, this is where noise gets initialized.
auto binScaling(const CKD& ckd,
                const BinningTable& binning_table,
                L1& l1_prod) -> void;

// Remove dark offset
auto darkOffset(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Get an initial noise estimate
auto noise(const CKD& ckd,
           const bool enabled,
           const BinningTable& binning_table,
           const double artificial_scaling,
           L1& l1_prod) -> void;

// Remove dark current. This is split from dark offset because noise
// calibration needs to happen in betwee.
auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Correct for nonlinearity
auto nonlinearity(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Correct for photoresponse non-uniformity and quantum efficiency
auto prnu(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Smooth over bad values using 1D cubic splines across rows. This is
// necessary for algorithms such as stray light correction which use
// all pixels.
auto removeBadValues(const CKD& ckd, L1& l1_prod) -> void;

// Correct for stray light
auto strayLight(const CKD& ckd,
                const BinningTable& binning_table,
                const int n_van_cittert,
                L1& l1_prod) -> void;

// Extract a set of spectra from one detector image. For exact_drawing
// option, see the mapping algorithm comments in IM.
auto mapFromDetector(const CKD& ckd,
                     const BinningTable& binning_table,
                     const int b_spline_order,
                     const bool exact_drawing,
                     L1& l1_prod) -> void;

// Mapping from the detector yields spectra on the intermediate
// wavelengths grid. Spectra need to be interpolated onto the CKD
// wavelength grids. This function is separate from mapFromDetector
// because the latter might not always run depending on the input data
// level.
auto changeWavelengthGrid(const CKD& ckd, L1& l1_prod) -> void;

// Radiometrically calibrate
auto radiometric(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

auto binL1B(const int bin, L1& l1_prod) -> void;

} // namespace tango
