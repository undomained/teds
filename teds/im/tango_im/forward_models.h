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

#include <string>
#include <vector>

namespace tango {

class BinningTable;
class CKD;
class L1;
class LinearSpline;

// Convolve each radiance spectrum with the ISRF. This significantly
// reduces data dimensions. The SGM file and starting ALT position are
// only used if the spectra are already not in memory. In that case
// they are read one by one for the convolution.
auto applyISRF(const CKD& ckd,
               const bool enabled,
               const double fwhm_gauss,
               const double shape,
               const std::string& sgm_filename,
               const int alt_beg,
               L1& l1_prod) -> void;

// Undo radiometric calibration
auto radiometric(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Map spectra to the detector. After this spectra will be deallocated
// and we work with the detector image. The exact algorithm draws each
// spectrum to the detector using the "up" and "down" pixels.
auto mapToDetector(const CKD& ckd,
                   const int b_spline_order,
                   const bool exact_drawing,
                   L1& l1_prod) -> void;

// Add stray light to the image
auto strayLight(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Multiply detector image with photoresponse nonuniformity (PRNU) factors
auto prnu(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Convert from the ideal to measured signal using the nonlinearity CKD
auto nonlinearity(const CKD& ckd,
                  const bool enabled,
                  const LinearSpline& nonlin_spline,
                  L1& l1_prod) -> void;

// Add dark current to the image
auto darkCurrent(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Add noise to the image
auto noise(const CKD& ckd,
           const bool enabled,
           const int seed,
           const double artificial_scaling,
           L1& l1_prod) -> void;

// Add dark offset to the image
auto darkOffset(const CKD& ckd, const bool enabled, L1& l1_prod) -> void;

// Bin all detector images using the binning table and binning ID from
// the configuration file. This is independent of the calibration
// level and can be applied always if the detector images
// exist.
// n_rows and n_cols are the dimensions of the unbinned detector
// image.
// scale_by_binsize determines whether to divide each signal value by
// the bin size. This is only relevant for L1A which must remained
// unscaled because the scaling will take place in the L1B processor.
auto binDetectorImages(const int n_rows,
                       const int n_cols,
                       const std::string& binning_filename,
                       const int binning_table_id,
                       const bool scale_by_binsize,
                       L1& l1_prod) -> void;

// Convert detector image to integer format. After this the detector
// image in floating point format is allocated and with work with the
// short integer format. No "enabled" argument necessary because this
// is the final process (set cal_level to raw or higher to skip this).
auto digitalToAnalog(const int nr_coadditions, L1& l1_prod) -> void;

// Compute and print the optimal coadding factor and exposure
// time. This is only informational and does not affect any
// observables.
auto estimateOptimalCoadd(const CKD& ckd,
                          const int FMC,
                          const double exposure_time,
                          const double f_sat,
                          const double full_well,
                          const double t_dead,
                          const L1& l1_prod) -> void;

} // namespace tango
