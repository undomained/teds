// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "spectrum.h"

#include "ckd.h"
#include "constants.h"

namespace tango {

auto Spectrum::extract(const CKD& ckd,
                       const std::vector<double>& image,
                       const std::vector<double>& image_stdev,
                       const std::vector<bool>& pixel_mask,
                       const int i_act) -> void
{
    //signal.resize(ckd.n_detector_cols, 0.0);
    //stdev.resize(ckd.n_detector_cols, 0.0);
    //mask.resize(ckd.n_detector_cols, false);
    //for (int i {}; i < ckd.n_detector_cols; ++i) {
    //    double total_weight {};
    //    for (int i_ind {}; i_ind < ckd.swath.n_indices; ++i_ind) {
    //        const int idx { i * ckd.swath.n_indices + i_ind };
    //        const int& ipix { ckd.swath.pix_indices[i_act][idx] };
    //        // Protection against ipix values that are outside the pixel_mask
    //        // dimension
    //        if ((ipix > 0) && (ipix < pixel_mask.size())) {
    //            if (!pixel_mask[ipix]) {
    //                const double& w { ckd.swath.weights[i_act][idx] };
    //                signal[i] += image[ipix] * w;
    //                const double n { image_stdev[ipix] * w };
    //                stdev[i] += n * n;
    //                total_weight += w;
    //            }
    //        }
    //    }
    //    if (total_weight > 0.0) {
    //        signal[i] /= total_weight;
    //        stdev[i] = std::sqrt(stdev[i] / (total_weight * total_weight));
    //    } else {
    //        mask[i] = true;
    //    }
    //}
}

auto Spectrum::calibrate(const CKD& ckd,
                         const double exposure_time,
                         const int i_act) -> void
{
    //const double exposure_time_inv { 1.0 / exposure_time };
    //for (int i {}; i < static_cast<int>(signal.size()); ++i) {
    //    if (mask[i]) {
    //        signal[i] = fill::d;
    //        stdev[i] = fill::d;
    //    } else {
    //        signal[i] *= ckd.rad.rad[i_act][i] * exposure_time_inv;
    //        stdev[i] *= ckd.rad.rad[i_act][i] * exposure_time_inv;
    //    }
    //}
}

} // namespace tango
