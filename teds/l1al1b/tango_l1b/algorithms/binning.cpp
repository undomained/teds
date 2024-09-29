// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "binning.h"

namespace tango {

//Binning::Binning() {
//}

//Binning::~Binning() {
//}

std::string Binning::getName() const {
    return std::string("Binning");
}

bool Binning::algoCheckInput(const CKD& ckd, L1& l1) {
    if (l1.nr_coadditions == 0) {
        spdlog::warn("nr of coadditions = 0, skipping coaddition step");
        return false;
    } else {
        return true;
    }
    
}

//void Binning::unloadData() {
//    spdlog::info("Binning unload fct still to be filled in");
//}

void Binning::algoExecute(L1& l1, const Dataset& input_data) {

    CKD const& ckd = input_data.get_container<CKD>("ckd");
    BinningTable const& binning = input_data.get_container<BinningTable>("binning");

    binning.bin(l1.image);   
    // For some reason the bin fct devides by binSize. Need to multiply by binsize here.
    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
//    for (int i {}; i < binning_table.nBins(); ++i) {
        l1.image[i] = l1.image[i] * binning.binSize(i);
    }

//        if (getModelType() == "IM"){
//            bin_id = input_data.get_dataset('binning_table_id', c_name='config', group='detector')
//        if (getModelType() == "L1B"){
//            bin_id = int(input_data.get_dataset('binning_table', c_name='measurement', group='image_attributes', kind='variable')[0])
//        table = f"Table_{bin_id}" 
//
//        count_table = input_data.get_dataset('count_table', c_name='binning', group=table, kind='variable')
//        binning_table = input_data.get_dataset('binning_table', c_name='binning', group=table, kind='variable')
//        binned_pixels = input_data.get_dataset('bins', c_name='binning', group=table, kind='dimension')
//        det_cols = input_data.get_dataset('detector_column', c_name='ckd', kind='dimension')
//        det_rows = input_data.get_dataset('detector_row', c_name='ckd', kind='dimension')
//        binned_rows = int(binned_pixels/det_cols)
//        binned_image = np.zeros((binned_rows, self._data.shape[1]))
//
//        # Note: nrs in binning-table are pixel numbers!
//        for det_row in range(binning_table.shape[0]):
//            pix_numbers = binning_table[det_row]
//            binned_row = pix_numbers[0]
//            # Taking a shortcut here. Assuming no binning in col direction.
//            # Binned row number for all columns should be the same.
//            if binned_row != 0:
//                binned_row /= det_cols
//            binned_image[int(binned_row),:] += self._data[det_row,:]
//
//
//    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
//        if (!l1.pixel_mask[i]) {
//            if (getModelType() == "L1B"){
//                l1.stdev[i] *= std::sqrt(l1.nr_coadditions);
//            } else if (getModelType() == "IM"){
//                l1.image[i] *= l1.nr_coadditions;
//            }
//        }
//    }
}

} // namespace tango
