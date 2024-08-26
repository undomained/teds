// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "prnu_correction.h"
#include "ckd.h"
#include "l1.h"

#include <cstdint>
#include <spdlog/spdlog.h>

namespace tango {

//PRNUCorrection::PRNUCorrection() {
//}

//PRNUCorrection::~PRNUCorrection() {
//}

std::string PRNUCorrection::getName() const {
    return std::string("PRNUCorrection");
}

void PRNUCorrection::algoCheckInput(const CKD& ckd, L1& l1)
{
    spdlog::info("PRNUCorrection algoCheckInput fct still to be filled in");
}

//void PRNUCorrection::unloadData() {
////    m_prnuBinned = NULL;
////    BaseAlgorithm::unloadData();
//    spdlog::info("PRNUCorrection unload fct still to be filled in");
//}

void PRNUCorrection::algoExecute(const CKD& ckd, const bool enabled, L1& l1)
{
    if (!enabled) {
        return;
    }
    spdlog::info("PRNUCorrection algoExecute fct still to be filled in");
//    for (int i {}; i < static_cast<int>(l1.image.size()); ++i) {
//        if (!l1.pixel_mask[i]) {
//            l1.image[i] /= ckd.prnu.prnu[i];
//            l1.stdev[i] /= ckd.prnu.prnu[i];
//        }
//    }
}

} // namespace tango
