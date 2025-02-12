// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "base_algo.h"

namespace tango {

class CKD;
class L1;

/// Dummy Correction class

class DrawOnDetector : public BaseAlgo {
public:

    /// Constructor.
    DrawOnDetector() = default;

    /// Destructor.
    ~DrawOnDetector() = default;

    /// Return the name of the class.
    std::string getName() const override;

    /// Retrieve the required datasets
    bool algoCheckInput(L1& l1, const Dataset& input_data) override;

//    /// Set all loaded data to null.
//    virtusl void unloadData();

    /// Perform the algorithm
    void algoExecute(L1& l1, const Dataset& input_data) override;


};
} // namespace tango
