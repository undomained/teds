// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "base_algo.h"

namespace tango {

class CKD;
class L1;

/// Dark Offset class

class DarkCurrent : public BaseAlgo {
public:

    /// Constructor.
    DarkCurrent() = default;

    /// Destructor.
    ~DarkCurrent() = default;

    /// Return the name of the class.
    std::string getName() const override;

    /// Retrieve the required datasets
    bool algoCheckInput(L1& l1, const Dataset& input_data) override;

//    /// Set all loaded data to null.
//    virtusl void unloadData() override;

    /// Perform the algorithm
    void algoExecute(L1& l1, const Dataset& input_data) override;

};
} // namespace tango
