// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "base_algo.h"

namespace tango {

class CKD;
class BinningTable;
class L1;

/// Binning class

class Binning : public BaseAlgo {
public:

    /// Constructor.
    Binning() = default;

    /// Destructor.
    ~Binning() = default;

    /// Return the name of the class.
    std::string getName() const override;

    /// Retrieve the required datasets
    bool algoCheckInput(L1& l1, const Dataset& input_data) override;

//    /// Set all loaded data to null.
//    virtusl void unloadData();

    /// Perform the algorithm
    void algoExecute(L1& l1, const Dataset& input_data) override;

    void binWavelength(L1& l1, BinningTable const& binning);

};
} // namespace tango
