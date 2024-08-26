// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "base_algo.h"

namespace tango {

class CKD;
class L1;

/// Dummy Correction class

class DummyCorrection : public BaseAlgo {
public:

    /// Constructor.
    DummyCorrection() = default;

    /// Destructor.
    ~DummyCorrection() = default;

    /// Return the name of the class.
    virtual std::string getName() const;

    /// Retrieve the required datasets
    virtual void algoCheckInput(const CKD& ckd, L1& l1);

//    /// Set all loaded data to null.
//    virtusl void unloadData();

    /// Perform the algorithm
    virtual void algoExecute(const CKD& ckd, const bool enabled, L1& l1);

};
} // namespace tango
