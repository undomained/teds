// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
/// Base class for the algorithms.

#pragma once

#include <string>
#include <map>
#include "../ckd.h"
#include "../l1.h"

namespace tango {

class CKD;
class L1;

/// Each algorithm will inherit from this class and only has to implement
/// the algoCheckInput function to check for necessary data and the algoExecute
/// function which performs the actual algorithm.
///

class BaseAlgo {
public:
    /// Constructor.
    BaseAlgo();

    /// Destructor.
    virtual ~BaseAlgo();

    /// Must be implemented in derived classes
    virtual std::string getName() const;

    /// Virtual function to be implemented by each algorithm which
    /// inherits from BaseAlgo. In this function algorithm-specific code
    /// should be placed which checks for the availability of necessary data
    /// needed by the algorithm.
    virtual void algoCheckInput(const CKD& ckd, L1& l1);

    /// Virtual function to be implemented by each algorithm which
    /// inherits from BaseAlgorithm. In this function algorithm-specific code
    /// should be placed which implements the intended correction.
    virtual void algoExecute(const CKD& ckd, L1& l1);

//    // TODO do we need this unload?
//    /// Set all data pointers to NULL.
//    virtual void unloadData() = 0;

};

} // namespace tango

