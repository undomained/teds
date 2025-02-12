// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
/// Base class for the algorithms.

#pragma once

#include <string>
#include <map>
#include "../datasets.h"
#include "../ckd.h"
#include "../l1.h"
#include "../binning_table.h"
#include <spdlog/spdlog.h>

namespace tango {

class Dataset;
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
    virtual bool algoCheckInput( L1& l1, const Dataset& input_data);

    /// Virtual function to be implemented by each algorithm which
    /// inherits from BaseAlgorithm. In this function algorithm-specific code
    /// should be placed which implements the intended correction.
    virtual void algoExecute(L1& l1, const Dataset& input_data);

//    // TODO do we need this unload?
//    /// Set all data pointers to NULL.
//    virtual void unloadData() = 0;

    /// Getter for model_type
    std::string getModelType() const {
        return model_type;
    }

    /// Setter for model_type
    void setModelType(const std::string& modelType) {
        if (modelType != "L1B" && modelType != "IM") {
            spdlog::error("Model type must be 'L1B' or 'IM'");
        } else {
            model_type = modelType;
        }
    }

private:
    std::string model_type {"L1B"};
};

} // namespace tango

