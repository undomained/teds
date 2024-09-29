// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing the calibration key data

#pragma once

#include "base_algo.h"

namespace tango {

class CKD;
class L1;


/// PRNU class

class PRNU: public BaseAlgo
{
public:

    /// Constructor.
    PRNU() = default;

    /// Destructor.
    ~PRNU() = default;

    /// Return the name of the class.
    std::string getName() const override;

    /// Retrieve the required datasets
    bool algoCheckInput(const CKD& ckd, L1& l1) override;

//    /// \brief Set all loaded data to null.
//    void unloadData() override;

    /// \brief Perform the algorithm
    /// \param [in] measurement The measurement.
    void algoExecute(L1& l1, const Dataset& input_data) override;

};

} // namespace tango
