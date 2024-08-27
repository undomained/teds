// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing the calibration key data

#pragma once

#include "base_algo.h"

namespace tango {

class CKD;
class L1;


/// PRNUCorrection class

class PRNUCorrection: public BaseAlgo
{
public:

    /// Constructor.
    PRNUCorrection() = default;

    /// Destructor.
    virtual ~PRNUCorrection() = default;

    /// Return the name of the class.
    virtual std::string getName() const;

    /// Retrieve the required datasets
    virtual void algoCheckInput(const CKD& ckd, L1& l1);

//    /// \brief Set all loaded data to null.
//    void unloadData() override;

    /// \brief Perform the algorithm
    /// \param [in] measurement The measurement.
    virtual void algoExecute(const CKD& ckd, const bool enabled, L1& l1);
};

} // namespace tango
