// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Class for storing the calibration key data

#pragma once

#include "base_algo.h"

namespace tango {

/// \class PRNUCorrection
/// \brief Correct the signal for the PRNU

class PRNUCorrection: public: BaseAlgo
{
public:

    /// \brief Constructor.
    PRNUCorrection();

    /// \brief Destructor.
    virtual ~PRNUCorrection();

    /// Return the name of the class.
    /// \return std::string Name of the class, i.e. "PRNUCorrection".
    virtual std::string getName() const;

    /// Retrieve the required datasets
    virtual void algoCheckInput(const CKD& ckd, const bool enabled, L1& l1);

    /// \brief Set all loaded data to null.
    virtual void unloadData();

    /// \brief Perform the algorithm
    /// \param [in] measurement The measurement.
    virtual void algoExecute(const CKD& ckd, const bool enabled, L1& l1);

} // namespace tango
