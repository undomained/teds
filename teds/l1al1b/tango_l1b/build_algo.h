// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.
//
/// Build class for the algorithms.

#pragma once

#include <string>

#include "base_algo.h"
#include "prnu_correction.h"
#include "dummy_correction.h"

namespace tango {

//class PRNUCorrection;
//class DummyCorrection;

class BuildAlgo
{
public:
    /// Constructor.
    BuildAlgo() = default;

    /// Destructor.
    ~BuildAlgo() = default;

    BaseAlgo* BuildMethod(){
        if(name == "DummyCorrection") {
            return new DummyCorrection();
        }
        else if(name == "PRNUCorrection") {
            return new PRNUCorrection();
        } else {
            return nullptr;
        }
    }

    BaseAlgo* CreateAlgo(std::string name) {
        this->name = name;
        BaseAlgo* ptr = this->BuildMethod();
        return ptr;
    }

    std::string name;
};

} // namespace tango


