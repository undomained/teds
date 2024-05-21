// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "header.h"
#include "logger.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace tango {

class Spexone_cal : public Logger {

    public:
    // Constructor.
    Spexone_cal(
        Logger *creator
    );
    ~Spexone_cal();

    // Execution of the program.
    int execute(
        string settings_file,
        bool foldsettings
    );

};

} // namespace tango
