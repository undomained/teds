// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef TANGO_CAL_H
#define TANGO_CAL_H

#include "header.h"
#include "logger.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

class Tango_cal : public Logger {

    public:
    // Constructor.
    Tango_cal(
        Logger *creator
    );
    ~Tango_cal();

    // Execution of the program.
    int execute(
        string settings_file,
        bool foldsettings
    );

};

#endif
