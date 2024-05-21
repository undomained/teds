#pragma once

#include "header.h"
#include "logger.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace tango {

class Program : public Logger {

    public:
    // Constructor.
    Program(
        Logger *creator
    );

    // Execution of the program.
    int execute(
        string settings_file,
        bool foldsettings
    );

};

} // namespace tango
