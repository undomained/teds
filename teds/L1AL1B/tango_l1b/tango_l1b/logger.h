// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#pragma once

#include "header.h"

namespace tango {

enum loglevel_t {
    log_verbose = 0, // Log on screen only when running with -v.
    log_trace = 1, // Progress, 'Doing this now'. Write on screen.
    log_debug = 2, // Write debugging data. Write just to the file.
    log_info = 3, // Write to screen and file.
    log_warning = 4, // Write to screen and file. Also count the warning.
    log_error = 5, // Write to screen and file with big ERROR. Also, conclude with unsucessful run at the end.
    // At the end, the warnings are counted and the success of the run is evaluated and logged with level info.
};

struct Percentage_struct {
    // Percentage progress logging.
    string percentage_message;
    int percentage;
};

class Logger {

    protected:
    // Constructor for derived objects.
    Logger(
        Logger *creator
    );
    public:
    // Constructor for raw new object.
    Logger(
        string a_filename, // File in which to log.
        bool a_verboseflag, // Flag to write out I/O statements.
        string a_timestamp // Time stamp for any NetCDF file that is created.
    );

    // MPI rank of this process
    int my_rank {};
    // Total number of MPI processes
    int num_procs { 1 };

    // Destructor.
    virtual ~Logger();

    // The main print routine. Only produce output if inactive=false (default).
    void writelog(
        loglevel_t level,
        string fmt,
        ...
    );

    // Sets the name of the logger, in principle, the name of the processor
    // that is running.
    void setName(
        string a_name
    );

    // Progress bar (percentage) routines.
    void percentagelog_open(
        string fmt,
        ...
    );
    void percentagelog_progress(
        size_t idx,
        size_t sz
    );
    void percentagelog_close(
        int percentage_end=100
    );
    auto setInactive(const bool inactive) -> void;

    private:
    Logger *trunk; // Pointer to the original Logger instance, the one that is created with the constructor that does not refer to another Logger.
    bool verboseflag;
    bool dead;
    size_t nwarn;
    string filename;
    string name;
    FILE *file;
    // Wether commands such as writelog produce any output. In a
    // parallel calculation this should be false only for process 0.
    bool inactive { false };

    // Administration of percentage logging.
    shared_ptr<vector<Percentage_struct>> perc;

    protected:
    string timestamp;

};

} // namespace tango
