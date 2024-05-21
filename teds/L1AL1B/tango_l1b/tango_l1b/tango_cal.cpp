// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "git.h"
#include "logger.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "ckd.h"
#include "l1b.h"
#include "tango_cal.h"

Tango_cal::Tango_cal( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
Tango_cal::~Tango_cal() {}

int Tango_cal::execute( // {{{
    string settings_file,
    bool foldsettings
)
{
    writelog(log_info,"%s %s (%s)", "tango_cal", git_tag.c_str(), git_date.c_str());
    writelog(log_debug,"Commit: %s",git_commit.c_str());
    if (num_procs > 1) {
        writelog(log_info, "Number of MPI processes: %i", num_procs);
    }
    // Log the settings into the log file.
    // This scope makes sure to clean any I/O-related stuff from the stack.
    // It would remain there until the end of the program.
    {
        struct stat settingsfile_info;
        check_error(stat(settings_file.c_str(),&settingsfile_info),"Error acquiring file information from '%s'.",settings_file.c_str());
        string modtime = ctime(&settingsfile_info.st_mtime);
        size_t found = modtime.find("\n");
        modtime = modtime.substr(0,found);
        writelog(log_debug,"Settings file: '%s' (Last modified %s).%s",settings_file.c_str(),modtime.c_str(),foldsettings?" {{{":"");
        ifstream setstream(settings_file.c_str());
        check_error(setstream.fail(),"Error opening settings file: '%s'.",settings_file.c_str());
        string line;
        while (getline(setstream,line)) {
            writelog(log_debug,line.c_str());
        }
        writelog(log_debug,"End of settings file.%s",foldsettings?" }}}":"");
    }

    // Read logistic (main) settings.
    Settings_main set_main(this);
    handle(set_main.init(settings_file));

    if (set_main.process.size() != 1){
        raise_error("Error: More than one process defined. This is not possible for Tango E2E \n");
    }
    string &processname = set_main.process[0];
    if (processname.compare("l1b") != 0) {
        raise_error("For Tango the process name can only be l1b! Process name %s not possible.",processname.c_str());
    }

    level_t level_ckd = LEVEL_L1B;
    bool write_ckd = false;

    // Construct CKD. Read actions not in constructor because of possible errors.
    CKD ckd(this);
    handle(ckd.read(set_main,level_ckd,write_ckd));

    unique_ptr<Processor> proc;
    proc = make_unique<L1B>(this,&ckd);

    handle(proc->execute(settings_file,&set_main));

    return 0;

} // }}}

