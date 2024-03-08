// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "git.h"
#include "logger.h"
#include "netcdf_object.h"
#include "settings_main.h"
#include "ckd.h"
#include "dimcal.h"
#include "darkcal.h"
#include "noisecal.h"
#include "nonlincal.h"
#include "prnucal.h"
#include "straycal.h"
#include "fovcal.h"
#include "swathcal.h"
#include "wavecal.h"
#include "radcal.h"
#include "l1b.h"
#include "l1c.h"
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

    // Verify the process list. It must be either:
    // 1. A contiguous block of CKD steps (need not be in right order).
    // 2. Just L1B.
    // 3. Just L1C.
    vector<bool> execution(nlevel,false);
    for (size_t iproc=0 ; iproc<set_main.process.size() ; iproc++) {
        string &processname = set_main.process[iproc];
        if (processname.compare("dim") == 0) execution[LEVEL_DIMCAL] = true;
        else if (processname.compare("dark") == 0) execution[LEVEL_DARKCAL] = true;
        else if (processname.compare("noise") == 0) execution[LEVEL_NOISECAL] = true;
        else if (processname.compare("nonlin") == 0) execution[LEVEL_NONLINCAL] = true;
        else if (processname.compare("prnu") == 0) execution[LEVEL_PRNUCAL] = true;
        else if (processname.compare("stray") == 0) execution[LEVEL_STRAYCAL] = true;
        else if (processname.compare("fov") == 0) execution[LEVEL_FOVCAL] = true;
        else if (processname.compare("swath") == 0) execution[LEVEL_SWATHCAL] = true;
        else if (processname.compare("wave") == 0) execution[LEVEL_WAVECAL] = true;
        else if (processname.compare("rad") == 0) execution[LEVEL_RADCAL] = true;
        else if (processname.compare("l1b") == 0) execution[LEVEL_L1B] = true;
        else if (processname.compare("l1c") == 0) execution[LEVEL_L1C] = true;
        else {
            raise_error("Error: Unrecognized process: %s.",processname.c_str());
        }
    }
    // Now we have the boolean array, check if it fulfills the desires.
    // L1B and L1C should be on their own.
    // For L1B, no CKD output should be written.
    // For L1C, no CKD should exist.
    bool write_ckd;
    if (execution[LEVEL_L1B] || execution[LEVEL_L1C]) {
        write_ckd = false; // For L1C, this does not matter.
        check_error(set_main.process.size() != 1,"Error: L1B and L1C processes cannot be combined with other processes.");
    } else write_ckd = true;
    bool turned_on = false;
    bool turned_off = false;
    level_t level_first = LEVEL_FILLVALUE;
    level_t level_last = LEVEL_FILLVALUE;
    for (size_t ilev=0 ; ilev<nlevel ; ilev++) {
        if (execution[ilev]) {
            check_error(turned_off,"Error: Non-contiguous batch of processor steps selected.");
            if (!turned_on) level_first = (level_t) ilev; // Save first action for CKD reading.
            turned_on = true;
            level_last = (level_t) ilev;
        } else {
            if (turned_on) turned_off = true;
        }
    }
    check_error(!turned_on,"Error: No processor selected at all.");
    check_error(level_first == LEVEL_FILLVALUE,"Program error: Variable level_first not properly initialized during search loop.");
    check_error(level_last == LEVEL_FILLVALUE,"Program error: Variable level_last not properly initialized during search loop.");

    level_t level_ckd;
    check_error(level_first != LEVEL_L1B && set_main.last_calibration_step.compare("") != 0,"Error: Reduced calibration steps only supported for L1B processor.");
    if (set_main.last_calibration_step.compare("") == 0) level_ckd = level_first;
    else if (set_main.last_calibration_step.compare("none") == 0) level_ckd = LEVEL_DIMCAL;
    else if (set_main.last_calibration_step.compare("dim") == 0) level_ckd = LEVEL_DARKCAL;
    else if (set_main.last_calibration_step.compare("dark") == 0) level_ckd = LEVEL_NOISECAL;
    else if (set_main.last_calibration_step.compare("noise") == 0) level_ckd = LEVEL_NONLINCAL;
    else if (set_main.last_calibration_step.compare("nonlin") == 0) level_ckd = LEVEL_PRNUCAL;
    else if (set_main.last_calibration_step.compare("prnu") == 0) level_ckd = LEVEL_STRAYCAL;
    else if (set_main.last_calibration_step.compare("stray") == 0) level_ckd = LEVEL_FOVCAL;
    else if (set_main.last_calibration_step.compare("fov") == 0) level_ckd = LEVEL_SWATHCAL;
    else if (set_main.last_calibration_step.compare("swath") == 0) level_ckd = LEVEL_WAVECAL;
    else if (set_main.last_calibration_step.compare("wave") == 0) level_ckd = LEVEL_RADCAL;
    else if (set_main.last_calibration_step.compare("rad") == 0) level_ckd = LEVEL_L1B;
    else {
        raise_error("Level %s not recognized.",set_main.last_calibration_step.c_str());
    }

    // Construct CKD. Read actions not in constructor because of possible errors.
    CKD ckd(this);
    // Do not read for L1C. If L1C is selected, it is the only process.
    if (!execution[LEVEL_L1C]) handle(ckd.read(set_main,level_ckd,write_ckd));

    for (int lev=(int)level_first ; lev<=(int)level_last ; lev++) {
        unique_ptr<Processor> proc;
        if (lev == LEVEL_DIMCAL) proc = make_unique<Dimcal>(this,&ckd);
        if (lev == LEVEL_DARKCAL) proc = make_unique<Darkcal>(this,&ckd);
        if (lev == LEVEL_NOISECAL) proc = make_unique<Noisecal>(this,&ckd);
        if (lev == LEVEL_NONLINCAL) proc = make_unique<Nonlincal>(this,&ckd);
        if (lev == LEVEL_PRNUCAL) proc = make_unique<Prnucal>(this,&ckd);
        if (lev == LEVEL_STRAYCAL) proc = make_unique<Straycal>(this,&ckd);
        if (lev == LEVEL_FOVCAL) proc = make_unique<Fovcal>(this,&ckd);
        if (lev == LEVEL_SWATHCAL) proc = make_unique<Swathcal>(this,&ckd);
        if (lev == LEVEL_WAVECAL) proc = make_unique<Wavecal>(this,&ckd);
        if (lev == LEVEL_RADCAL) proc = make_unique<Radcal>(this,&ckd);
        if (lev == LEVEL_L1B) proc = make_unique<L1B>(this,&ckd);
        if (lev == LEVEL_L1C) proc = make_unique<L1C>(this,(CKD*)NULL);
        handle(proc->execute(settings_file,&set_main));
    }
    return 0;

} // }}}

