#include "header.h"
#include "git.h"
#include "functions.h"
#include "logger.h"
#include "netcdf_object.h"

// Constructor, only taking over the logging.
NetCDF_object::NetCDF_object( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}

// NetCDF structure that uses destructor to mention that the file is closed.
// The destructor of the NcFile itself does the actual file closing.
NetCDF_object::~NetCDF_object( // {{{
)
{
    if (close(false) != 0) {
        writelog(log_error,"Error during destructor of NetCDF_object.");
        writelog(log_error,"this error escapes the regular error handling procedure.");
        error_message("Destructor error raised");
    }
} // }}}

int NetCDF_object::open( // {{{
    string &a_filename,
    NcFile::FileMode fmode
)
{
    if (fmode == NcFile::replace || fmode == NcFile::newFile) {
        // Add time stamp to file.
        size_t dotpos = a_filename.find(".");
        if (dotpos == string::npos) {
            filename = a_filename + timestamp;
        } else {
            filename = a_filename.substr(0,dotpos) + timestamp + a_filename.substr(dotpos,a_filename.length()-dotpos);
        }
    } else filename = a_filename;
    writelog(log_verbose,"Opening NetCDF file: '%s'.",filename.c_str());
    netcdf_check(this,ncid = make_unique<NcFile>(filename.c_str(),fmode));
    writeaccess = fmode != NcFile::read;
    opened = true;
    if (writeaccess) {
        netcdf_check(this,ncid->putAtt("Conventions","CF-1.6"));
        netcdf_check(this,ncid->putAtt("institution","SRON Netherlands Institute for Space Research"));
        netcdf_check(this,ncid->putAtt("instrument","TANGO"));
        netcdf_check(this,ncid->putAtt("processing_version",git_tag));
        netcdf_check(this,ncid->putAtt("product_name",filename));
        netcdf_check(this,ncid->putAtt("project","TANGO Project"));
        //netcdf_check(this,ncid->putAtt("creator_email","@sron.nl"));
        netcdf_check(this,ncid->putAtt("creator_name","SRON/Earth Science"));
        netcdf_check(this,ncid->putAtt("creator_url","https://www.sron.nl/missions-earth"));
        //netcdf_check(this,ncid->putAtt("publisher_email","@sron.nl"));
        netcdf_check(this,ncid->putAtt("publisher_name","SRON/Earth Science"));
        netcdf_check(this,ncid->putAtt("publisher_url","https://www.sron.nl/missions-earth"));
        netcdf_check(this,ncid->putAtt("date_created",now_timestring()));
        //netcdf_check(this,ncid->putAtt("git_commit",git_commit));
    }
    return 0;
} // }}}

int NetCDF_object::close( // {{{
    bool strict
)
{
    check_error(strict && !opened,"Error: Closing already closed file: '%s'.",filename.c_str());
    if (opened) {
        if (writeaccess) {
            netcdf_check(this,ncid->putAtt("date_finalized",now_timestring()));
        }
        writelog(log_verbose,"Closing NetCDF file: '%s'.",filename.c_str());
        ncid->close();
        opened = false;
    }
    return 0;
} // }}}

