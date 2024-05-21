// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "netcdf_object.h"

#include <numeric>

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
    writeaccess = fmode != NcFile::read;
    if (writeaccess && my_rank > 0) {
        // First process creates the correct output files. Ohter
        // processes create dummy files which are later deleted. The
        // reason is that all processes need to refer to variables
        // such as var_time but netCDF variables cannot be created
        // without first creating and opening a file.
        filename += std::to_string(my_rank);
    }
    netcdf_check(this,ncid = make_unique<NcFile>(filename.c_str(),fmode));
    opened = true;
    if (writeaccess) {
        netcdf_check(this,ncid->putAtt("Conventions","CF-1.6"));
        netcdf_check(this,ncid->putAtt("institution","SRON Netherlands Institute for Space Research"));
        netcdf_check(this,ncid->putAtt("instrument","TANGO"));
        netcdf_check(this,ncid->putAtt("product_name",filename));
        netcdf_check(this,ncid->putAtt("project","TANGO Project"));
        //netcdf_check(this,ncid->putAtt("creator_email","SPEXone-MPC@sron.nl"));
        netcdf_check(this,ncid->putAtt("creator_name","SRON/Earth Science"));
        netcdf_check(this,ncid->putAtt("creator_url","https://www.sron.nl/missions-earth"));
        //netcdf_check(this,ncid->putAtt("publisher_email","SPEXone-MPC@sron.nl"));
        netcdf_check(this,ncid->putAtt("publisher_name","SRON/Earth Science"));
        netcdf_check(this,ncid->putAtt("publisher_url","https://www.sron.nl/missions-earth"));
        netcdf_check(this,ncid->putAtt("date_created",now_timestring()));
        //netcdf_check(this,ncid->putAtt("git_commit",git_commit));
    }
    return 0;
} // }}}

auto NetCDF_object::copyAtts(const netCDF::NcVar& in,
                             netCDF::NcVar& out) -> void
{
    void* value_buf[32768] {};
    for (const auto& [name, att] : in.getAtts()) {
        att.getValues(value_buf);
        out.putAtt(name, att.getType(), att.getAttLength(), value_buf);
    }
}

auto NetCDF_object::copyVar(const netCDF::NcVar& var,
                            netCDF::NcGroup& group,
                            const std::vector<netCDF::NcDim>& dims)
  -> netCDF::NcVar
{
    auto var_out { group.addVar(
          var.getName(), var.getType(), dims.empty() ? var.getDims() : dims) };
    copyAtts(var, var_out);
    std::vector<size_t> dims_int {};
    for (const auto& dim : var.getDims()) {
        dims_int.push_back(dim.getSize());
    }
    size_t total_size { std::accumulate(dims_int.cbegin(),
                                        dims_int.cend(),
                                        static_cast<size_t>(1),
                                        std::multiplies<size_t>()) };
    void* values[var.getType().getSize() * total_size];
    var.getVar(values);
    var_out.putVar(values);
    return var_out;
}

int NetCDF_object::close( // {{{
    bool strict
)
{
    check_error(strict && !opened,"Error: Closing already closed file: '%s'.",filename.c_str());
    if (opened) {
        if (writeaccess) {
            netcdf_check(this,ncid->putAtt("date_finalized",now_timestring()));
            if (my_rank > 0) {
                std::remove(filename.c_str());
            }
        }
        writelog(log_verbose,"Closing NetCDF file: '%s'.",filename.c_str());
        ncid->close();
        opened = false;
    }
    return 0;
} // }}}

