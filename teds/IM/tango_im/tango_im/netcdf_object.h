#pragma once

#include "header.h"
#include "logger.h"

// Performs a NetCDF action and if it fails, raise error with NetCDF error message.
#define netcdf_check(netcdf_instance,act) \
    try { \
        act; \
    } \
    catch (const NcException &e) { \
        string errormessage = e.what(); \
        size_t netcdffound = errormessage.find("NetCDF: "); \
        if (netcdffound != string::npos) errormessage = errormessage.substr(netcdffound+8); \
        size_t breakfound = errormessage.find("\n"); \
        errormessage = errormessage.substr(0,breakfound); \
        raise_error("NetCDF error: '%s' occurred for file '%s'.",errormessage.c_str(),(netcdf_instance)->filename.c_str()); \
    }

namespace tango {

// NetCDF structure that uses destructor to mention that the file is closed.
// The destructor of the NcFile itself does the actual file closing.
class NetCDF_object : public Logger {
    private:
    bool opened = false;
    bool writeaccess;
    public:
    string filename; // Public because it is needed inside a macro.
    unique_ptr<NcFile> ncid;
    NetCDF_object(
        Logger *creator
    );
    ~NetCDF_object();
    int open(
        string &a_filename,
        NcFile::FileMode fmode
    );
    int close(
        bool strict=true
    );
};

} // namespace tango
