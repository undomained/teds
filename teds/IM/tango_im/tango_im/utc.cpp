#include "header.h"
#include "logger.h"
#include "netcdf_object.h"
#include "utc.h"

UTC::UTC( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}

UTC::~UTC() {}

int UTC::init( // {{{
    string &utcfile
)
{

    // Read UTC file.
    NetCDF_object nc_utc(this);
    handle(nc_utc.open(utcfile,NcFile::read));
    netcdf_check(&nc_utc,dim_utc = nc_utc.ncid->getDim("sz").getSize());
    mjd_utc.resize(dim_utc);
    tdiff_utc.resize(dim_utc);
    netcdf_check(&nc_utc,nc_utc.ncid->getVar("mjd").getVar(mjd_utc.data()));
    netcdf_check(&nc_utc,nc_utc.ncid->getVar("timediff").getVar(tdiff_utc.data()));

    return 0;

} // }}}

