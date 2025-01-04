# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
from netCDF4 import Dataset

import logging
_logger = logging.getLogger(__name__)


def write_parameters_to_nc(parameter_orbit,
                           varname,
                           lat_orbit,
                           lon_orbit,
                           alt_orbit,
                           julday_orbit,
                           npixels,
                           file_name):
    nvar = len(varname)

    _logger.info("Writing parameters {} to file {}".format(varname, file_name))

    da = Dataset(file_name, 'w', format='NETCDF4')

    da.createDimension('Npixels_orbit', npixels)

    da.createVariable('lon', 'f4', ('Npixels_orbit'))
    da.createVariable('lat', 'f4', ('Npixels_orbit'))
    da.createVariable('alt', 'f4', ('Npixels_orbit'))
    da.createVariable('julday', 'f4', ('Npixels_orbit'))

    da.variables['lat'][:] = lat_orbit
    da.variables['lon'][:] = lon_orbit
    da.variables['alt'][:] = alt_orbit
    da.variables['julday'][:] = julday_orbit

    for ivarname in range(0, nvar):
        variablename = varname[ivarname]
        da.createVariable(variablename, 'f4', ('Npixels_orbit'))
        da.variables[variablename][:] = parameter_orbit[ivarname, :]

    da.close()


def write_regridded_parameters_to_nc(parameter_orbit,
                                     varname,
                                     lat_orbit,
                                     lon_orbit,
                                     npixels,
                                     file_name):
    nvar = len(varname)

    da = Dataset(file_name, 'w', format='NETCDF4')

    da.createDimension('Npixels_orbit', npixels)

    da.createVariable('lon', 'f4', ('Npixels_orbit'))
    da.createVariable('lat', 'f4', ('Npixels_orbit'))

    da.variables['lat'][:] = lat_orbit
    da.variables['lon'][:] = lon_orbit

    for ivarname in range(0, nvar):
        variablename = varname[ivarname]
        da.createVariable(variablename, 'f4', ('Npixels_orbit'))
        da.variables[variablename][:] = parameter_orbit[ivarname, :]

    da.close()


def write_pixel_parameter_to_nc(parameter,
                                parameter_name,
                                lat_orbit,
                                lon_orbit,
                                npixels,
                                file_name):
    _logger.info("Writing paramter {} to file {}".format(parameter_name,
                                                         file_name))

    da = Dataset(file_name, 'w', format='NETCDF4')

    da.createDimension('Npixels_orbit', npixels)

    da.createVariable('lon', 'f4', ('Npixels_orbit'))
    da.createVariable('lat', 'f4', ('Npixels_orbit'))
    da.createVariable(parameter_name, 'f4', ('Npixels_orbit'))

    da.variables['lat'][:] = lat_orbit
    da.variables['lon'][:] = lon_orbit
    da.variables[parameter_name][:] = parameter

    da.close()
