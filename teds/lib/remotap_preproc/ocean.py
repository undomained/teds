# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
from netCDF4 import Dataset
import numpy as np
import time as tm

from .collocation_algorithm import interpolate_to_orbit
from .exceptions import ProcessError
import logging

_logger = logging.getLogger(__name__)


class Ocean(object):

    def __init__(self, path):
        self.source = "MiPrep processing from 1*1 degree resolution ocean data"
        self.description = "Ocean data collocation"
        self.chl_filename = path + "chl_regrid.nc"
        self.ws_filename = path + "windspeed.nc"
        self.output_ocean_prefix = 'xocean'

        self.varname_ocean = ['Xchl', 'uwnd', 'vwnd']
        self.npar = 3
        self.ocean_para_in = None
        self.lons_ocean = None
        self.lats_ocean = None

        self.ocean_orbit = None
        self.lat_orbit = None
        self.lon_orbit = None
        self.n_pixels = None

    def read_file(self):

        _logger.info("Reading {}".format(self.chl_filename))
        root_grp = Dataset(self.chl_filename)

        self.lons_ocean = root_grp.variables["longitude"][:]
        self.lats_ocean = root_grp.variables["latitude"][:]
        xchl = root_grp.variables["xchl"][:][:]  # (xchl) 180(lat)*360(lon)
        nlat_ocean, nlon_ocean = xchl.shape
        self.ocean_para_in = np.zeros((self.npar, nlat_ocean, nlon_ocean))
        # (self.ocean_para_in[0,:,:]) 180(lat)*360(lon)
        self.ocean_para_in[0, :, :] = xchl[:, :]
        root_grp.close()

        root_grp = Dataset(self.ws_filename)
        # (self.ocean_para_in[1,:,:]) 180(lat)*360(lon)
        self.ocean_para_in[1, :, :] = root_grp.variables['uwnd'][:][:]
        # (self.ocean_para_in[2,:,:]) 180(lat)*360(lon)
        self.ocean_para_in[2, :, :] = root_grp.variables['vwnd'][:][:]
        root_grp.close()

    def collocate(self, julday_orbit, lat, lon, n_pixels):
        self.lat_orbit = lat
        self.lon_orbit = lon
        self.n_pixels = n_pixels

        if self.ocean_para_in is None:
            self.read_file()

        _logger.info("Collocation of ocean data.")

        self.ocean_orbit = interpolate_to_orbit(self.lat_orbit,
                                                self.lon_orbit,
                                                self.n_pixels,
                                                self.lats_ocean,
                                                self.lons_ocean,
                                                self.ocean_para_in,
                                                self.npar)

    def write_output(self, output_dir, alt_orbit, julday_orbit, add_grid=True):

        if self.ocean_orbit is None:
            raise ProcessError("Collocation on orbit has not been successful.")

        _logger.info('Write (orbital) xocean to ncdf')

        # write Ocean orbit data to netcdf file
        oceanPara_file = output_dir + self.output_ocean_prefix + '.nc'
        root_grp = Dataset(oceanPara_file, 'w', format='NETCDF4')

        root_grp.source = ""
        root_grp.title = "Collocated ocean data"
        root_grp.institution = "SRON Netherlands Institute for Space Research"
        root_grp.date_created = tm.ctime(tm.time())

        root_grp.createDimension('Npixels_orbit', self.n_pixels)

        root_grp.createVariable('Xchl', 'f4', ('Npixels_orbit'))
        root_grp.variables['Xchl'][:] = self.ocean_orbit[0, :]

        root_grp.createVariable('uwnd', 'f4', ('Npixels_orbit'))
        root_grp.variables['uwnd'][:] = self.ocean_orbit[1, :]

        root_grp.createVariable('vwnd', 'f4', ('Npixels_orbit'))
        root_grp.variables['vwnd'][:] = self.ocean_orbit[2, :]

        if add_grid:
            root_grp.createVariable('lon', 'f4', ('Npixels_orbit'))
            root_grp.createVariable('lat', 'f4', ('Npixels_orbit'))
            root_grp.variables['lat'][:] = self.lat_orbit
            root_grp.variables['lon'][:] = self.lon_orbit

        root_grp.close()

    def write_to_group(self, group, start=0, end=None):
        _logger.info('Write ocean parameter to ncdf')

        for ivar in range(self.npar):
            var = self.varname_ocean[ivar]
            group.variables[var][start:end] = self.ocean_orbit[ivar, :]
