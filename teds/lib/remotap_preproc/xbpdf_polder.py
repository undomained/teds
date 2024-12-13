from netCDF4 import Dataset
import numpy as np

from .exceptions import ProcessError
from .collocation_algorithm import interpolate_single_parameter_to_orbit
from .write_module import write_pixel_parameter_to_nc

import logging
_logger = logging.getLogger(__name__)


class XbpdfPolder(object):
    def __init__(self, path, year, month):
        self.source = ""
        self.description = ""
        self.path = path
        self.year = year
        self.month = month

        self.lats_xbpdf = None
        self.lons_xbpdf = None
        self.xbpdf_seasonal_mean_in = None

        self.lat_orbit = None
        self.lon_orbit = None
        self.n_pixels = None

    def read_file(self):
        i_month = int(self.month)
        if 3 >= i_month >= 1:
            month_polder_xbpdf = ['01', '02', '03']
        elif 6 >= i_month >= 4:
            month_polder_xbpdf = ['04', '05', '06']
        elif 9 >= i_month >= 7:
            month_polder_xbpdf = ['07', '08', '09']
        elif 12 >= i_month >= 10:
            month_polder_xbpdf = ['10', '11', '12']
        else:
            raise ProcessError('Something wrong with date_id_short for polder xbpdf')

        # year_polder_xbpdf = '2006'
        fname_Xbpdf = self.path + self.year + '/SRON_parasol_gridded' + self.year + \
            month_polder_xbpdf[0] + '.nc'

        # xbpdf: 2*2 degree resolution data
        # Take seasonal mean.

        _logger.info('Reading Xbpdf polder file.')

        root_grp = Dataset(fname_Xbpdf)
        xbpdf_tmp = root_grp.variables['xbpdf'][:][:][:]  # (xbpdf_tmp) 30(days)*180(lon)*90(lat)
        Ndays, Nlon_xbpdf, Nlat_xbpdf = xbpdf_tmp.shape

        xbpdf = np.zeros((3, Nlon_xbpdf, Nlat_xbpdf))
        self.lats_xbpdf = root_grp.variables['lat_center'][:]
        self.lons_xbpdf = root_grp.variables['lon_center'][:]
        xbpdf[0, :, :] = np.nanmean(xbpdf_tmp, axis=0)  # (xbpdf[0,:,:]) 180(lon)*90(lat)
        root_grp.close()

        for i in range(1, len(month_polder_xbpdf)):
            fname_Xbpdf = self.path + self.year + '/SRON_parasol_gridded' + self.year + \
                month_polder_xbpdf[i] + '.nc'
            root_grp = Dataset(fname_Xbpdf)
            xbpdf_tmp = root_grp.variables['xbpdf'][:][:][:]  # (xbpdf_tmp) 31(days)*180(lon)*90(lat)
            xbpdf[i, :, :] = np.nanmean(xbpdf_tmp, axis=0)  # (xbpdf[1,:,:]) 180(lon)*90(lat)
            root_grp.close()

        self.xbpdf_seasonal_mean_in = np.nanmean(xbpdf[:, :, :], axis=0)  # (Xbpdf_seasonal_mean) 180(lon)*90(lat)
        self.xbpdf_seasonal_mean_in = np.transpose(self.xbpdf_seasonal_mean_in)  # (Xbpdf_seasonal_mean) 90(lat)*180(lon)

    def collocate(self, julday_orbit, lat, lon, n_pixels):

        if self.xbpdf_seasonal_mean_in is None:
            self.read_file()

        _logger.info("Collocation of xbpdf parameters.")
        self.lat_orbit = lat
        self.lon_orbit = lon
        self.n_pixels = n_pixels
        self.xbpdf_orbit = interpolate_single_parameter_to_orbit(lat, lon, n_pixels, self.lats_xbpdf, self.lons_xbpdf, self.xbpdf_seasonal_mean_in)

    def write_output(self, output_dir, alt_orbit, julday):
        _logger.info("write Xbpdf orbit data to netcdf file")
        fname_xbpdf_prefix = 'xbpdf'
        xbpdf_file = output_dir + fname_xbpdf_prefix+'.nc'

        write_pixel_parameter_to_nc(self.xbpdf_orbit, "Xbpdf", self.lat_orbit, self.lon_orbit, self.n_pixels, xbpdf_file)

    def write_to_group(self, group, start=0, end=None):

        group.variables['Xbpdf'][start:end] = self.xbpdf_orbit
