# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
# =============================================================================
#     geophysical scene generation module for different E2E simulator profiles
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================

from netCDF4 import Dataset
import numpy as np

from .collocation_algorithm import interpolate_to_orbit
import logging

_logger = logging.getLogger(__name__)

# This module uses a directory containing the following files:
# 'brdf_joost_jan.nc'
# 'brdf_joost_apr.nc'
# 'brdf_joost_july.nc'
# 'brdf_joost_oct.nc'


class XbrdfGomeModis(object):

    def __init__(self, path, month, combined_with_sentinel=False):
        self.source = "MiPrep"
        self.description = "XBRDF from GOME and MODIS"
        self.brdf_gome_modis_path = path
        self.month = month
        self.combined_with_sentinel = combined_with_sentinel

        self.fname_xbrdf_wave_prefix = 'xbrdf_wave'
        self.fname_xbrdf_rli_prefix = 'xbrdf_rli'
        self.fname_xbrdf_wave_combine_prefix = 'xbrdf_wave_combine'

        self.brdf_filename = None

        self.waves_gome2_bands = None
        self.xbrdf_wave_in = None

        self.waves_combined_bands = None
        self.xbrdf_wave_combine_in = None
        self.xbrdf_rli_in = None

        self.lats_surf = None
        self.lons_surf = None
        self.nbands_GOME2 = None

        self.lat_orbit = None
        self.lon_orbit = None
        self.n_pixels = None

        self.npar_rli = 2

        self.xbrdf_wave_orbit = None
        self.xbrdf_wave_combine_orbit = None
        self.xbrdf_rli_orbit = None

    def read_file(self):
        # Select the correct file based on the specified month.
        if 3 >= int(self.month) >= 1:
            self.brdf_filename = (
                self.brdf_gome_modis_path + 'brdf_joost_jan.nc')
        elif 6 >= int(self.month) >= 4:
            self.brdf_filename = (
                self.brdf_gome_modis_path + 'brdf_joost_apr.nc')
        elif 9 >= int(self.month) >= 7:
            self.brdf_filename = (
                self.brdf_gome_modis_path + 'brdf_joost_july.nc')
        elif 12 >= int(self.month) >= 10:
            self.brdf_filename = (
                self.brdf_gome_modis_path + 'brdf_joost_oct.nc')
        else:
            raise Exception(FileNotFoundError)

        root_grp = Dataset(self.brdf_filename)

        self.lons_surf = root_grp.variables['longitude'][:]
        self.lats_surf = root_grp.variables['latitude'][:]

        # Select bands for GOME2
        # (waves_gome2) 21(wave)
        self.waves_gome2_bands = root_grp.variables['bands_gome2'][:]
        # Select data for GOME2
        # (xbrdf_wave) 21(wave)*180(lat)*360(lon)
        self.xbrdf_wave_in = root_grp.variables['min_ler'][:][:][:]
        # Determine the shape
        self.nbands_GOME2, nlat_surf, nlon_surf = self.xbrdf_wave_in.shape

        ###################
        # (waves_gome2) 21(wave)
        self.waves_combined_bands = root_grp.variables['bands_gome2'][:]
        # (waves_modis) 4(wave)
        self.waves_combined_bands = np.append(
            self.waves_combined_bands, root_grp.variables['bands_modis'][1])
        self.waves_combined_bands = np.append(
            self.waves_combined_bands, root_grp.variables['bands_modis'][4:7])

        # (xbrdf_wave_combine) 25(wave)*180(lat)*360(lon)
        self.xbrdf_wave_combine_in = root_grp.variables['min_ler'][:][:][:]
        xbrdf_wave_tmp = np.zeros((4, nlat_surf, nlon_surf))
        xbrdf_wave_tmp[0, :, :] = root_grp.variables['brdf_par1_band2'][:][:]
        xbrdf_wave_tmp[1, :, :] = root_grp.variables['brdf_par1_band5'][:][:]
        xbrdf_wave_tmp[2, :, :] = root_grp.variables['brdf_par1_band6'][:][:]
        xbrdf_wave_tmp[3, :, :] = root_grp.variables['brdf_par1_band7'][:][:]
        self.xbrdf_wave_combine_in = np.append(
            self.xbrdf_wave_combine_in, xbrdf_wave_tmp, axis=0)

        # Nbands_combine = len(self.waves_combined_bands)
        # ('Nbands_combine = ', Nbands_combine)
        _logger.info("BRDF read for different bands.")

        ###################
        self.xbrdf_rli_in = np.zeros((self.npar_rli, nlat_surf, nlon_surf))
        # (Xbrdf_rli[0,:,:]) 180(lat)*360(lon)
        self.xbrdf_rli_in[0, :, :] = (
            root_grp.variables['brdf_LiSparse_band4'][:][:])
        # (Xbrdf_rli[1,:,:]) 180(lat)*360(lon)
        self.xbrdf_rli_in[1, :, :] = (
            root_grp.variables['brdf_RossThick_band4'][:][:])

        root_grp.close()

    def get_nbands(self):
        if self.xbrdf_wave_combine_in is None:
            self.read_file()
        return len(self.waves_combined_bands)

    def collocate(self,
                  julday_orbit,
                  lat,
                  lon,
                  n_pixels,
                  add_xbrdf_wave=True,
                  add_xbdrf_wave_combine=True,
                  add_xbrdf_rli=True):

        if self.xbrdf_wave_combine_in is None:
            self.read_file()

        _logger.info("Collocation of xbrdf file started.")

        self.lat_orbit = lat
        self.lon_orbit = lon
        self.n_pixels = n_pixels

        if add_xbrdf_rli:
            self.xbrdf_rli_orbit = interpolate_to_orbit(self.lat_orbit,
                                                        self.lon_orbit,
                                                        n_pixels,
                                                        self.lats_surf,
                                                        self.lons_surf,
                                                        self.xbrdf_rli_in,
                                                        self.npar_rli)

        if add_xbrdf_wave:
            self.xbrdf_wave_orbit = interpolate_to_orbit(self.lat_orbit,
                                                         self.lon_orbit,
                                                         n_pixels,
                                                         self.lats_surf,
                                                         self.lons_surf,
                                                         self.xbrdf_wave_in,
                                                         self.nbands_GOME2)

        if add_xbdrf_wave_combine:
            self.xbrdf_wave_combine_orbit = interpolate_to_orbit(
                self.lat_orbit,
                self.lon_orbit,
                n_pixels,
                self.lats_surf,
                self.lons_surf,
                self.xbrdf_wave_combine_in,
                len(self.waves_combined_bands))

    def write_output(self,
                     output_dir,
                     alt_orbit,
                     julday_orbit,
                     add_xbrdf_wave=True,
                     add_xbdrf_wave_combine=True,
                     add_xbrdf_rli=True):

        if add_xbrdf_rli:
            _logger.debug('now Write (orbital) xbrdf to ncdf')
            xbrdf_rli_file = output_dir + self.fname_xbrdf_rli_prefix + '.nc'
            self.write_to_nc_xbrdf(self.xbrdf_rli_orbit,
                                   self.lat_orbit,
                                   self.lon_orbit,
                                   self.n_pixels,
                                   xbrdf_rli_file)

        if add_xbrdf_wave:
            _logger.debug('now Write (orbital) xbrdf_wave (GOME2) to ncdf')
            xbrdf_file = output_dir + self.fname_xbrdf_wave_prefix + '.nc'

            self.write_to_nc_xbrdfwave(self.xbrdf_wave_orbit,
                                       self.lat_orbit,
                                       self.lon_orbit,
                                       self.n_pixels,
                                       self.nbands_GOME2,
                                       self.waves_gome2_bands,
                                       xbrdf_file)

        if add_xbdrf_wave_combine:
            _logger.debug('now Write (orbital) xbrdf_wave_combine (GOME2 and '
                          'MODIS) to ncdf')
            xbrdf_combine_file = (
                output_dir + self.fname_xbrdf_wave_combine_prefix + '.nc')
            self.write_to_nc_xbrdfwave_combine(self.xbrdf_wave_combine_orbit,
                                               self.lat_orbit,
                                               self.lon_orbit,
                                               self.n_pixels,
                                               len(self.waves_combined_bands),
                                               self.waves_combined_bands,
                                               xbrdf_combine_file)

    def write_to_nc_xbrdfwave(self,
                              Xbrdf_wave_orbit,
                              lat_orbit,
                              lon_orbit,
                              Npixels,
                              Nbands_GOME2,
                              waves_gome2,
                              xbrdf_file):
        da = Dataset(xbrdf_file, 'w', format='NETCDF4')

        da.createDimension('Nbands_GOME2', Nbands_GOME2)
        da.createDimension('Npixels_orbit', Npixels)

        da.createVariable('Bands_GOME2', 'f4', ('Nbands_GOME2'))
        da.variables['Bands_GOME2'][:] = waves_gome2

        da.createVariable(
            'Xbrdf_wave_GOME2', 'f4', ('Nbands_GOME2', 'Npixels_orbit'))

        da.variables['Xbrdf_wave_GOME2'][:] = Xbrdf_wave_orbit

        da.createVariable('lon', 'f4', ('Npixels_orbit'))
        da.createVariable('lat', 'f4', ('Npixels_orbit'))
        da.variables['lat'][:] = lat_orbit
        da.variables['lon'][:] = lon_orbit

        da.close()

    def write_to_nc_xbrdfwave_combine(self,
                                      Xbrdf_wave_orbit,
                                      lat_orbit,
                                      lon_orbit,
                                      Npixels,
                                      Nbands_combine,
                                      waves_combine,
                                      xbrdf_file):
        da = Dataset(xbrdf_file, 'w', format='NETCDF4')

        da.createDimension('Nbands_combine', Nbands_combine)
        da.createDimension('Npixels_orbit', Npixels)

        da.createVariable('Bands_combine', 'f4', ('Nbands_combine'))
        da.variables['Bands_combine'][:] = waves_combine

        da.createVariable(
            'Xbrdf_wave_combine', 'f4', ('Nbands_combine', 'Npixels_orbit'))

        da.variables['Xbrdf_wave_combine'][:] = Xbrdf_wave_orbit

        da.createVariable('lon', 'f4', ('Npixels_orbit'))
        da.createVariable('lat', 'f4', ('Npixels_orbit'))
        da.variables['lat'][:] = lat_orbit
        da.variables['lon'][:] = lon_orbit

        da.close()

    # write Xbrdf data to netcdf file
    def write_to_nc_xbrdf(self,
                          Xbrdf_rli_orbit,
                          lat_orbit,
                          lon_orbit,
                          Npixels,
                          xbrdf_rli_file):
        da = Dataset(xbrdf_rli_file, 'w', format='NETCDF4')

        da.createDimension('Npixels_orbit', Npixels)

        da.createVariable('Lisparse', 'f4', ('Npixels_orbit'))
        da.variables['Lisparse'][:] = Xbrdf_rli_orbit[0, :]

        da.createVariable('Rossthick', 'f4', ('Npixels_orbit'))
        da.variables['Rossthick'][:] = Xbrdf_rli_orbit[1, :]

        da.createVariable('lon', 'f4', ('Npixels_orbit'))
        da.createVariable('lat', 'f4', ('Npixels_orbit'))
        da.variables['lat'][:] = lat_orbit
        da.variables['lon'][:] = lon_orbit

        da.close()

    def write_to_group(self, group, start=0, end=None):

        if not self.combined_with_sentinel:

            group.variables['Bands_combine'][:] = self.waves_combined_bands
            group.variables['Xbrdf_wave_combine'][:, start:end] = (
                self.xbrdf_wave_combine_orbit[:])

        group.variables['Lisparse'][start:end] = self.xbrdf_rli_orbit[0, :]
        group.variables['Rossthick'][start:end] = self.xbrdf_rli_orbit[1, :]
        group.variables['water_frac_inland'][start:end] = (
            0.0 * np.ones(self.n_pixels))
        group.variables['water_wind'][start:end] = 7.0 * np.ones(self.n_pixels)
        group.variables['snow_kernel_coef'][start:end] = (
            0.0 * np.ones(self.n_pixels))
