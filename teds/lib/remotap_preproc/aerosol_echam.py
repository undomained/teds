from netCDF4 import Dataset
import numpy as np
import os
from .exceptions import ProcessError
from .tool_algorithm import convert2julian7
from .collocation_algorithm import interpolate_to_orbit

import logging
_logger = logging.getLogger(__name__)


# Very different from the aerosol echam in MiPrep
class AerosolEcham(object):

    def __init__(self, path, year, month, day, refr_mode_coeffcients):
        self.source = ""
        self.description = ""
        self.path = path
        self.year = year
        self.month = month
        self.day = day
        self.refr_mode_coeffcients = refr_mode_coeffcients
        self.fname_echam_old = self.path + 'ECHAM_aerosol_new_scaled.nc'
        self.fname_echam = self.path + 'ECHAM_aerosol_new_scaled_DEMheight.nc'

        self.varname_echam = [
            'ALP_AI', 'ALP_AS', 'ALP_CI', 'ALP_CS', 'ALP_KI', 'ALP_KS', 'ALP_NS',
                                'num_AI', 'num_AS', 'num_CI', 'num_CS', 'num_KI', 'num_KS', 'num_NS',
                                'reff_AI', 'reff_AS', 'reff_CI', 'reff_CS', 'reff_KI', 'reff_KS', 'reff_NS',
                                'vfrac_BC_AI', 'vfrac_BC_AS', 'vfrac_BC_CI', 'vfrac_BC_CS', 'vfrac_BC_KI', 'vfrac_BC_KS', 'vfrac_BC_NS',
                                'vfrac_DU_AI', 'vfrac_DU_AS', 'vfrac_DU_CI', 'vfrac_DU_CS', 'vfrac_DU_KI', 'vfrac_DU_KS', 'vfrac_DU_NS',
                                'vfrac_H2O_AI', 'vfrac_H2O_AS', 'vfrac_H2O_CI', 'vfrac_H2O_CS', 'vfrac_H2O_KI', 'vfrac_H2O_KS', 'vfrac_H2O_NS',
                                'vfrac_OC_AI', 'vfrac_OC_AS', 'vfrac_OC_CI', 'vfrac_OC_CS', 'vfrac_OC_KI', 'vfrac_OC_KS', 'vfrac_OC_NS',
                                'vfrac_SO4_AI', 'vfrac_SO4_AS', 'vfrac_SO4_CI', 'vfrac_SO4_CS', 'vfrac_SO4_KI', 'vfrac_SO4_KS', 'vfrac_SO4_NS',
                                'vfrac_SS_AI', 'vfrac_SS_AS', 'vfrac_SS_CI', 'vfrac_SS_CS', 'vfrac_SS_KI', 'vfrac_SS_KS', 'vfrac_SS_NS',
                                'ALH_AI', 'ALH_AS', 'ALH_CI', 'ALH_CS', 'ALH_KI', 'ALH_KS', 'ALH_NS',
                                'pressure_surface_echam', 'elevation_surface_echam'
        ]

        # cham data: 1.875(192 lon)*?(96 lat,degree)*1(365 days) degree resolution data
        self.varname_echam_part3 = ['veff_AI', 'veff_AS', 'veff_CI', 'veff_CS', 'veff_KI', 'veff_KS', 'veff_NS']
        self.varname_echam_part4 = ['mr_550nm_AI', 'mr_550nm_AS', 'mr_550nm_CI', 'mr_550nm_CS', 'mr_550nm_KI', 'mr_550nm_KS', 'mr_550nm_NS']
        self.ivar_echam_part4_start = 21
        self.varname_echam_part5 = ['mi_550nm_AI', 'mi_550nm_AS', 'mi_550nm_CI', 'mi_550nm_CS', 'mi_550nm_KI', 'mi_550nm_KS',
                                    'mi_550nm_NS']

        self.ivar_echam_part5_start = 21

        self.fname_echam_prefix = 'pars_echam'

        self.lats_echam = None
        self.lons_echam = None
        self.echam_full_in = None

        self.lat_orbit = None
        self.lon_orbit = None
        self.n_pixels = None
        self.parameter_echam_orbit = None

        self.parameter_echam_part3_orbit = None
        self.parameter_echam_part4_orbit = None
        self.parameter_echam_part5_orbit = None

    def read_file(self):

        _logger.info("Reading echam file {}".format(self.fname_echam))

        iday_echam = int(convert2julian7(int(self.year), int(self.month), int(self.day), 0., 0., 0., 0.))
        
        if not os.path.isfile(self.fname_echam):
            _logger.warning("Newer ECHAM file {} does not exist. Defaulting to {}".format(self.fname_echam, self. fname_echam_old))
            self.fname_echam = self.fname_echam_old

        _logger.info("Reading echam file {}".format(self.fname_echam))


        root_grp = Dataset(self.fname_echam)

        self.lats_echam = root_grp.variables['lat'][:]   # lats 1d (96,)
        self.lons_echam = root_grp.variables['lon'][:]   # lons 1d (192,)
        times_echam = root_grp.variables['time'][:]   # times 1d (365,)

        _logger.debug('echam lats_echam.shape = {}'.format(self.lats_echam.shape))
        _logger.debug('echam lons_echam.shape = {}'.format(self.lons_echam.shape))
        _logger.debug('echam times_echam.shape = {}'.format(times_echam.shape))

        nlon_echam = len(self.lons_echam)
        nlat_echam = len(self.lats_echam)
        nvariables_echam = len(self.varname_echam)

        self.echam_full_in = np.zeros((nvariables_echam, nlat_echam, nlon_echam))

        for ivar_echam in range(0, nvariables_echam):
            _logger.debug('Read in variables {}'.format(self.varname_echam[ivar_echam]))

            self.echam_full_in[ivar_echam, :, :] = root_grp.variables[self.varname_echam[ivar_echam]][iday_echam][:][:]  # ,reff_CS, 96(lat)*192(lon)

        root_grp.close()

    def collocate(self, julday_orbit, lat, lon, n_pixels):

        _logger.info("Collocation of aerosol_echam file started.")
        # Check mode_coeff format.
        self.collocate_part1(lat, lon, n_pixels)
        self.collocate_part3(lat, lon, n_pixels)
        self.collocate_part4(n_pixels, self.refr_mode_coeffcients)
        self.collocate_part5(n_pixels, self.refr_mode_coeffcients)
        _logger.info("Collocation of aerosoal_echam file done.")

    def collocate_part1(self, lat, lon, n_pixels):

        self.lat_orbit = lat
        self.lon_orbit = lon
        self.n_pixels = n_pixels

        if self.echam_full_in is None:
            self.read_file()

        self.parameter_echam_orbit = interpolate_to_orbit(self.lat_orbit, self.lon_orbit, self.n_pixels,
                                                          self.lats_echam, self.lons_echam, self.echam_full_in,
                                                          len(self.varname_echam))

    def collocate_part3(self, lat, lon, n_pixels):

        _logger.debug('now for echam_part3 aerosol prameters (mr) in orbit')

        if self.echam_full_in is None:
            self.read_file()

        nvariables_echam_part3 = len(self.varname_echam_part3)

        self.parameter_echam_part3_orbit = np.zeros((nvariables_echam_part3, n_pixels))

        for ivar_echam_part3 in range(0, nvariables_echam_part3):  # (0,Nvariables_echam_part3)

            if (self.varname_echam_part3[ivar_echam_part3] == 'veff_CI') or (self.varname_echam_part3[ivar_echam_part3] == 'veff_CS'):
                self.parameter_echam_part3_orbit[ivar_echam_part3, :] = 0.62
            else:
                self.parameter_echam_part3_orbit[ivar_echam_part3, :] = 0.24

    def collocate_part4(self, n_pixels, mode_coeff):

        if self.echam_full_in is None:
            self.read_file()

        _logger.debug('now for echam_part4 aerosol prameters (mr) in orbit')

        # Check type of mode coefficients.

        # interpolation to orbit data
        nvariables_echam_part4 = len(self.varname_echam_part4)
        self.parameter_echam_part4_orbit = np.zeros((nvariables_echam_part4, n_pixels))

        for ivar_echam_part4 in range(0, nvariables_echam_part4):  # (0, Nvariables_echam_part4)

            self.parameter_echam_part4_orbit[ivar_echam_part4, :] = \
                self.parameter_echam_orbit[self.ivar_echam_part4_start + nvariables_echam_part4 * 0 + ivar_echam_part4, :] * mode_coeff.mr_BC + \
                self.parameter_echam_orbit[self.ivar_echam_part4_start + nvariables_echam_part4 * 1 + ivar_echam_part4, :] * mode_coeff.mr_DUST + \
                self.parameter_echam_orbit[self.ivar_echam_part4_start + nvariables_echam_part4 * 2 + ivar_echam_part4, :] * mode_coeff.mr_H2O + \
                self.parameter_echam_orbit[self.ivar_echam_part4_start + nvariables_echam_part4 * 3 + ivar_echam_part4, :] * mode_coeff.mr_OC + \
                self.parameter_echam_orbit[self.ivar_echam_part4_start + nvariables_echam_part4 * 4 + ivar_echam_part4, :] * mode_coeff.mr_IO + \
                self.parameter_echam_orbit[self.ivar_echam_part4_start + nvariables_echam_part4 * 5 + ivar_echam_part4, :] * mode_coeff.mr_IO

    def collocate_part5(self, n_pixels, mode_coeff):

        _logger.debug('now for echam_part5 aerosol prameters (mr) in orbit')

        # Check format of mode_coeff.

        if self.echam_full_in is None:
            self.read_file()

        nvariables_echam_part5 = len(self.varname_echam_part5)

        # interpolation to orbit data
        self.parameter_echam_part5_orbit = np.zeros((nvariables_echam_part5, n_pixels))

        for ivar_echam_part5 in range(0, nvariables_echam_part5):  # (0, Nvariables_echam_part5)

            self.parameter_echam_part5_orbit[ivar_echam_part5, :] = \
                self.parameter_echam_orbit[self.ivar_echam_part5_start + nvariables_echam_part5 * 0 + ivar_echam_part5, :] * mode_coeff.mi_BC + \
                self.parameter_echam_orbit[self.ivar_echam_part5_start + nvariables_echam_part5 * 1 + ivar_echam_part5, :] * mode_coeff.mi_DUST + \
                self.parameter_echam_orbit[self.ivar_echam_part5_start + nvariables_echam_part5 * 2 + ivar_echam_part5, :] * mode_coeff.mi_H2O + \
                self.parameter_echam_orbit[self.ivar_echam_part5_start + nvariables_echam_part5 * 3 + ivar_echam_part5, :] * mode_coeff.mi_OC + \
                self.parameter_echam_orbit[self.ivar_echam_part5_start + nvariables_echam_part5 * 4 + ivar_echam_part5, :] * mode_coeff.mi_IO + \
                self.parameter_echam_orbit[self.ivar_echam_part5_start + nvariables_echam_part5 * 5 + ivar_echam_part5, :] * mode_coeff.mi_IO

    def write_output(self, output_dir):

        # Check if all data is available etc....
        echam_file = output_dir + self.fname_echam_prefix+'.nc'
        _logger.info("Writing echam aerosols to {}".format(echam_file))

        nvar_echam = len(self.varname_echam)
        nvar_echam_part3 = len(self.varname_echam_part3)
        nvar_echam_part4 = len(self.varname_echam_part4)
        nvar_echam_part5 = len(self.varname_echam_part5)

        da = Dataset(echam_file, 'w', format='NETCDF4')
        da.createDimension('Npixels_orbit', self.n_pixels)

        da.createVariable('lon', 'f4', ('Npixels_orbit'))
        da.createVariable('lat', 'f4', ('Npixels_orbit'))
        da.variables['lat'][:] = self.lat_orbit
        da.variables['lon'][:] = self.lon_orbit

        for ivarname_echam in range(0, nvar_echam):
            _logger.debug("Writing echam variables")
            variablename = self.varname_echam[ivarname_echam]
            da.createVariable(variablename, 'f4', ('Npixels_orbit'))
            da.variables[variablename][:] = self.parameter_echam_orbit[ivarname_echam, :]

        for ivarname_echam_part3 in range(0, nvar_echam_part3):
            _logger.debug("Writing echam 'part 3' variables")
            variablename = self.varname_echam_part3[ivarname_echam_part3]
            da.createVariable(variablename, 'f4', ('Npixels_orbit'))
            da.variables[variablename][:] = self.parameter_echam_part3_orbit[ivarname_echam_part3, :]

        for ivarname_echam_part4 in range(0, nvar_echam_part4):
            _logger.debug("Writing echam 'part 4' variables")
            variablename = self.varname_echam_part4[ivarname_echam_part4]
            da.createVariable(variablename, 'f4', ('Npixels_orbit'))
            da.variables[variablename][:] = self.parameter_echam_part4_orbit[ivarname_echam_part4, :]

        for ivarname_echam_part5 in range(0, nvar_echam_part5):
            _logger.debug("Writing echam 'part 5' variables")
            variablename = self.varname_echam_part5[ivarname_echam_part5]
            da.createVariable(variablename, 'f4', ('Npixels_orbit'))
            da.variables[variablename][:] = self.parameter_echam_part5_orbit[ivarname_echam_part5, :]

        da.close()

    def get_orbit_parameter(self, parameter_name, julday_orbit, lat, lon, n_pixels):

        if self.parameter_echam_orbit is None:
            self.collocate(julday_orbit, lat, lon, n_pixels)

        if parameter_name in self.varname_echam:
            return self.parameter_echam_orbit[self.varname_echam.index(parameter_name)]
        if parameter_name in self.varname_echam_part3:
            return self.parameter_echam_part3_orbit[self.varname_echam_part3.index(parameter_name)]
        if parameter_name in self.varname_echam_part4:
            return self.parameter_echam_part4_orbit[self.varname_echam_part4.index(parameter_name)]
        if parameter_name in self.varname_echam_part5:
            return self.parameter_echam_part5_orbit[self.varname_echam_part5.index(parameter_name)]

        raise ProcessError("Parameter {} not found in {}".format(parameter_name, self.fname_echam))

    def scale_to_grid(self, coarse_grid, fine_grid):

        _logger.info("Scaling the 'num_AI','num_AS','num_CI','num_CS','num_KI','num_KS','num_NS' parameters to the fine and coarse grid.")
        scale_f = fine_grid / (self.parameter_echam_orbit[7, :] + self.parameter_echam_orbit[8, :])  # AI, AS
        self.parameter_echam_orbit[7, :] = scale_f * self.parameter_echam_orbit[7, :]
        self.parameter_echam_orbit[8, :] = scale_f * self.parameter_echam_orbit[8, :]
        scale_c = coarse_grid / (self.parameter_echam_orbit[9, :] + self.parameter_echam_orbit[10, :])  # CI, CS
        self.parameter_echam_orbit[9, :] = scale_c * self.parameter_echam_orbit[9, :]
        self.parameter_echam_orbit[10, :] = scale_c * self.parameter_echam_orbit[10, :]
        

    def write_to_group(self, group, start=0, end=None):

        varname_aerosol_out = [
            'sphere_mode4', 'sphere_mode5', 'sphere_mode6', 'sphere_mode7', 'sphere_mode1', 'sphere_mode2',
            'sphere_mode3',
            'num_mode4', 'num_mode5', 'num_mode6', 'num_mode7', 'num_mode1', 'num_mode2', 'num_mode3',
            'reff_mode4', 'reff_mode5', 'reff_mode6', 'reff_mode7', 'reff_mode1', 'reff_mode2', 'reff_mode3',
            'vfrac_species1_mode4', 'vfrac_species1_mode5', 'vfrac_species1_mode6', 'vfrac_species1_mode7',
            'vfrac_species1_mode1', 'vfrac_species1_mode2', 'vfrac_species1_mode3',
            'vfrac_species2_mode4', 'vfrac_species2_mode5', 'vfrac_species2_mode6', 'vfrac_species2_mode7',
            'vfrac_species2_mode1', 'vfrac_species2_mode2', 'vfrac_species2_mode3',
            'vfrac_species3_mode4', 'vfrac_species3_mode5', 'vfrac_species3_mode6', 'vfrac_species3_mode7',
            'vfrac_species3_mode1', 'vfrac_species3_mode2', 'vfrac_species3_mode3',
            'vfrac_species4_mode4', 'vfrac_species4_mode5', 'vfrac_species4_mode6', 'vfrac_species4_mode7',
            'vfrac_species4_mode1', 'vfrac_species4_mode2', 'vfrac_species4_mode3',
            'vfrac_species5_mode4', 'vfrac_species5_mode5', 'vfrac_species5_mode6', 'vfrac_species5_mode7',
            'vfrac_species5_mode1', 'vfrac_species5_mode2', 'vfrac_species5_mode3',
            'vfrac_species6_mode4', 'vfrac_species6_mode5', 'vfrac_species6_mode6', 'vfrac_species6_mode7',
            'vfrac_species6_mode1', 'vfrac_species6_mode2', 'vfrac_species6_mode3',
            'ALH_mode4', 'ALH_mode5', 'ALH_mode6', 'ALH_mode7', 'ALH_mode1', 'ALH_mode2', 'ALH_mode3',
            'veff_mode4', 'veff_mode5', 'veff_mode6', 'veff_mode7', 'veff_mode1', 'veff_mode2', 'veff_mode3',
            'pressure_surface_aerosol', 'elevation_surface_aerosol']

        # get spherical fraction.
        self.parameter_echam_orbit[0:7, :] = 1.0 - self.parameter_echam_orbit[28:35, :]

        for i, var in enumerate(varname_aerosol_out):
            if i >= len(self.varname_echam) - 2:
                if i < len(self.varname_echam) + len(self.varname_echam_part3) - 2:
                    group.variables[var][start:end] = self.parameter_echam_part3_orbit[i-(len(self.varname_echam) - 2)][:]
                else:
                    group.variables[var][start:end] = self.parameter_echam_orbit[i-len(self.varname_echam_part3)][:]
            else:
                group.variables[var][start:end] = self.parameter_echam_orbit[i][:]
