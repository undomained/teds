import logging
import random
import numpy as np
from netCDF4 import Dataset
from .tool_algorithm import determine_season_ids
from .collocation_algorithm import interpolate_to_orbit
from .write_module import write_parameters_to_nc

_logger = logging.getLogger(__name__)


class AerosolPolder(object):

    def __init__(self, path, year, month, refr_mode_coeff=None):

        self.source = ""
        self.description = ""
        self.path = path
        self.year = year
        self.month = month
        self.refr_mode_coeff = refr_mode_coeff

        #year_polder = '2006'
        self.varname = [
            'AOT_440nm', 'AOT_490nm', 'AOT_563nm', 'AOT_670nm', 'AOT_865nm', 'AOT_1020nm',
            'SSA_440nm', 'SSA_490nm', 'SSA_563nm', 'SSA_670nm', 'SSA_865nm', 'SSA_1020nm',
            'N_fine', 'N_coarse', 'reff_fine', 'veff_fine', 'reff_coarse', 'veff_coarse',
            'm_r_fine', 'm_r_coarse', 'm_i_fine', 'm_i_coarse', 'sphere_frac_fine', 'sphere_frac_coarse',
            'number_of_points_land', 'number_of_points_ocean', 'psurf',

        ]
        self.var_coeff = ['Coef_fine_Inorg', 'Coef_fine_BC', 'Coef_coarse_Inorg', 'Coef_coarse_Dust']
        self.fname_aero_prefix = 'aerosol_pars'
        self.fname_coeff_prefix = 'aerosol_coeff'

        self.lats = None
        self.lons = None

        self.dataset_full = None
        self.polderdata_seasonal_mean_in = None

        self.orbit_polderdata = None
        self.orbit_coeff = None

    def read_orbit_data(self, filepath):

        root_grp = Dataset(filepath)

        self.lats = root_grp.variables['lat_center'][:]  # lats 1d (180,)
        self.lons = root_grp.variables['lon_center'][:]  # lons 1d (360,)

        AOT_440 = root_grp.variables['AOT_440nm'][:][:][:]  # AOT_440:   30(days)*360(lon)*180(lat)
        ndays, nlon, nlat = AOT_440.shape

        root_grp.close()

        return nlon, nlat

    def read_file(self):

        case_fnames = determine_season_ids(self.month)

        fname_orbit_prefix = self.path + 'L1C3/GLOBAL/GRIDDED14/'+self.year+'/SRON_parasol_gridded'+self.year
        fname_orbit_suffix = '.nc'

        _logger.info("Reading aeosol polder file {}".format(self.fname_orbit_prefix))

        fname_orbit_tmp = fname_orbit_prefix + case_fnames[0] + fname_orbit_suffix
        Nlon, Nlat = self.read_orbit_data(fname_orbit_tmp)
        icount = 0

        nvariables = len(self.varname)
        self.dataset_full = np.zeros((nvariables, Nlon, Nlat, len(case_fnames)))

        #########################
        for case_fname in case_fnames:

            fname_orbit_tmp = fname_orbit_prefix + case_fname + fname_orbit_suffix

            #########################
            # Polder data: 1*1 degree resolution data
            root_grp = Dataset(fname_orbit_tmp)

            for ivar in range(0, nvariables):
                self.dataset_full[ivar, :, :, icount] = np.nanmean(
                    root_grp.variables[self.varname[ivar]][:][:][:], axis=0)  # e.g., AOT_440:   30(days)*360(lon)*180(lat)

            root_grp.close()

            self.polderdata_seasonal_mean_in = np.nanmean(self.dataset_full, axis=3)  # Polderdata_seasonal_mean:   28(variables)*360(lon)*180(lat)

    def collocate(self, julday_orbit, lat, lon, n_pixels):

        _logger.info("Collocation of aerosol polder {}".format(self.fname_orbit_prefix))

        if self.polderdata_seasonal_mean_in is None:
            self.read_file()

        self.collocate_data(lat, lon, n_pixels)
        if self.refr_mode_coeff is not None:
            self.add_mode_coeffs()

    def collocate_data(self, lat, lon, n_pixels):

        self.lat_orbit = lat
        self.lon_orbit = lon
        self.n_pixels = n_pixels

        self.orbit_polderdata = interpolate_to_orbit(self.lat_orbit, self.lon_orbit, n_pixels, self.lats, self.lons,
                                                     self.polderdata_seasonal_mean_in, len(self.varname), transpose=True)

    def add_mode_coeffs(self):

        self.orbit_coeff = np.zeros((len(self.var_coeff), self.n_pixels))

        for ipixel in range(0, self.n_pixels):

            ###############################################
            random.seed(ipixel)

            maxtry = 50
            step_fac = 0.9

            #############
            orbit_mr_f = self.orbit_polderdata[18, ipixel]
            orbit_mi_f = self.orbit_polderdata[20, ipixel]

            orbit_mr_c = self.orbit_polderdata[19, ipixel]
            orbit_mi_c = self.orbit_polderdata[21, ipixel]

            ###############################################

            a = np.array([[self.refr_mode_coeff.mr_IO, self.refr_mode_coeff.mr_BC], [self.refr_mode_coeff.mi_IO, self.refr_mode_coeff.mi_BC]])
            b = np.array([orbit_mr_f, orbit_mi_f])

            coef_fine = np.linalg.solve(a, b)

            coef_fine[0] = min(max(coef_fine[0], 1e-5), 1.0)
            coef_fine[1] = min(max(coef_fine[1], 1e-5), 1.0)

            #############

            c = np.array([[self.refr_mode_coeff.mr_IO, self.refr_mode_coeff.mr_DUST], [self.refr_mode_coeff.mi_IO, self.refr_mode_coeff.mi_DUST]])
            d = np.array([orbit_mr_c, orbit_mi_c])

            for itry in range(0, maxtry):
                Coef_coarse = np.linalg.solve(c, d)
                if (1. >= Coef_coarse[0] >= 0.) and (1. >= Coef_coarse[1] >= 0.):
                    break
                d[1] = d[1] * step_fac

            Coef_coarse[0] = min(max(Coef_coarse[0], 1e-5), 1.0)
            Coef_coarse[1] = min(max(Coef_coarse[1], 1e-5), 1.0)

            ###############################################
            ###############################################
            ###############################################
            self.orbit_coeff[0, ipixel] = coef_fine[0]  # fraction of Inorg for fine mode
            self.orbit_coeff[1, ipixel] = coef_fine[1]  # fraction of black carbon for fine mode
            self.orbit_coeff[2, ipixel] = Coef_coarse[0]  # fraction of Inorg for coarse mode
            self.orbit_coeff[3, ipixel] = Coef_coarse[1]  # fraction of dust for coarse mode

            # parameter_orbit[Nvariables, ipixel] = coef_fine[0]  # fraction of Inorg for fine mode
            # parameter_orbit[Nvariables + 1, ipixel] = coef_fine[1]  # fraction of black carbon for fine mode
            # parameter_orbit[Nvariables + 2, ipixel] = Coef_coarse[0]  # fraction of Inorg for coarse mode
            # parameter_orbit[Nvariables + 3, ipixel] = Coef_coarse[1]  # fraction of dust for coarse mode

    def get_input_parameter(self, parameter):
        if self.polderdata_seasonal_mean_in is None:
            self.read_file()

            ipar = self.varname.index(parameter)
            return self.polderdata_seasonal_mean_in[ipar, :]

    def get_orbit_parameter(self, parameter, julday_orbit, lat, lon, npixels):
        if self.orbit_polderdata is None:
            self.collocate(julday_orbit, lat, lon, npixels)

        ipar = self.varname.index(parameter)
        return self.orbit_polderdata[ipar, :]

    def write_output(self, output, alt_orbit, julday):
        # write polder and add the coefficients.
        main_file = output + self.fname_aero_prefix+'.nc'

        _logger.info("Writing polder output to: {}".format(main_file))

        write_parameters_to_nc(self.orbit_polderdata, self.varname, self.lat_orbit, self.lon_orbit, alt_orbit,
                               julday, self.n_pixels, main_file)

        coeff_file = output + self.fname_coeff_prefix+'.nc'

        _logger.info("Writing coefficients to: {}".format(coeff_file))

        write_parameters_to_nc(self.orbit_coeff, self.var_coeff, self.lat_orbit, self.lon_orbit, alt_orbit,
                               julday, self.n_pixels, coeff_file)

    def write_to_group(self, group):
        pass
