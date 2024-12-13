import sys
import numpy as np
from netCDF4 import Dataset
import scipy.interpolate as interp

from .exceptions import ProcessError
from .collocation_algorithm import interpolate_single_parameter_to_orbit, interpolate_to_orbit


import logging
_logger = logging.getLogger(__name__)


##################################################################################
def julday2hours(julday, year):
    hour_base = 0
    for iyear in range(1900, year):
        if (iyear % 4 == 0):
            hour_base = hour_base + 366
        else:
            hour_base = hour_base + 365

    hour_base = hour_base * 24
    hours = hour_base + julday * 24.0

    return hours


##################################################################################

##################################################################################
##################################################################################
grav = 9.80665                  # Standard gravity at ground [m/s2]
NAvog = 6.02214076E23            # Avogadro number [#/mol]
Mdry_air = 28.9647E-3           # Molar mass dry air [kg/mol]
avoga = NAvog
air_m = Mdry_air
sp = NAvog/(Mdry_air*grav)*1.E2  # [#/m^2 * 1/hPa] air column above P is P*NA/Mair/g from p = m*g/area

mol_mass_co2 = 44.01  # g/mol
mol_mass_ch4 = 16.04  # g/mol
mol_mass_h2o = 18.01528  # g/mol
mol_mass_no2 = 46.0055  # g/mol
mol_mass_o3 = 48.0  # g/mol
mol_mass_dry = 28.97  # g/mol

year = 2017  # must be 2017


hybrid_a = np.array([0.0e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01,
                     3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
                     7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, 1.091817e+02, 1.189586e+02,
                     1.286959e+02, 1.429100e+02, 1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
                     2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02,
                     2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
                     7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01, 1.979160e+01, 9.292942e+00,
                     4.076571e+00, 1.650790e+00, 6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02,
                     0.000000e+00])  # [hpa]

hybrid_a = np.flip(hybrid_a)  # now from bottom to top

hybrid_b = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
                     8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
                     7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01,
                     5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
                     2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
                     6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                     0.000000e+00])  # [-]

hybrid_b = np.flip(hybrid_b)  # now from bottom to top


def get_number_of_levels():
    return len(hybrid_a)
#####################################################################################


class WRFCO2M(object):
    def __init__(self, path):

        self.path = path
        self.filename = path + 'wrfout_d03_20130721.nc4'
        self.julday = 181 + 21
        self.year = 2013
        # use AFGL profile to calculate corresponding altitude leveles
        self.afgl_file = path + 'prof.AFGL.US.std'
        self.xco2_in = None
        self.xco2 = None
        self.co2_final = None
        self.mtime = None  # Time selection.

        self.lat_in = None
        self.lon_in = None

        self.lat = None
        self.lon = None

        self.nlat = None
        self.nlon = None
        self.nlev = None
        self.press_final = None
        self.press_final_orbit = None

    def read_file(self):

        data = Dataset(self.filename)

        # xtarget = 2500   # UNUSED
        # ytarget = 4000

        nlev = data["co2"][:, 0, 0].size
        nlat = data["co2"][0, :, 0].size
        nlon = data["co2"][0, 0, :].size

        co2 = data["co2"]    # co2 dry air mixing ratio (nlime, nlev, nlat, nlon)
        spress = data["ps"]  # surface pressure Pa ->hPa (ntime, nlat, nlon)
        self.lat_in = data["lat"]    # latitude, northern degree (ntime, nlat, nlon)
        self.lon_in = data["lon"]    # longitude, degree east (ntime, nlat, nlon)
        eta = data["eta"]    # fraction of surface pressure, with ptop at 50hPa (ntime, nlev+1)

        # ntime = co2[:, 0, 0, 0].size
        nlev = co2[0, :, 0, 0].size

        self.lat = co2[0, 0, :, 0]
        self.lon = co2[0, 0, 0, :]
        self.nlat = co2[0, 0, :, 0].size
        self.nlon = co2[0, 0, 0, :].size

        # reverse level/layer index (bottom-> top)
        eta = np.flip(eta, 1)
        co2 = np.flip(co2, 1)

        # select time 12 h mtime = 12 (?)
        self.mtime = 12
        press = np.empty([nlev+1, nlat, nlon])

        for mlat in range(self.nlat):
            for mlon in range(self.nlon):
                press[:, mlat, mlon] = spress[self.mtime, mlat, mlon] * eta[self.mtime, :]

        # From Pa to HPa
        press = press * 1.E-2
        press[0, :, :] = 50  # see decription of upper pressure boundary in netcdf file

        # Read AFGL filepri
        atm_afgl = np.genfromtxt(self.afgl_file, skip_header=2)

        press_afgl = atm_afgl[:, 1]           # pressure [hPa]
        air_afgl = atm_afgl[:, 3]              # air number density [#/cm^3]
        co2_afgl = atm_afgl[:, 7]/air_afgl     # co2 number density -> mole fraction [-]

        # To expand the wrf profile, we define two indices index1 and index2. index1 indicates
        # the layer where afgl and wrf profiles are clued together. Between index1 and index2
        # we scale the afgl to have a smooth transition at level index1. Above index2, the scaling
        # depends on layer index such that at TOA we obtain the original afgl value.
        # Here, index2 is chosen such that afgl mixing ratio is nearly constant at altitudes below.

        ind = (np.abs(press_afgl - 50.)).argmin()
        if (press_afgl[ind] <= 50):
            index1 = ind
        else:
            index1 = ind-1

        index2 = 10

        self.press_final = np.empty([nlev+1+index1, self.nlat, self.nlon])
        self.co2_final = np.empty([nlev+index1, self.nlat, self.nlon])
        scale = np.empty([nlev+index1])
        self.xco2_in = np.empty([self.nlat, self.nlon])

        for mlat in range(self.nlat):
            for mlon in range(self.nlon):

                self.press_final[:, mlat, mlon] = np.append(press_afgl[:index1], press[:, mlat, mlon])
                self.co2_final[:, mlat, mlon] = np.append(co2_afgl[:index1]*1.E6, co2[self.mtime, :, mlat, mlon])

                # define a scaling to get a smooth transition between afgl and wrf
                fact = (co2[self.mtime, 0, mlat, mlon] / co2_afgl[index1+1]*1.E-6)-1

                scale[index1:] = 1.              # range of the wrf profile
                scale[index2:index1] = 1.+fact  # constant scaling of afgl profile for smooth extension of wrf
                scale[:index2] = 1. + fact*np.arange(index2)/float(index2)  # to get upper value of afgl at TOA

                self.co2_final[:, mlat, mlon] = scale * self.co2_final[:, mlat, mlon]

                # calculate dry air column mixing ratio

                vc_air = sp * self.press_final[:, mlat, mlon]  # air column [#/m^2] above pressure level

                self.xco2_in[mlat, mlon] = np.sum((vc_air[1:] - vc_air[0:vc_air.size-1]) * self.co2_final[:, mlat, mlon]) / vc_air[vc_air.size-1]

        self.nlev = nlev + index1

    def collocate(self, julday_orbit, lat_orbit, lon_orbit, n_pixels):

        if self.lat is None:
            _logger.debug("Reading WRF data")
            self.read_file()

        _logger.debug("Collocation of wrf co2.")
        self.lat = lat_orbit
        self.lon = lon_orbit
        self.n_pixels = n_pixels
        self.xco2 = interpolate_single_parameter_to_orbit(lat_orbit, lon_orbit, n_pixels,
                                                          self.lat_in[self.mtime, :, :], self.lon_in[self.mtime, :, :],
                                                          self.xco2_in, lonlat_gridded=True, method='linear')

        self.press_final_orbit = interpolate_to_orbit(lat_orbit, lon_orbit, n_pixels, self.lat_in[self.mtime, :, :], self.lon_in[self.mtime, :, :],
                                                      self.press_final, self.nlev+1, lonlat_gridded=True, method='linear')

        self.co2_final_orbit = interpolate_to_orbit(lat_orbit, lon_orbit, n_pixels, self.lat_in[self.mtime, :, :], self.lon_in[self.mtime, :, :],
                                                    self.co2_final, self.nlev, lonlat_gridded=True, method='linear')

    def write_output(self, output_dir):
        # Write data to file.

        file_name = output_dir + '/wrf_data.nc'
        da = Dataset(file_name, 'w', format='NETCDF4')

        da.createDimension('Npixels_orbit', self.n_pixels)
        da.createDimension('nlevel', self.nlev)
        da.createDimension('nlevel_p', self.nlev + 1)

        da.createVariable('lon', 'f4', ('Npixels_orbit'))
        da.createVariable('lat', 'f4', ('Npixels_orbit'))
        da.createVariable('pressure', 'f4', ('nlevel_p', 'Npixels_orbit'))
        da.createVariable('xco2', 'f4', ('Npixels_orbit'))
        da.createVariable('co2_final', 'f4', ('nlevel', 'Npixels_orbit'))

        da.variables['lat'][:] = self.lat
        da.variables['lon'][:] = self.lon
        da.variables['pressure'][:] = self.press_final_orbit
        da.variables['xco2'][:] = self.xco2
        da.variables['co2_final'][:] = self.co2_final_orbit

        da.close()

    def write_to_group(self, group, start=0, end=None):
        group['CMR_co2'][start:end] = self.xco2
        group['CMR_co2'].long_name = "dry air column mixing ratio total, top down"


class CAMSEGGCO2M(object):

    def __init__(self, path, year, month, day):
        self.path = path

        self.tcco2_cams = None
        self.tcch4_cams = None

        self.tcco2_orbit = None
        self.tcch4_orbit = None
        self.day = day
        self.month = month
        self.year = year

        self.lon_cams = None
        self.lat_cams = None

    def read_file(self, julday):
        hours = julday2hours(julday, int(self.year))
        fname_co2 = 'CAMS_EGG_3hourly.nc'
        fname_co2 = self.path + fname_co2
        f_r = Dataset(fname_co2, 'r')
        self.lon_cams = f_r.variables['longitude'][:]
        idx = np.where(self.lon_cams > 180)
        self.lon_cams[idx] = self.lon_cams[idx] - 360
        self.lat_cams = f_r.variables['latitude'][:]
        time_cams = f_r.variables['time'][:]
        time_s = time_cams[0]
        time_d = time_cams[1] - time_cams[0]
        idx_time = int(hours[int(0.5 * len(hours))] - time_s) / time_d
        self.tcco2_cams = f_r.variables['tcco2'][idx_time][:][:]  # [ppm] dry air
        self.tcch4_cams = f_r.variables['tcch4'][idx_time][:][:]  # [ppb] dry air

        f_r.close()

    def collocate(self, julday_orbit, lat, lon, npixels):

        if self.tcco2_cams is None:
            self.read_file(julday_orbit)

        # nearest_linear.
        self.tcco2_orbit = interpolate_single_parameter_to_orbit(lat, lon, npixels, self.lat_cams, self.lon_cams, self.tcco2_cams)
        self.tcch4_orbit = interpolate_single_parameter_to_orbit(lat, lon, npixels, self.lat_cams, self.lon_cams, self.tcch4_cams)


    def write_output(self, output_dir):
        pass


class CAMSEAC4CO2M(object):
    def __init__(self, path, year, month, day):

        self.path = path
        self.year = year
        self.month = month
        self.day = day

        self.lat_cams = None
        self.lon_cams = None

        self.tcno2_cams = None
        self.tco3_cams = None
        self.tch2o_cams = None

        self.tcno2_orbit = None
        self.tco3_orbit = None
        self.tch2o_orbit = None

    def _read_file(self, julday):

        _logger.info("Reading CAMS_EAC4 file.")

        hours = julday2hours(julday, int(self.year))
        fname_no2 = 'CAMS_EAC4_3hourly.nc'
        fname_no2 = self.path + fname_no2
        f_r = Dataset(fname_no2, 'r')
        self.lon_cams = f_r.variables['longitude'][:]
        idx = np.where(self.lon_cams > 180)
        self.lon_cams[idx] = self.lon_cams[idx] - 360
        self.lat_cams = f_r.variables['latitude'][:]
        time_cams = f_r.variables['time'][:]
        time_s = time_cams[0]
        time_d = time_cams[1] - time_cams[0]
        idx_time = int(hours[int(0.5 * len(hours))] - time_s) / time_d

        self.tco3_cams = f_r.variables['gtco3'][idx_time][:][:]  # [kg m-2]
        self.tcno2_cams = f_r.variables['tcno2'][idx_time][:][:]  # [kg m-2]
        self.tch2o_cams = f_r.variables['tcwv'][idx_time][:][:]  # [kg m-2]

        f_r.close()

    def collocate(self, julday_orbit, lat, lon, npixels):

        if self.tch2o_cams is None:
            self._read_file(julday_orbit)

        _logger.info("Collocation of CAMS_EAC4 file started.")

        self.tcno2_orbit = interpolate_single_parameter_to_orbit(lat, lon, npixels, self.lat_cams, self.lon_cams, self.tcno2_cams)

        self.tco3_orbit = interpolate_single_parameter_to_orbit(lat, lon, npixels, self.lat_cams, self.lon_cams, self.tco3_cams)

        self.tch2o_orbit = interpolate_single_parameter_to_orbit(lat, lon, npixels, self.lat_cams, self.lon_cams, self.tch2o_cams)

    def write_output(self):
        pass

    def write_to_group(self, start=0, end=None):
        pass


class CO2ME2EAUX(object):
    def __init__(self, path, year, month, day, include_wrf=False):
        self.path = path
        self.fname_meteo = path + 'prof.AFGL.US.std'
        self.fname_echam = path + 'pars_echam.nc'
        self.year = year
        self.month = month
        self.day = day

        self.camsEGG = CAMSEGGCO2M(path, year, month, day)
        self.camsEAC4 = CAMSEAC4CO2M(path, year, month, day)
        self.camsEGGref = CAMSEGGCO2M(path, 2013, 7, 21)
        self.include_wrf = include_wrf

        if self.include_wrf:
            self.wrf = WRFCO2M(path)  # year, month, day

        self.n_pixels = None

        self.z_in = None  # [km] => [m]
        self.p_in = None
        self.t_in = None

        self.xo3_in = None
        self.xh2o_in = None
        self.xco2_in = None
        self.xno2_in = None
        self.xch4_in = None

        self.pixel_altitude = None
        self.psurf = None
        self.pixel_altitude_meteo = None

        self.tlay = None
        self.tropopause = None
        self.xo3 = None
        self.xh2o = None
        self.xco2 = None
        self.xco2W = None
        self.xno2 = None
        self.xch4 = None

    def read_meteo_file(self):

        _logger.info("Reading of AFGL profiles file started.")

        cols = np.loadtxt(self.fname_meteo, skiprows=2, delimiter=None)

        self.z_in = np.array(cols[:, 0]) * 1.0e3  # [km] => [m]
        self.p_in = np.array(cols[:, 1])
        self.t_in = np.array(cols[:, 2])

        dvair_in = np.array(cols[:, 3])  # [cm-3]
        do3_in = np.array(cols[:, 4])  # [cm-3]
        dh2o_in = np.array(cols[:, 6])  # [cm-3]
        dco2_in = np.array(cols[:, 7])  # [cm-3]
        dno2_in = np.array(cols[:, 8])  # [cm-3]
        dch4_in = np.array(cols[:, 9]) * dvair_in * 1.0e-6  # [ppm] => [cm-3]

        dvair_in = dvair_in * 1.e5  # dz [cm] => [km], [pars km-1 cm-2]
        do3_in = do3_in * 1.e5
        dh2o_in = dh2o_in * 1.e5
        dco2_in = dco2_in * 1.e5
        dno2_in = dno2_in * 1.e5
        dch4_in = dch4_in * 1.e5

        # convert to dry VMR, afgl
        self.xo3_in = do3_in / dvair_in  # [molar / molar]
        self.xh2o_in = dh2o_in / dvair_in  # [molar / molar]
        self.xco2_in = dco2_in / dvair_in  # [molar / molar]
        self.xno2_in = dno2_in / dvair_in  # [molar / molar]
        self.xch4_in = dch4_in / dvair_in  # [molar / molar]

    def read(self):
        _logger.info("Starting reading of meteo and CAMS files")

        self.read_meteo_file()

    def collocate(self, julday_orbit, lat_orbit, lon_orbit, n_pixels):

        self.n_pixels = n_pixels

        if self.z_in is None:
            self.read()

        _logger.info("Starting collocation of CAMS EGGfile")
        self.camsEGG.collocate(julday_orbit, lat_orbit, lon_orbit, n_pixels)
        _logger.info("Starting collocation of CAMS EAC4 file")

        self.camsEAC4.collocate(julday_orbit, lat_orbit, lon_orbit, n_pixels)
        if self.include_wrf:
            _logger.info("Starting collocation of WRF file")

            self.wrf.collocate(julday_orbit, lat_orbit, lon_orbit, n_pixels)
            self.camsEGGref.collocate(self.wrf.julday * np.ones(julday_orbit.shape), lat_orbit, lon_orbit, n_pixels)

    def post_process(self, n_pixels, psurf, pixel_altitude):

        self.psurf = psurf
        self.pixel_altitude = pixel_altitude

        # pressure
        _logger.debug("Calculating pressure")
        nlev = len(hybrid_a)
        nlay = nlev - 1

        plev = np.zeros((n_pixels, nlev))
        for ilev in range(nlev):
            plev[:, ilev] = hybrid_a[ilev] + hybrid_b[ilev] * self.psurf

        play = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            play[:, ilay] = 0.5 * (plev[:, ilay] + plev[:, ilay + 1])

        self.play = play

        # pixel_altitude_meteo, surface, tropopause
        _logger.debug("Calculating pixel altitude meteo")
        f = interp.interp1d(np.log(self.p_in), self.z_in, fill_value='extrapolate')
        self.pixel_altitude_meteo = np.zeros(n_pixels)
        self.pixel_altitude_meteo = f(np.log(self.psurf))

        zlay = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            zlay[:, ilay] = f(np.log(play[:, ilay]))

        # temp
        _logger.debug("Calculating temperature")
        f = interp.interp1d(np.log(self.p_in), self.t_in, fill_value='extrapolate')
        self.tlay = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            self.tlay[:, ilay] = f(np.log(play[:, ilay]))

        _logger.debug("Calculating tropopause")
        tlev_tmp = np.zeros(nlev)
        dtemp = np.zeros(nlay)
        zlev_tmp = np.zeros(nlev)

        self.tropopause = np.zeros(n_pixels)
        for ipix in range(n_pixels):
            tlev_tmp[1:nlev] = f(np.log(plev[ipix, 1:nlev]))
            tlev_tmp[0] = 2 * tlev_tmp[1] - tlev_tmp[2]
            zlev_tmp[1:nlev - 1] = 0.5 * (zlay[ipix, 0:nlev - 2] + zlay[ipix, 1:nlev - 1])
            zlev_tmp[0] = 2 * zlev_tmp[1] - zlev_tmp[2]
            zlev_tmp[nlev - 1] = 2 * zlev_tmp[nlev - 2] - zlev_tmp[nlev - 3]
            dtemp = (tlev_tmp[1:nlev] - tlev_tmp[0:nlev - 1]) / (zlev_tmp[1:nlev] - zlev_tmp[0:nlev - 1])
            idx5k = np.argmin(np.abs(zlay[ipix, :] - 5000))
            idx20k = np.argmin(np.abs(zlay[ipix, :] - 20000))
            idx = np.argmin(np.abs(dtemp[idx5k + 1:idx20k:-1])) + idx20k
            self.tropopause[ipix] = zlay[ipix, idx]

        # interp to get dry VMR
        f = interp.interp1d(np.log(self.p_in), self.xo3_in, fill_value='extrapolate')
        xo3 = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            xo3[:, ilay] = f(np.log(play[:, ilay]))

        f = interp.interp1d(np.log(self.p_in), self.xh2o_in, fill_value='extrapolate')
        xh2o = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            xh2o[:, ilay] = f(np.log(play[:, ilay]))

        f = interp.interp1d(np.log(self.p_in), self.xco2_in, fill_value='extrapolate')
        xco2 = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            xco2[:, ilay] = f(np.log(play[:, ilay]))

        wrfxco2 = np.zeros((n_pixels, nlay))
        if self.include_wrf:
            # Note that the original wrf data is tramsposed.
            press_center = np.zeros((self.wrf.nlev, n_pixels))
            for lay in range(self.wrf.nlev):
                press_center[lay, :] = 0.5 * (self.wrf.press_final_orbit[lay, :] + self.wrf.press_final_orbit[lay + 1, :])
            for pid in range(n_pixels):
                f = interp.interp1d(np.log(press_center[:, pid]), self.wrf.co2_final_orbit[:, pid], fill_value='extrapolate')
                # Here the data will have the same shape as the rest of the data.
                for ilay in range(nlay):
                    wrfxco2[pid, ilay] = f(np.log(play[pid, ilay]))

        f = interp.interp1d(np.log(self.p_in), self.xno2_in, fill_value='extrapolate')
        xno2 = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            xno2[:, ilay] = f(np.log(play[:, ilay]))

        f = interp.interp1d(np.log(self.p_in), self.xch4_in, fill_value='extrapolate')
        xch4 = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            xch4[:, ilay] = f(np.log(play[:, ilay]))

        # get dvair dry, rough with gravity independent on lat and alt
        dvair = np.zeros((n_pixels, nlay))
        for ilay in range(nlay):
            dvair[:, ilay] = (plev[:, ilay + 1] - plev[:, ilay]) * avoga * 1.e-2 / \
                             (air_m * grav *
                              (1. + xh2o[:, ilay] / 1.60855))  # Mair/Mwater=1.60855

        _logger.debug("XO3 scaling")
        # scale vmol
        # o3
        # vmr=> dvmol [par /km /cm-2]=> dmass[g / km / cm-2] =>[kg / km / m-2]
        tmp = xo3 * dvair / avoga * mol_mass_o3 * 1e1
        scale = self.camsEAC4.tco3_orbit / np.sum(tmp, axis=1)
        for ilay in range(nlay):
            xo3[:, ilay] = xo3[:, ilay] * scale

        _logger.debug("XO3 scaling min/max {}, {}".format(np.min(scale), np.max(scale)))

        _logger.debug("H2O scaling")
        # h2o
        # vmr=> dvmol [par /km /cm-2]=> dmass[g / km / cm-2] =>[kg / km / m-2]
        tmp = xh2o * dvair / avoga * mol_mass_h2o * 1e1
        scale = self.camsEAC4.tch2o_orbit / np.sum(tmp, axis=1)
        for ilay in range(nlay):
            xh2o[:, ilay] = xh2o[:, ilay] * scale

        _logger.debug("H2O scaling min/max {}, {}".format(np.min(scale), np.max(scale)))

        _logger.debug("NO2 scaling")
        # no2
        # vmr=> dvmol [par /km /cm-2]=> dmass[g / km / cm-2] =>[kg / km / m-2]
        tmp = xno2 * dvair / avoga * mol_mass_no2 * 1e1
        scale = self.camsEAC4.tcno2_orbit / np.sum(tmp, axis=1)
        for ilay in range(nlay):
            xno2[:, ilay] = xno2[:, ilay] * scale

        _logger.debug("NO2 scaling min/max {}, {}".format(np.min(scale), np.max(scale)))

        _logger.debug("CO2 scaling")
        # co2
        tmp = xco2 * dvair
        tmp = 1e6 * np.sum(tmp, axis=1) / np.sum(dvair, axis=1)
        scale = self.camsEGG.tcco2_orbit / tmp
        for ilay in range(nlay):
            xco2[:, ilay] = xco2[:, ilay] * scale

        _logger.debug("CO2 scaling min/max {}, {}".format(np.min(scale), np.max(scale)))

        _logger.debug("WRF CO2 scaling")
        # Correct from 2013 avg to 2017 (or other selection).
        if self.include_wrf:
            # Scale the co2 values with the CAMS and mixing ratio.
            scale = self.camsEGGref.tcco2_orbit / (self.wrf.xco2 * 1e6)
            for ilay in range(nlay):
                wrfxco2[:, ilay] = wrfxco2[:, ilay] * scale

            _logger.debug("WRF CO2 scaling min/max {}, {}".format(np.min(scale), np.max(scale)))

        # ch4
        _logger.debug("CH4 scaling")
        tmp = xch4 * dvair
        tmp = 1e9 * np.sum(tmp, axis=1) / np.sum(dvair, axis=1)
        scale = self.camsEGG.tcch4_orbit / tmp
        for ilay in range(nlay):
            xch4[:, ilay] = xch4[:, ilay] * scale

        _logger.debug("CH4 scaling min/max {}, {}".format(np.min(scale), np.max(scale)))

        # get MMR, dry
        _logger.debug("Determining the dry MMR values")
        self.xo3 = xo3 * mol_mass_o3 / mol_mass_dry
        self.xno2 = xno2 * mol_mass_no2 / mol_mass_dry
        self.xco2 = xco2 * mol_mass_co2 / mol_mass_dry

        if self.include_wrf:
            self.xco2W = wrfxco2 * mol_mass_co2 / mol_mass_dry
            indices = np.argwhere(np.isnan(self.xco2W))
            for ilay in range(nlay):
                indices = np.where(~np.isnan(self.xco2W[:, ilay]))
                self.xco2[indices, ilay] = self.xco2W[indices, ilay]

        self.xch4 = xch4 * mol_mass_ch4 / mol_mass_dry
        xh2o = xh2o * mol_mass_h2o / mol_mass_dry
        self.xh2o = xh2o / (1 + xh2o)  # dry MMR to specific humidity (moist MMR)

    def write_output(self, output_dir):
        if self.psurf is None:
            raise ProcessError("Post processing step, which needs aerosol echam variables not executed.")

        self.wrf.write_output(output_dir)

        fname_w = output_dir + "aux.nc"

        nlev = len(hybrid_a)
        nlay = nlev - 1

        f_w = Dataset(fname_w, 'w', format='NETCDF4')
        f_w.createDimension('Npix', self.n_pixels)
        f_w.createDimension('nlev', nlev)
        f_w.createDimension('nlay', nlay)

        f_w.createVariable('pixel_id', 'i4', ('Npix'))
        pixel_id = np.arange(self.n_pixels) + 1
        pixel_id.astype(int)
        f_w.variables['pixel_id'][:] = pixel_id

        co2_obj = f_w.createVariable('co2', 'f4', ('Npix', 'nlay'))
        co2_obj[:] = self.xco2_int

        press = f_w.createVariable('press', 'f4', ('Npix', 'nlay'))
        press[:] = self.play

        self.write_to_group(f_w)

        f_w.close()

    def write_to_group(self, group, start=0, end=None):

        nlev = len(hybrid_a)
        nlay = nlev - 1

        group['surface_pressure_meteo'][start:end] = self.psurf
        group["surface_elevation_dem"][start:end] = self.pixel_altitude
        group['hybrid_a'][:] = hybrid_a
        group['hybrid_b'][:] = hybrid_b
        group['surface_elevation_meteo'][start:end] = self.pixel_altitude_meteo
        group['surface_elevation_std'][start:end] = np.zeros(self.n_pixels)
        group['tropopause_height'][start:end] = self.tropopause
        group['temperature'][start:end, :] = self.tlay
        group['MMR_o3'][start:end, :] = self.xo3
        group['MMR_co2'][start:end, :] = self.xco2
        group['MMR_ch4'][start:end, :] = self.xch4
        group['MMR_no2'][start:end, :] = self.xno2
        group["specific_humidity"][start:end, :] = self.xh2o
