#==============================================================================
#   Tools to generate model atmosphere
#   This source code is licensed under the 3-clause BSD license found in
#   the LICENSE file in the root directory of this project.
#==============================================================================

import os
import sys
import numpy as np
import netCDF4 as nc
from scipy.interpolate import griddata
from tqdm import tqdm
from copy import deepcopy
from .libMeteo import readmeteodata
from . import constants
from .libNumTools import convolution_2d, TransformCoords, getconvolutionparams
# import matplotlib.pyplot as plt
###########################################################


class Emptyclass:
    """Data container."""

    pass


def shrink_extend_domain(gm_x, gm_y, x, y):
    """Extend or shrink domain.
    
    Shrink or extend domain of meteo data to GM domain.

    Parameters
    ----------
    gm_x : Matrix
        latitude GM
    gm_y : Matrix
        longitude GM
    x : Matrix
        Existing meteo data latitudes
    y : Matrix
        Existing meteo data longitudes
    
    """
    # get bounds of the GM domain
    _x = np.array([gm_x[0,0], gm_x[0,-1], gm_x[-1,0], gm_x[-1,-1]])
    _y = np.array([gm_y[0,0], gm_y[0,-1], gm_y[-1,0], gm_y[-1,-1]])

    # compute dx, dy
    dx = (x[1] - x[0])
    dy = (y[1] - y[0])
    
    # find the padding bounds and min and max index for GM domain
    pad_x1, pad_x2, pad_y1, pad_y2 = 0, 0, 0, 0  # default padding is zero
    ix1, ix2, iy1, iy2 = 0, x.size, 0, y.size
    if _x.min() <= x[0]:
        pad_x1 = np.int_((x[0] - _x.min())/dx + 1)
    else:
        ix1 = np.searchsorted(x, _x.min()) - 1
    if _x.max() > x[-1]:
        pad_x2 = np.int_((_x.max() - x[-1])/dx + 1)
    else:
        ix2 = np.searchsorted(x, _x.max())
    if _y.min() <= y[0]:
        pad_y1 = np.int_((y[0] - _y.min())/dy + 1)
    else:
        iy1 = np.searchsorted(y, _y.min()) - 1
    if _y.max() >= y[-1]:
        pad_y2 = np.int_((_y.max() - y[-1])/dy + 1)
    else:
        iy2 = np.searchsorted(y, _y.max())

    # compute new x and y and then latitude and longitude
    x_left = np.arange(-pad_x1, 0)*dx + x[0]
    x_right = np.arange(pad_x2)*dx + dx + x[-1]
    x_new = np.concatenate((x_left, x[ix1:ix2], x_right))

    y_left = np.arange(-pad_y1, 0)*dy + y[0]
    y_right = np.arange(pad_y2)*dy + dy + y[-1]
    y_new = np.concatenate((y_left, y[iy1:iy2], y_right))

    return x_new, y_new, (ix1, ix2), (iy1, iy2), (pad_x1, pad_x2), (pad_y1, pad_y2)


def flip_zyx2yxz(conc):
    """Flip indices from zyx to yxz.

    Parameters
    ----------
    conc : 3D array
        Three dimensional array

    """
    sh = conc.shape
    conc1 = np.zeros((sh[1], sh[2], sh[0]))
    for i in range(sh[0]):
        conc1[:,:,i] = conc[i,:,:]
    return conc1


def get_atmosphericdata_new(gm_lat, gm_lon, meteo_settings):
    """Get meterological data to same domain as input lat, lon.

    Read meterological data and extend to given domain.

    Parameters
    ----------
    gm_lat : Array (m,n)
        Input GM Latitude
    gm_lon : Array (m,n)
        Input GM longitude
    meteo_settings : Dict
        Dict containing meteo data settings
    """
    # read data
    # data is already in molceules/m^2
    print('Getting meteo data ...')
    data = readmeteodata(meteo_settings["path_data"], meteo_settings['gases'], meteo_settings["filesuffix"])
    # shrink or extend domain according to the input

    src = data.__getattribute__(meteo_settings['gases'][0]+"_source")
    # # if the mesh is regular (see Topology section in https://www.xdmf.org/index.php/XDMF_Model_and_Format)
    # if data.gridtype == "3DCoRectMesh":
    # create a transform method 
    trans = TransformCoords(src[1:])
    # convert lat-lon of gm to x-y and get bounds

    gm_x, gm_y = trans.latlon2xymts(gm_lat, gm_lon)

    x_new, y_new, idx, idy, padx, pady = shrink_extend_domain(gm_x, gm_y, data.x, data.y)
    XX, YY = np.meshgrid(x_new, y_new)
    lat_new, lon_new = trans.xymts2latlon(XX, YY)
    # else:
    #     print("the grids are not regular")
    #     exit()

    # create a new class to have meteo data
    meteo_data = Emptyclass()
    meteo_data.__setattr__('lat', lat_new)
    meteo_data.__setattr__('lon', lon_new)
    meteo_data.__setattr__('zlay', data.z)
    meteo_data.__setattr__('zlev', data.znodes)
    meteo_data.__setattr__("dx", data.dx)
    meteo_data.__setattr__("dy", data.dy)
    meteo_data.__setattr__("dz", data.z[1] - data.z[0])
    meteo_data.__setattr__("x_new", x_new)
    meteo_data.__setattr__("y_new", y_new)
    meteo_data.__setattr__("gm_x", gm_x)
    meteo_data.__setattr__("gm_y", gm_y)
    # Add emission and source to meteo data
    for gas in meteo_settings['gases']:
        meteo_data.__setattr__(gas+"_emission_in_kgps", data.__getattribute__(gas+"_emission_in_kgps"))
        meteo_data.__setattr__(gas+"_source", data.__getattribute__(gas+"_source"))
        conc = data.__getattribute__(gas)
        conc_new = np.pad(conc[:, idy[0]:idy[1], idx[0]:idx[1]], pad_width=((0,0), pady, padx), constant_values=0)
        meteo_data.__setattr__(gas+"_raw", flip_zyx2yxz(conc_new))
    print('                     ...done')
    return meteo_data


def get_atmosphericdata(s2_lat, s2_lon,  meteo_settings, kernel_settings):
    """Get meterological data.

    Read meterological data, convolve and interpolate it givel lat-lon.

    Parameters
    ----------
    s2_lat : Array (m,n)
        Input Latitude
    s2_lon : Array (m,n)
        Input longitude
    meteo_settings : Dict
        Dict containing meteo data settings
    kernel_settings : Dict
        Kernel settings for convolution

    """
    # read data
    # data is already in molceules/m^2
    print('Getting meteo data ...')
    data = readmeteodata(meteo_settings["path_data"], meteo_settings['gases'], meteo_settings["filesuffix"])
    print('                     ...done')

    # create a new class to have meteo data
    dim_alt, dim_act = s2_lat.shape
    meteo_data = Emptyclass()
    meteo_data.__setattr__('lat', s2_lat)
    meteo_data.__setattr__('lon', s2_lon)
    meteo_data.__setattr__('zlay', data.z)
    meteo_data.__setattr__('zlev', data.znodes)
    meteo_data.__setattr__("dx", data.dx)
    meteo_data.__setattr__("dy", data.dy)
    # Add emission and source to meteo data
    for gas in meteo_settings['gases']:
        meteo_data.__setattr__(gas+"_emission_in_kgps", data.__getattribute__(gas+"_emission_in_kgps"))
        meteo_data.__setattr__(gas+"_source", data.__getattribute__(gas+"_source"))
    # convolution of co2 field and interpolation to S2 grid
    conv_settings = getconvolutionparams(kernel_settings, data.dx, data.dy)
    for gas in meteo_settings['gases']:
        concgas = data.__getattribute__(gas)
        conv_gas = np.zeros_like(concgas)
                
        for iz in range(data.z.size):
            conv_gas[iz, :, :] = convolution_2d(concgas[iz, :, :], conv_settings)
        data.__setattr__("conv_"+gas, conv_gas)

    # Interpolate data to a given lat-lon grid
    print('Interpolating data to S2 mesh...')
    # interpolate for all given gases
    for gas in meteo_settings['gases']:
        interpdata = np.zeros([dim_alt, dim_act, data.z.size])
        conv_data = data.__getattribute__("conv_"+gas)

        dxdy = np.column_stack((data.lat.ravel(), data.lon.ravel()))
        for iz in tqdm(range(data.z.size)):
            interpdata[:, :, iz] = griddata(dxdy, conv_data[iz, :,:].ravel(), (s2_lat, s2_lon), fill_value=0.0)
        meteo_data.__setattr__(gas, interpdata)
    print('                     ...done')

    return meteo_data


def get_AFGL_atm_homogenous_distribution(AFGL_path, nlay, dzlay, xco2_ref=405, xch4_ref=1800., xh2o_ref=1.E4):

    # Vertical layering assuming equidistance gridding in geometrical distance
    nlev = nlay + 1  # number of levels
    # counting from top to bottom
    zlay = (np.arange(nlay-1, -1, -1.)+0.5)*dzlay  # altitude of layer midpoint
    zlev = np.arange(nlev-1, -1, -1.)*dzlay  # altitude of layer interfaces = levels
    psurf = 101300.0  # Pa
    atm = atmosphere_data(zlay, zlev, psurf)
    atm.get_data_AFGL(AFGL_path)

    xco2 = np.sum(atm.CO2) / np.sum(atm.air) * 1.E6
    xch4 = np.sum(atm.CH4) / np.sum(atm.air) * 1.E9
    xh2o = np.sum(atm.H2O) / np.sum(atm.air) * 1.E6
    air = np.sum(atm.air)

    atm.CO2 = xco2_ref/xco2 * atm.CO2
    atm.CH4 = xch4_ref/xch4 * atm.CH4
    atm.H2O = xh2o_ref/xh2o * atm.H2O

    # xco2 = np.sum(atm.CO2) / np.sum(atm.air) * 1.E6
    # xch4 = np.sum(atm.CH4) / np.sum(atm.air) * 1.E9
    # xh2o = np.sum(atm.H2O) / np.sum(atm.air) * 1.E6
    # print('XCO2 ', xco2)
    # print('XCH4 ', xch4)
    # print('XH2O ', xh2o)
    # print('H2O col',np.sum(atm.H2O))
    return(atm)


class atmosphere_data:
    """
    The atmosphere_data class collects data of
    the thermodynamic state and the composition of
    the atmosphere.

    CONTAINS

    method __init__(self,zlay,zlev) |
    method get_data_AFGL(self,filename) |
        method get_data_ECMWF_ads_egg4(self,filename, month, longitude, latitude)
    """
    ###########################################################

    def __init__(self, zlay, zlev, psurf):
        """
        init class

        Parameters
        ----------
        zlay: array of vertical height layers, mid-points [nlay] [m]
        zlev: array of vertical height levels, boundaries [nlev=nlay+1] [m]
        psurf: scalar of surface pressure [Pa]

        Returns
        -------
        .zlev: array of vertical height levels, boundaries [nlev=nlay+1] [m]
        .zlay: array of vertical height layers [nlay] [m]
        """
        self.__setattr__('zlay', zlay)
        self.__setattr__('dzlay', zlev[0:zlev.size-1]-zlev[1:zlev.size])
        self.__setattr__('zlev', zlev)
        self.__setattr__('psurf', psurf)

    ###########################################################

    def get_data_AFGL(self, filename):
        """
        Read atmospheric data from AFGL database
        file, interpolate on output height grid

        Parameters
        ----------
        filename: str
                file with AFGL data

        Returns
        -------
        .tlev: temperature level profile [nlev] [K]
        .tlay: temperature layer profile [nlev] [K]
        .plev: pressure level profile [nlev] [hPa]
        .play: pressure layer profile [nlay] [hPa]
        .air:  air partial column profile [nlay] [#/m^2]
        .O3:   o3 partial column profile [nlay] [#/m^2]
        .H2O:  h2o prtial column profile [nlay] [#/m^2]
        .CO2:  co2 partial column profile [nlay] [#/m^2]
        .NO2:  no2 partial column profile [nlay] [#/m^2]
        .O2:   o2 partial column profile [nlay] [#/m^2]
        .CH4:  CH4 partial column profile [nlay] [#/m^2]
        """
        # check whether input is in range
        while True:
            if os.path.exists(filename):
                break
            else:
                sys.exit("ERROR! atmosphere_data.get_data_AFGL: file does not exist.")

        # Read AFGL file
        atm_in = np.genfromtxt(filename, skip_header=2)

        zalt_in = atm_in[:, 0]*1.E3         # height [km] -> [m]
        press_in = atm_in[:, 1]*1.E2        # pressure [Pa]
        # print("pressure",press_in)
        temp_in = atm_in[:, 2]              # temperature [K]
        # print("temperature",temp_in)
        air_in = atm_in[:, 3]*1E6           # air number density [#/m^3]
        # print("air_in",air_in)
        o3_in = atm_in[:, 4]*1E6/air_in     # o3 number density -> mole fraction [-]
        o2_in = atm_in[:, 5]*1E6/air_in     # o2 number density -> mole fraction [-]
        h2o_in = atm_in[:, 6]*1E6/air_in    # h2o number density -> mole fraction [-]
        # print("h2o",h2o_in)
        co2_in = atm_in[:, 7]*1E6/air_in    # co2 number density -> mole fraction [-]
        no2_in = atm_in[:, 8]*1E6/air_in    # no2 number density -> mole fraction [-]
        nlev_in = zalt_in.size              # number of input levels

        # [#/m^2 * 1/hPa] air column above P is P*NA/constants.MDRYAIR/g from p = m*kg/area
        sp = constants.NA/(constants.MDRYAIR*constants.g0)

        # truncate or extrapolate the AFGL profile depending on psurf
        # print(press_in[press_in.size-1], self.psurf)
        if press_in[press_in.size-1] < self.psurf:
            # extrapolation required
            dz = np.log(self.psurf/press_in[press_in.size-1])*constants.Rgas * \
                temp_in[temp_in.size-1]/(constants.grav*constants.MDRYAIR)
            press_in = np.append(press_in, self.psurf)
            zalt_in = np.append(zalt_in+dz, 0.)
            temp_in = np.append(temp_in, temp_in[temp_in.size-1])
            air_in = np.append(air_in, press_in[press_in.size-1]*constants.NA /
                               (constants.Rgas * temp_in[temp_in.size-1]))
            o3_in = np.append(o3_in, o3_in[o3_in.size-1])
            o2_in = np.append(o2_in, o2_in[o2_in.size-1])
            h2o_in = np.append(h2o_in, h2o_in[h2o_in.size-1])
            co2_in = np.append(co2_in, co2_in[co2_in.size-1])
            no2_in = np.append(no2_in, no2_in[no2_in.size-1])
            nlev_in = nlev_in-1
        elif press_in[press_in.size-1] > self.psurf:
            # interpolation required
            # self.psurf is in the interval [press_in[intv], press_in[intv-1]]
            intv = np.searchsorted(press_in, self.psurf)
            press_in = np.append(press_in[0:intv], self.psurf)
            temp_in = temp_in[0:intv+1]
            air_in = np.append(air_in[0:intv], press_in[press_in.size-1] *
                               constants.NA / (constants.Rgas * temp_in[temp_in.size-1]))
            o3_in = o3_in[0:intv+1]
            o2_in = o2_in[0:intv+1]
            h2o_in = h2o_in[0:intv+1]
            co2_in = co2_in[0:intv+1]
            no2_in = no2_in[0:intv+1]
            zalt_in = zalt_in[0:intv+1]
            dz = np.log(press_in[press_in.size-1]/press_in[press_in.size-2])*constants.Rgas * \
                temp_in[temp_in.size-1]/(constants.grav*constants.MDRYAIR)
            zalt_in = np.append(zalt_in[0:intv]-zalt_in[intv-1]+dz, 0)

        # Interpolate temperature [K] on output layers
        # Flip arrays because our heights are descending
        # (from top to bottom), while np.interp expects ascending order
        self.__setattr__('tlev', np.flip(np.interp(np.flip(self.zlev), np.flip(zalt_in), np.flip(temp_in))))
        self.__setattr__('tlay', np.flip(np.interp(np.flip(self.zlay), np.flip(zalt_in), np.flip(temp_in))))
        # print("temp",self.tlay,self.tlev)
        # Calculate pressure [hPa] on output levels and layers
        self.__setattr__('plev', np.flip(np.interp(np.flip(self.zlev), np.flip(zalt_in), np.flip(press_in))))
        self.__setattr__('play', np.flip(np.interp(np.flip(self.zlay), np.flip(zalt_in), np.flip(press_in))))
        # print("press",self.play,self.plev)
        # Calculate the vertical column of air above pressure level
        # and use this to calculate the partial vertical air columns per layer [#/m^2].
        # Partial columns have the advantage that multiplication with cross sections
        # yields optical depth.
        nlev = len(self.zlev)
        nlay = len(self.zlay)
        # [#/m^2 * 1/hPa] air column above P is P*NA/constants.MDRYAIR/g from p = m*g/area
        sp = constants.NA/(constants.MDRYAIR*constants.g0)
        vc_air = sp*self.plev  # air column [#/m^2] above pressure level
        self.air = (vc_air[1:nlev]-vc_air[0:nlev-1])  # [#/m^2]
        self.air[0] = vc_air[0]  # [#/m^2] uppermost layer extends to infinity in terms of number of molecules
        # print(self.air)
        # Interpolate mole fractions on output height grid
        # and then calculate partial columns per layer [#/m^2]
        # ozone
        self.__setattr__('O3', np.flip(np.interp(np.flip(self.zlay),
                                                 np.flip(zalt_in), np.flip(o3_in)))*self.air)
        # water vapor
        self.__setattr__('H2O', np.flip(np.interp(np.flip(self.zlay),
                                                  np.flip(zalt_in), np.flip(h2o_in)))*self.air)
        # co2
        self.__setattr__(
            'CO2', np.flip(np.interp(np.flip(self.zlay), np.flip(zalt_in), np.flip(co2_in)))*self.air)
        # no2
        self.__setattr__('NO2', np.flip(np.interp(np.flip(self.zlay),
                                                  np.flip(zalt_in), np.flip(no2_in)))*self.air)
        # o2 use a constant mixing ratio
        self.__setattr__('O2', constants.XO2*self.air)
        # ch4, it is not in the AFGL record so we assume a constant mixing ratio with altitude.
        self.__setattr__('CH4', constants.XCH4*self.air)

        return

    ###########################################################

    def get_data_ECMWF_ads_egg4(self, filename, month, longitude, latitude):
        """
        Read atmospheric data provided by ECMWF ADS EGG4 run:
        https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-ghg-reanalysis-egg4-monthly

        Download script:

        import cdsapi
        c = cdsapi.Client()
        c.retrieve(
             'cams-global-ghg-reanalysis-egg4-monthly',
             {
                 'variable': [
                     'carbon_dioxide', 'geopotential', 'methane',
                     'relative_humidity', 'temperature',
                 ],
                 'pressure_level': [
                     '1', '2', '3',
                     '5', '7', '10',
                     '20', '30', '50',
                     '70', '100', '150',
                     '200', '250', '300',
                     '400', '500', '600',
                     '700', '800', '850',
                     '900', '925', '950',
                     '1000',
                 ],
                 'year': '2016',
                 'month': [
                     '01', '02', '03',
                     '04', '05', '06',
                     '07', '08', '09',
                     '10', '11', '12',
                 ],
                 'product_type': 'monthly_mean',
                 'format': 'netcdf',
             },
             'download.nc'

        Parameters
        ----------
        filename: str
                filepath to ECMWF netcdf file
        month: int
                [01 ... 12]
        longitude: int
                [0 ... 360 degree]
        latitude: int
                 [-90 ... 90 degree]
        Returns
        -------
        temp: temperature profile [nlay] [K]
        plev: pressure level profile [nlev] [hPa]
        atmo[tlev]: temperature level profile [nlev] [K]
        atmo[tlay]: temperature layer profile [nlev] [K]
        atmo[plev]: pressure level profile [nlev] [hPa]
        atmo[play]: pressure layer profile [nlay] [hPa]
        atmo[AIR]: air partial column profile [nlay] [#/m^2]
        atmo[H2O]: h2o partial column profile [nlay] [#/m^2]
        atmo[CO2]: co2 partial column profile [nlay] [#/m^2]
        atmo[CH4]: no2 partial column profile [nlay] [#/m^2]
        atmo[O2]: o2 partial column profile [nlay] [#/m^2]
        """
        # check whether input is in range
        while True:
            if os.path.exists(filename) and 1 <= month <= 12 and -90. <= latitude <= 90. and 0. <= longitude <= 360.:
                break
            elif not os.path.exists(filename):
                print("ERROR! read_ecmwf_ads_egg4: filename does not exist.")
                raise StopExecution
            else:
                print("ERROR! read_ecmwf_ads_egg4: input out of range.")
                raise StopExecution

        # Open netcdf file
        ds = nc.Dataset(filename)
        # print(ds.variables)

        # Select month index, latitude/longitude index (next neighbour)
        itime = int(month-1)
        ilat = np.argmin(abs(ds['latitude'][:] - latitude))
        ilon = np.argmin(abs(ds['longitude'][:] - longitude))

        # ECMWF: Geopotential [m2 s-2] converted to height [m], approximate use of g0
        zalt_in = np.array([d/constants.g0 for d in ds['z'][itime, :, ilat, ilon]])
        nlev_in = zalt_in.size
        # ECMWF: Pressure [hPa]
        press_in = ds['level'][:]
        # ECMWF: Temperature [K]
        temp_in = ds['t'][itime, :, ilat, ilon]
        # ECMWF: Humidity [%] converted to water vapor mole fraction [mol/mol] via Clausius-Clapeyron equation
        pS = [clausius_clapeyron(Ti) for Ti in temp_in]
        # ECMWF: Mole fraction is partial pressure over dry total pressure, partial pressure is rel. hum. * sat. vapor pressure
        h2o_in = np.array([d/100.*pSi/(pi - d/100.*pSi)
                          for d, pSi, pi in zip(ds['r'][itime, :, ilat, ilon], pS, press_in)])
        # ECMWF: Carbon dioxide mass mixing ratio [kg kg-1] converted to mole fraction [mol/mol]
        co2_in = np.array([d/MCO2*MDRYAIR for d in ds['co2'][itime, :, ilat, ilon]])
        # ECMWF: Methane mass mixing ratio [kg kg-1] converted to mole fraction [mol/mol]
        ch4_in = np.array([d/MCH4*MDRYAIR for d in ds['ch4'][itime, :, ilat, ilon]])
        ds.close

        ################################################################################################################
        # [#/m^2 * 1/hPa] air column above P is P*NA/constants.MDRYAIR/g from p = m*g/area
        sp = constants.NA/(constants.MDRYAIR*constants.g0)*1.E2
        # air_in=sp
        # truncate or extrapolate the AFGL profile depending on psurf
        # print(press_in[press_in.size-1], self.atmo['psurf'])

        if press_in[press_in.size-1] < self.atmo['psurf']:
            # extrapolation required
            dz = np.log(self.atmo['psurf']/press_in[press_in.size-1])*constants.Rgas * \
                temp_in[temp_in.size-1]/(constants.grav*constants.MDRYAIR)
            press_in = np.append(press_in, self.atmo['psurf'])
            zalt_in = np.append(zalt_in-dz, 0.)
            temp_in = np.append(temp_in, temp_in[temp_in.size-1])
            # air_in   = np.append(air_in,press_in[press_in.size-1]*constants.NA*1.E-4/(constants.Rgas *temp_in[temp_in.size-1]))
            # o3_in    = np.append(o3_in, o3_in[o3_in.size-1])
            # o2_in    = np.append(o2_in, o2_in[o2_in.size-1])
            h2o_in = np.append(h2o_in, h2o_in[h2o_in.size-1])
            co2_in = np.append(co2_in, co2_in[co2_in.size-1])
            ch4_in = np.append(ch4_in, ch4_in[ch4_in.size-1])
            # no2_in   = np.append(no2_in, no2_in[no2_in.size-1])
            nlev_in = nlev_in-1
        elif press_in[press_in.size-1] > self.atmo['psurf']:
            # interpolation required
            # self.atmo['psurf'] is in the interval [press_in[intv], press_in[intv-1]]
            intv = np.searchsorted(press_in, self.atmo['psurf'])
            press_in = np.append(press_in[0:intv], self.atmo['psurf'])
            temp_in = temp_in[0:intv+1]
            # air_in   = np.append(air_in[0:intv],press_in[press_in.size-1]*constants.NA*1.E-4/(constants.Rgas *temp_in[temp_in.size-1]))
            # o3_in    = o3_in[0:intv+1]
            # o2_in    = o2_in[0:intv+1]
            h2o_in = h2o_in[0:intv+1]
            co2_in = co2_in[0:intv+1]
            ch4_in = ch4_in[0:intv+1]
            # no2_in   = no2_in[0:intv+1]
            zalt_in = zalt_in[0:intv+1]
            dz = np.log(press_in[press_in.size-1]/press_in[press_in.size-2])*constants.Rgas * \
                temp_in[temp_in.size-1]/(constants.grav*constants.MDRYAIR)
            zalt_in = np.append(zalt_in[0:intv]-zalt_in[intv-1]+dz, 0)

        # Interpolate temperature [K] on output layers
        # Flip arrays because our heights are descending
        # (from top to bottom), while np.interp expects ascending order
        self.atmo['tlay'] = np.flip(np.interp(np.flip(self.atmo['zlay']), np.flip(zalt_in), np.flip(temp_in)))
        self.atmo['tlev'] = np.flip(np.interp(np.flip(self.atmo['zlev']), np.flip(zalt_in), np.flip(temp_in)))

        # Calculate pressure [hPa]
        self.atmo['play'] = np.flip(np.interp(np.flip(self.atmo['zlay']), np.flip(zalt_in), np.flip(press_in)))
        self.atmo['plev'] = np.flip(np.interp(np.flip(self.atmo['zlev']), np.flip(zalt_in), np.flip(press_in)))

        # Calculate the vertical column of air above pressure level
        # and use this to calculate the partial vertical air columns per layer [#/m^2].
        # Partial columns have the advantage that multiplication with cross sections
        # yields optical depth.
        nlev = len(self.atmo['zlev'])
        sp = NA/(MDRYAIR*g0)*1.E2  # [#/m^2 * 1/hPa] air column above P is P*NA/Mair/g from p = m*g/area
        vc_air = sp*self.atmo['plev']  # air column [#/m^2] above pressure level
        self.atmo['AIR'] = (vc_air[1:nlev]-vc_air[0:nlev-1])  # [#/m^2]
        self.atmo['AIR'][0] = vc_air[0]  # [#/m^2] uppermost layer extends to infinity in terms of number of molecules

        # Interpolate mole fractions on output height grid
        # and then calculate partial columns per layer [#/m^2]
        # water vapor
        self.atmo['H2O'] = np.flip(np.interp(np.flip(self.atmo['zlay']),
                                             np.flip(zalt_in), np.flip(h2o_in)))*self.atmo['AIR']
        # co2
        self.atmo['CO2'] = np.flip(np.interp(np.flip(self.atmo['zlay']),
                                             np.flip(zalt_in), np.flip(co2_in)))*self.atmo['AIR']
        # no2
        self.atmo['CH4'] = np.flip(np.interp(np.flip(self.atmo['zlay']),
                                             np.flip(zalt_in), np.flip(ch4_in)))*self.atmo['AIR']
        # self.atmo['CH4'] = XCH4*self.atmo['AIR']
        # o2 use a constant mixing ratio
        self.atmo['O2'] = constants.XO2*self.atmo['AIR']

    ###########################################################

    def get_albedo_flat(self, alb_prior):
        """
        Generate spectrally flat albedo array

        Parameters
        ----------
        alb_prior: albedo value to be used throughout spectral range
        Returns
        -------
        alb: constant albedo array [wavelength]
        """
        # check whether input is in range
        while True:
            if 0. <= alb_prior <= 1.:
                break
            else:
                print("ERROR! surface_prop.get_albedo_flat: albedo prior needs to be in [0...1].")
                raise StopExecution

        # Read data from file
        self.alb = np.array([alb_prior for wi in self.wave])

    ###########################################################
    def get_albedo_poly(self, wave, c_lst):
        """
         Parameters
        ----------
        lambda : wavelenght [nm].
        c_lst: list of constants

        Returns
        -------
        albedo polynomial dependent on wavelenght

        """

        albedo = np.zeros(len(wave))

        for i in range(0, len(c_lst)):
            albedo = albedo + c_lst[i]*(wave)**(i)

        self.alb = albedo

    ###########################################################
    def get_albedo_CUSTOM(self, filename, sfname):
        """
        Read albedo from custom database. This is generic typical
        data. For a comprehensive albedo database, have a look
        at: https://speclib.jpl.nasa.gov/

        Parameters
        ----------
        filename: str
                file with albedo database
        sfname: str
                name of surface type
                [sand,soil,snow,vegetation,water]
        Returns
        -------
        alb: albedo array interpolated to wavelength [wavelength]
        """
        # check whether input is in range
        sftypes = ['sand', 'soil', 'snow', 'vegetation', 'water']
        while True:
            if os.path.exists(filename) and sfname in sftypes:
                break
            else:
                print("ERROR! surface_prop.get_albedo_CUSTOM: input out of range.")
                raise StopExecution

        # Read data from file
        raw = np.genfromtxt(filename, skip_header=15)  # read file into numpy array, skip 15 header lines

        # Index of surface type
        isf = sftypes.index(sfname)

        # Interpolate albedo to wavelength array
        wave_org = raw[:, 0]
        alb_org = raw[:, isf+1]
        self.alb = np.interp(self.wave, wave_org, alb_org)

    ###########################################################
    def get_albedo_ECOSTRESS(self, filename):
        """
        Read albedo from ECOSTRESS database
        at: https://speclib.jpl.nasa.gov/

        Parameters
        ----------
        filename: file with albedo database
        Returns
        -------
        alb: albedo array interpolated to wavelength [wavelength]
        """
        # Check whether input is in range
        while True:
            if os.path.exists(filename):
                break
            else:
                print("ERROR! surface_prop.get_albedo_ECOSTRESS: input out of range.")
                raise StopExecution

        # Read data from file
        raw = np.genfromtxt(filename, skip_header=21, unpack=True)

        wv_in = np.array([a*1E3 for a in raw[0, :]])  # wavelength [nm]
        alb_in = np.array([a/1E2 for a in raw[1, :]])  # albedo [0...1]
        # Check if wavelength in ascending order. If not, flip arrays.
        if wv_in[0] > wv_in[-1]:
            # Interpolate albedo to wavelength array
            self.alb = np.interp(self.wave, np.flip(wv_in), np.flip(alb_in))
        else:
            self.alb = np.interp(self.wave, wv_in, alb_in)


def combine_meteo_standard_atm(meteodata, atm_std, gases, suffix=""):
    nalt, nact = meteodata.lat.shape
    # index of standard atmosphere that corresp
    print('Combining microHH and AFGL model aonds to TOA of microHH')
    idx = (np.abs(atm_std.zlev - np.max(meteodata.zlev))).argmin()

    if (atm_std.zlev[idx] != np.max(meteodata.zlev)):
        sys.exit('vertical grids are not aligned in libATM.combine_microHH_standard_atm')
    #    zlev_tmp = np.concatenate()
    zmeteo = meteodata.zlay
    zatm = atm_std.zlay
    nlay = zatm.size

    # define a mask for vertical integration of microHH data
    ztop = atm_std.zlev[idx: nlay]
    zbot = atm_std.zlev[idx+1: nlay+1]

    midx = {}
    for ilay in range(ztop.size):
        midx[ilay] = np.where((zmeteo < ztop[ilay]) & (zmeteo >= zbot[ilay]))

    # integrate microHH and use it to replace the lowest four values of atm_std
    atm = np.ndarray((nalt, nact), np.object_)

    # for NO we assume no contribution above upper model boundary of microHH
    for ialt in tqdm(range(nalt)):
        for iact in range(nact):

            atm[ialt, iact] = deepcopy(atm_std)
            if "NO" in gases:
                atm[ialt, iact].__setattr__('NO', np.zeros(nlay))
            if "NO2" in gases:
                atm[ialt, iact].__setattr__('NO2', np.zeros(nlay))

            for gas in gases:
                tmp = atm[ialt, iact].__getattribute__(gas.upper())
                meteo = meteodata.__getattribute__(gas+suffix)
                for ilay in range(ztop.size):
                    tmp[idx+ilay] += np.sum(meteo[ialt, iact, midx[ilay]])

    return atm


# def read_microHH_data(filename, keys_to_read, time=2):
#     """
#     Read the data

#     Parameters
#     ----------
#     filename : TYPE
#         DESCRIPTION.
#     keys_to_read : TYPE
#         DESCRIPTION.
#     time : Integer, optional
#         This represents the time index corresponding to a time. The default is 2.

#     Returns
#     -------
#     data container at a time. CO2 is in moles/m^3

#     """

#     # longitude :   Zonal location
#     # latitude :   Meridional location
#     # ps :   Surface air pressure
#     # zsurf: Surface elevation
#     # p :    Air pressure
#     # ph :   air pressure at cell edge
#     # zh :   height above surface at cell edge [m]
#     # z :    layer height above surface [m]
#     # ta :   Absolute temperature
#     # hus :  Specific humidity
#     # ua :   Zonal wind
#     # va :   Meridional wind
#     # wa :   Vertical wind
#     # CO2_PP_L :     CO2 tracer mole fraction low release
#     # CO2_PP_M :     CO2 tracer mole fraction medium release
#     # CO2_PP_H :     CO2 tracer mole fraction high release
#     # CO2_PP_M_ECW : CO2 tracer mole fraction medium release, from individual tower groups

#     # time: size of 97, starts at 2018-05-23 00:00:00, time step is 15mins
#     #       So 0 is 2018-05-23 00:00:00
#     #          1 is 2018-05-23 00:00:15
#     #          97 is 2018-05-24 00:00:00

#     fl = nc.Dataset(filename, 'r')
#     # Outside these indices in lat-lon data is masked so no use reading it
#     i1 = 36
#     i2 = 264
#     j1 = 8
#     j2 = 252

#     # this is dict
#     data = Emptyclass()

#     # read all data
#     for ky in keys_to_read:
#         # if "CO2" in ky:
#         #     co2_tmp = (fl[ky][time, :, i1:i2, j1:j2]).filled(fill_value=0)
#         #     if 'p' not in data.__dict__.keys():
#         #         data.__setattr__('p', fl['p'][time, :, i1:i2, j1:j2])
#         #     if 'ta' not in data.__dict__.keys():
#         #         data.__setattr__("ta", fl['ta'][time, :, i1:i2, j1:j2])
#         #     # pV = nRT => n/V = p/RT [N/m2  mol K/J  1/K]  = [mol/ m3]
#         #     # air_density = data.p/(data.ta * constants.Rgas)

#         #     data.__setattr__(ky, (co2_tmp*data.p)/(data.ta*constants.Rgas).data)
#         if ky in ['longitude', 'latitude']:
#             data.__setattr__(ky, (fl[ky][i1:i2, j1:j2]).data)
#         elif ky in ['ps', 'zsurf']:
#             data.__setattr__(ky, (fl[ky][time, i1:i2, j1:j2]).data)
#         else:
#             data.__setattr__(ky, fl[ky][time, :, i1:i2, j1:j2])
#     fl.close()

#     return data
