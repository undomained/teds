# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Tools to generate model atmosphere."""
from dataclasses import dataclass
from pyproj import Transformer
from scipy.interpolate import interpn
from tqdm import tqdm
from typing import Self
import numpy as np
import numpy.typing as npt

from .types import Gas
from .types import Meteo
from teds import log
from teds.lib import constants


@dataclass
class Atmosphere:
    """Thermodynamic state and composition of the atmosphere."""
    zlay: npt.NDArray[np.float64]
    dzlay: npt.NDArray[np.float64]
    zlev: npt.NDArray[np.float64]
    psurf: float
    tlev: npt.NDArray[np.float64]
    tlay: npt.NDArray[np.float64]
    plev: npt.NDArray[np.float64]
    play: npt.NDArray[np.float64]
    air: npt.NDArray[np.float64]
    gases: list[Gas]

    @classmethod
    def from_empty(cls) -> Self:
        return cls(np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   0.0,
                   np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   [])

    @classmethod
    def from_file(cls,
                  zlay: npt.NDArray[np.float64],
                  zlev: npt.NDArray[np.float64],
                  psurf: float,
                  filename: str) -> Self:
        """Read atmospheric data from AFGL database file interpolate
        to a height grid.

        tlev: temperature level profile [nlev] [K]
        tlay: temperature layer profile [nlev] [K]
        plev: pressure level profile [nlev] [hPa]
        play: pressure layer profile [nlay] [hPa]
        air:  air partial column profile [nlay] [#/m^2]
        O3:   o3 partial column profile [nlay] [#/m^2]
        H2O:  h2o prtial column profile [nlay] [#/m^2]
        CO2:  co2 partial column profile [nlay] [#/m^2]
        NO2:  no2 partial column profile [nlay] [#/m^2]
        O2:   o2 partial column profile [nlay] [#/m^2]
        CH4:  CH4 partial column profile [nlay] [#/m^2]

        Parameters
        ----------
        zlay
            Vertical height layers, mid-points [nlay] [m]
        zlev
            array of vertical height levels, boundaries [nlev=nlay+1] [m]
        psurf
            Scalar of surface pressure [Pa]
        filename
            file with AFGL data

        """
        atm = cls.from_empty()
        # Read AFGL file
        atm_in = np.genfromtxt(filename, skip_header=2)

        # Height [km] -> [m]
        zalt_in = atm_in[:, 0] * 1e3
        # Pressure [Pa]
        press_in = atm_in[:, 1] * 1e2
        # Temperature [K]
        temp_in = atm_in[:, 2]
        # Air number density [#/cm^3]
        air_in = atm_in[:, 3]
        # O3 number density -> mole fraction [-]
        o3_in = atm_in[:, 4] / air_in
        # O2 number density -> mole fraction [-]
        o2_in = atm_in[:, 5] / air_in
        # H2O number density -> mole fraction [-]
        h2o_in = atm_in[:, 6] / air_in
        # CO2 number density -> mole fraction [-]
        co2_in = atm_in[:, 7] / air_in
        # NO2 number density -> mole fraction [-]
        no2_in = atm_in[:, 8] / air_in
        # Number of input levels
        nlev_in = zalt_in.size

        # [#/m^2 * 1/hPa] air column above P is
        # P*NA/constants.MDRYAIR/g from p = m*kg/area
        sp = constants.NA/(constants.MDRYAIR*constants.g0)

        atm.zlay = zlay
        atm.dzlay = zlev[0:zlev.size-1] - zlev[1:zlev.size]
        atm.zlev = zlev
        atm.psurf = psurf

        # Truncate or extrapolate the AFGL profile depending on
        # surface pressure.
        if press_in[-1] < atm.psurf:
            # Extrapolation required
            dz = (np.log(atm.psurf / press_in[-1]) * constants.Rgas
                  * temp_in[-1] / (constants.grav * constants.MDRYAIR))
            press_in = np.append(press_in, atm.psurf)
            zalt_in = np.append(zalt_in + dz, 0.0)
            temp_in = np.append(temp_in, temp_in[-1])
            o3_in = np.append(o3_in, o3_in[-1])
            o2_in = np.append(o2_in, o2_in[-1])
            h2o_in = np.append(h2o_in, h2o_in[-1])
            co2_in = np.append(co2_in, co2_in[-1])
            no2_in = np.append(no2_in, no2_in[-1])
            nlev_in = nlev_in - 1
        elif press_in[-1] > atm.psurf:
            # Interpolation required
            # atm.psurf is in the interval [press_in[intv], press_in[intv-1]]
            intv = np.searchsorted(press_in, atm.psurf)
            press_in = np.append(press_in[0:intv], atm.psurf)
            temp_in = temp_in[0:intv+1]
            o3_in = o3_in[0:intv+1]
            o2_in = o2_in[0:intv+1]
            h2o_in = h2o_in[0:intv+1]
            co2_in = co2_in[0:intv+1]
            no2_in = no2_in[0:intv+1]
            zalt_in = zalt_in[0:intv+1]
            dz = (np.log(press_in[-1] / press_in[press_in.size-2])
                  * constants.Rgas * temp_in[-1]
                  / (constants.grav * constants.MDRYAIR))
            zalt_in = np.append(zalt_in[0:intv]-zalt_in[intv-1]+dz, 0)

        # Interpolate temperature [K] on output layers. Flip arrays
        # because our heights are descending (from top to bottom),
        # while np.interp expects ascending order.
        atm.tlev = np.flip(
            np.interp(np.flip(atm.zlev), np.flip(zalt_in), np.flip(temp_in)))
        atm.tlay = np.flip(
            np.interp(np.flip(atm.zlay), np.flip(zalt_in), np.flip(temp_in)))
        # Calculate pressure [hPa] on output levels and layers
        atm.plev = np.flip(
            np.interp(np.flip(atm.zlev), np.flip(zalt_in), np.flip(press_in)))
        atm.play = np.flip(
            np.interp(np.flip(atm.zlay), np.flip(zalt_in), np.flip(press_in)))
        # Calculate the vertical column of air above pressure level
        # and use this to calculate the partial vertical air columns
        # per layer [#/m^2]. Partial columns have the advantage that
        # multiplication with cross sections yields optical depth.
        nlev = len(atm.zlev)
        # [#/m^2 * 1/hPa] air column above P is
        # P*NA/constants.MDRYAIR/g from p = m*g/area.
        sp = constants.NA / (constants.MDRYAIR * constants.g0)
        vc_air = sp * atm.plev  # air column [#/m^2] above pressure level
        atm.air = vc_air[1:nlev] - vc_air[0:nlev-1]  # [#/m^2]
        # [#/m^2] uppermost layer extends to infinity in terms of
        # number of molecules.
        atm.air[0] = vc_air[0]
        # Interpolate mole fractions on output height grid and then
        # calculate partial columns per layer [#/m^2] ozone.
        atm.gases = [Gas('O3',
                         None,
                         None,
                         np.flip(np.interp(np.flip(atm.zlay),
                                           np.flip(zalt_in),
                                           np.flip(o3_in))) * atm.air),
                     # Water vapor
                     Gas('H2O',
                         None,
                         None,
                         np.flip(np.interp(np.flip(atm.zlay),
                                           np.flip(zalt_in),
                                           np.flip(h2o_in))) * atm.air),
                     # CO2
                     Gas('CO2',
                         None,
                         None,
                         np.flip(np.interp(np.flip(atm.zlay),
                                           np.flip(zalt_in),
                                           np.flip(co2_in))) * atm.air),
                     # NO2
                     Gas('NO2',
                         None,
                         None,
                         np.flip(np.interp(np.flip(atm.zlay),
                                           np.flip(zalt_in),
                                           np.flip(no2_in))) * atm.air),
                     # O2 use a constant mixing ratio
                     Gas('O2',
                         None,
                         None,
                         constants.XO2 * atm.air),
                     # CH4, it is not in the AFGL record so we assume
                     # a constant mixing ratio with altitude.
                     Gas('CH4',
                         None,
                         None,
                         constants.XCH4 * atm.air),
                     ]
        return atm

    def get_gas(self, name: str) -> Gas:
        res = list(filter(lambda x: x.name.lower() == name.lower(),
                          self.gases))
        if not res:
            log.error(f'object does not contain {name}')
            exit(1)
        return res[0]


def get_AFGL_atm_homogenous_distribution(
        AFGL_path: str,
        nlay: int,
        dzlay: float,
        xco2_ref: float = 405,
        xch4_ref: float = 1800.0,
        xh2o_ref: float = 1.0e4) -> Atmosphere:
    # Vertical layering assuming equidistant gridding in geometrical
    # distance
    nlev = nlay + 1  # number of levels
    # Counting from top to bottom, altitude of layer midpoint
    zlay = (np.arange(nlay - 1, -1, -1) + 0.5) * dzlay
    # Altitude of layer interfaces = levels
    zlev = np.arange(nlev - 1, -1, -1) * dzlay
    psurf = 101300.0  # Pa
    atm = Atmosphere.from_file(zlay, zlev, psurf, AFGL_path)

    for gas in atm.gases:
        if gas.name == 'CO2':
            xco2 = np.sum(gas.concentration) / np.sum(atm.air)
            gas.concentration = xco2_ref * 1e-6 / xco2 * gas.concentration
        elif gas.name == 'CH4':
            xch4 = np.sum(gas.concentration) / np.sum(atm.air)
            gas.concentration = xch4_ref * 1e-9 / xch4 * gas.concentration
        elif gas.name == 'H2O':
            xh2o = np.sum(gas.concentration) / np.sum(atm.air)
            gas.concentration = xh2o_ref * 1e-6 / xh2o * gas.concentration

    return atm


def rotate_grid(x: npt.NDArray[np.float64],
                y: npt.NDArray[np.float64],
                x_origin: float,
                y_origin: float,
                angle: float) -> tuple[npt.NDArray[np.float64],
                                       npt.NDArray[np.float64]]:
    """Rotate xy grid around an origin by an angle."""
    x0 = x - x_origin
    y0 = y - y_origin
    x_rot = x0 * np.cos(angle) - y0 * np.sin(angle)
    y_rot = x0 * np.sin(angle) + y0 * np.cos(angle)
    return x_rot + x_origin, y_rot + y_origin


def get_atmospheric_data(gm_lat: npt.NDArray[np.float64],
                         gm_lon: npt.NDArray[np.float64],
                         crs: str,
                         microhh_files: list[str]) -> Meteo:
    """Get meteorological data to same domain as input lat, lon.

    Read meteorological data and extend to given domain.

    Parameters
    ----------
    gm_lat
        Input GM Latitude
    gm_lon
        Input GM longitude
    crs
        MicroHH data is interpolated onto the GM lat/lon using a local
        UTM coordinate system given by this CRS.
    microhh_files
        List of files with plume data

    Returns
    -------
        Meteorological data (trace gas plumes).

    """
    # Read data which is already in molceules/m^2
    log.info('Gathering meteorological data')
    if microhh_files and gm_lat.shape[0] > 1:
        meteo = Meteo.from_files(microhh_files)
    else:
        meteo = Meteo.from_empty()
        meteo.lat = gm_lat
        meteo.lon = gm_lon
        return meteo
    # Rotate both GM and MicroHH grids such that they align with x and y axes
    transformer = Transformer.from_crs('EPSG:4326', crs)
    # MicroHH coordinates in UTM
    meteo_x, meteo_y = transformer.transform(meteo.lat, meteo.lon)
    # Rotation angle to align x and y axes
    phi = np.arctan2(
        meteo_y[0, -1] - meteo_y[0, 0], meteo_x[0, -1] - meteo_x[0, 0])
    # Rotate the grid around the lower left corner
    meteo_x_rot, meteo_y_rot = rotate_grid(
        meteo_x, meteo_y, meteo_x[0, 0], meteo_y[0, 0], -phi)
    # Convert to UTM and rotate the GM grid around the same point
    gm_x, gm_y = transformer.transform(gm_lat, gm_lon)
    gm_x_rot, gm_y_rot = rotate_grid(
        gm_x, gm_y, meteo_x[0, 0], meteo_y[0, 0], -phi)
    # Target grid for interpolation
    target_grid = np.column_stack((gm_x_rot.ravel(), gm_y_rot.ravel()))
    for gas in meteo.gases:
        log.info(f'  Interpolating {gas.name} from MicroHH to GM grid')
        nz = gas.concentration.shape[2]
        gas_interp = np.zeros((gm_x.shape[0], gm_x.shape[1], nz))
        for i_z in tqdm(range(nz)):
            # Take the first x and y coordinates assuming a rectilinear grid
            gas_interp[:, :, i_z] = interpn(
                (meteo_x_rot[0, :], meteo_y_rot[:, 0]),
                gas.concentration[:, :, i_z],
                target_grid,
                method='cubic',
                bounds_error=False,
                fill_value=0.0).reshape(gm_x_rot.shape)
        gas.concentration = gas_interp
    # Update MicroHH grid information
    meteo.crs = crs
    meteo.lat = gm_lat
    meteo.lon = gm_lon
    meteo.x = gm_x
    meteo.y = gm_y
    meteo.dx = (gm_x[0, -1] - gm_x[0, 0]) / gm_x.shape[1]
    meteo.dy = (gm_y[-1, 0] - gm_y[0, 0]) / gm_y.shape[0]
    # z dimension is unchanged
    return meteo


def combine_meteo_standard_atm(meteo: Meteo,
                               atm: Atmosphere,
                               config: dict) -> None:
    """Extend the atmosphere by combining it with MicroHH plumes."""
    gases_microHH = [gas.name for gas in meteo.gases]
    gases_afgl = [gas.upper() for gas in config['afgl_gases']]
    gases_all = list(set(gases_microHH) | set(gases_afgl))

    n_alt, n_act = meteo.lat.shape
    n_lay = atm.zlay.size

    # Copy standard atmosphere columns to all ground points
    for gas in atm.gases:
        gas.concentration = np.tile(
            gas.concentration, n_alt * n_act).reshape(
                (n_alt, n_act, n_lay))

    if meteo.gases:
        # Use AFGL above top boundary of microHH and add microHH at
        # lower altitudes.
        log.info('Combining microHH and AFGL model')
        # Index of standard atmosphere that corresponds to the highest
        # meteo.znodes. We shall work with lower altitudes because
        # that's where MicroHH has data. Note that atm.zlev is in
        # descending order.
        idx = len(atm.zlev) - np.searchsorted(
            atm.zlev[::-1], meteo.znodes[-1], side='right')
        if atm.zlev[idx] != np.max(meteo.znodes):
            log.error('Atmosphere and MicroHH vertical grids not aligned. '
                      f'Set layer_thickness to a multiple of {meteo.dz}')
            exit(1)
        # Define a mask for vertical integration of MicroHH data
        ztop = atm.zlev[idx: n_lay]
        zbot = atm.zlev[idx+1: n_lay+1]
        midx = {}
        for i in range(ztop.size):
            midx[i] = np.where((meteo.z >= zbot[i]) & (meteo.z < ztop[i]))[0]
        # Generate fields for MicroHH gases
        for name in gases_microHH:
            gas = meteo.get_gas(name)
            gas_averaged = np.zeros((n_alt, n_act, ztop.size))
            for i in range(ztop.size):
                gas_averaged[:, :, i] = np.sum(
                    gas.concentration[:, :, midx[i]], axis=2)
            if list(filter(lambda x: x.name == name, atm.gases)):
                atm_gas = atm.get_gas(name)
                atm_gas.concentration[:, :, idx:] += gas_averaged[:, :, :]
                atm_gas.source = gas.source
                atm_gas.emission_in_kgps = gas.emission_in_kgps
            else:
                # Not among the AFGL gases so create a new gas for the
                # atmosphere.
                gas = Gas(name, None, None, np.zeros((n_alt, n_act, n_lay)))
                gas.concentration[:, :, idx:] = gas_averaged
                atm.gases.append(gas)
    # Dry air column mixing ratios
    col_air = np.sum(atm.air)
    for name in gases_all:
        gas = atm.get_gas(name)
        atm.gases.append(Gas('X'+gas.name,
                             gas.source,
                             gas.emission_in_kgps,
                             np.sum(gas.concentration, axis=2) / col_air))
