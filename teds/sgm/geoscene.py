# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Geophysical scene generation module for different E2E simulator profiles."""
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from xarray import DataArray
import numpy as np
import numpy.typing as npt

from .atmosphere import combine_meteo_standard_atm
from .atmosphere import get_AFGL_atm_homogenous_distribution
from .atmosphere import get_atmospheric_data
from .io import write_atmosphere
from .s2 import read_albedo
from teds import log
from teds.gm.io import read_geometry
from teds.lib.io import merge_config_with_default
from teds.lib.io import print_heading
from teds.lib.io import print_system_info


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config
        Path of YAML configuration file.

    """
    for entry in ('geometry', 'afgl'):
        filename = config['io_files'][entry]
        if not Path(filename).is_file():
            log.error(f'[io_files][{entry}] ({filename}) not found')
            exit(1)


def interp_sentinel2_albedo(s2_albedos: list[DataArray],
                            lat: npt.NDArray[np.float64],
                            lon: npt.NDArray[np.float64],) -> list[DataArray]:
    if lat.shape[0] == 1:
        # Nothing to do here because there is no grid
        return s2_albedos
    s2_albedos_regridded = []
    for s2_albedo in s2_albedos:
        log.info(f'Sentinel 2 band {s2_albedo.band_label}:')

        # Change coordinate system to WGS84
        log.info('  Projecting to WSG84')
        s2_albedo = s2_albedo.rio.reproject('EPSG:4326')
        s2_albedo = s2_albedo.rename({'x': 'lon', 'y': 'lat'})

        # Extract data on target grid. Define an interpolating
        # function interp such that interp(lat,lon) is an interpolated
        # value.
        log.info('  Interpolating to geometry grid')
        interp = RegularGridInterpolator(
            (s2_albedo.lat, s2_albedo.lon), s2_albedo.values, method='linear')

        target_points = np.array(list(zip(lat.ravel(), lon.ravel())))

        res = interp(target_points).reshape(lat.shape)

        crs = s2_albedo.rio.crs
        s2_albedo = DataArray(res,
                              dims=('y', 'x'),
                              coords={
                                  'lat': (['y', 'x'], lat),
                                  'lon': (['y', 'x'], lon)
                              },
                              attrs={
                                  'gsd': s2_albedo.gsd,
                                  'band_label': s2_albedo.band_label,
                              })
        s2_albedo.rio.write_crs(crs, inplace=True)

        # Add additional metadata
        central_wavelengths = {
            'B01': 442.1, 'B02': 492.4, 'B03': 559.8, 'B04': 664.6,
            'B05': 704.1, 'B06': 740.5, 'B07': 782.8, 'B08': 832.8,
            'B8A': 864.7, 'B09': 945.1, 'B10': 1372.5, 'B11': 1613.7,
            'B12': 2202.4,
        }
        bandwidths = {
            'B01': 21.0, 'B02': 66.0, 'B03': 36.0, 'B04': 31.0, 'B05': 15.0,
            'B06': 15.0, 'B07': 20.0, 'B08': 106.0, 'B8A': 21.0, 'B09': 20.0,
            'B10': 31.0, 'B11': 91.0, 'B12': 175.0
        }
        s2_albedo.attrs['central_wavelength'] = (
            central_wavelengths[s2_albedo.band_label])
        s2_albedo.attrs['bandwidth'] = (
            bandwidths[s2_albedo.band_label])
        s2_albedos_regridded.append(s2_albedo)

    return s2_albedos_regridded


def geoscene_generation(config_user: dict) -> None:
    """Generate a geophysical scene.

    Parameters
    ----------
    config
        Configuration dictionary

    """
    print_heading('Tango geoscene generation', empty_line=False)
    print_system_info()
    print(flush=True)

    config = merge_config_with_default(config_user, 'teds.sgm')
    check_config(config)

    geometry = read_geometry(config['io_files']['geometry'], 'extended')

    # Assume the same standard atmosphere for all pixels of the granule
    atm = get_AFGL_atm_homogenous_distribution(
        config['io_files']['afgl'],
        config['atmosphere']['n_layers'],
        config['atmosphere']['layer_thickness'],
        config['scale_gas']['xco2'],
        config['scale_gas']['xch4'],
        config['scale_gas']['xh2o'])

    # Individual spectra and single swath
    if geometry.sza.shape[0] == 1:
        # Use this profile where you want to study the retrieval
        # dependence as a function of one parameter.
        n_alt, n_act = geometry.sza.shape
        if len(config['scene_spec']['albedo']) != n_act:
            log.error(
                'input error in sgm, albedo dimension not consistent with gm')
            exit(1)
        s2_albedo = DataArray(
            np.reshape(config['scene_spec']['albedo'], (n_alt, n_act)),
            dims=('y', 'x')
        )
        s2_albedo.attrs['band_label'] = 'B11'
        s2_albedo.attrs['central_wavelength'] = 1620
        s2_albedo.attrs['bandwidth'] = '60'
        s2_albedo.rio.write_crs('EPSG:4326', inplace=True)
        s2_albedos = [s2_albedo]
    else:
        # Read Sentinel 2 albedo maps for different wavelengths. The
        # S2 data is read first because that's where we get the CRS of
        # the local UTM coordinate system which is used for
        # interpolating the MicroHH data.
        s2_albedos = read_albedo(config['io_files']['s2_albedo'])

    # Meteorological data
    meteo = get_atmospheric_data(geometry.lat,
                                 geometry.lon,
                                 s2_albedos[0].rio.crs,
                                 config['io_files']['meteo'])
    s2_albedos = interp_sentinel2_albedo(
        s2_albedos, geometry.lat, geometry.lon)

    combine_meteo_standard_atm(meteo, atm, config)

    log.info('Writing output')
    write_atmosphere(
        config['io_files']['atmosphere'], atm, meteo, s2_albedos, geometry)

    print_heading('Success')
