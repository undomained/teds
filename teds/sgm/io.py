# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
from netCDF4 import Dataset
from xarray import DataArray
import numpy as np
import numpy.typing as npt
import yaml

from .atmosphere import Atmosphere
from .types import Gas
from .types import Meteo
from teds.gm.io import nc_write_geometry
from teds.gm.types import Geometry
from teds.l1al1b.types import L1


def write_atmosphere(filename: str,
                     atm: Atmosphere,
                     meteo: Meteo,
                     s2_albedos: list[DataArray],
                     geometry: Geometry) -> None:
    """Write geophysical scene data to output."""
    default_fill = -32767
    n_alt, n_act = geometry.lat.shape
    n_lay = len(atm.zlay)
    n_lev = n_lay + 1

    nc = Dataset(filename, mode='w')
    nc.title = 'Tango Carbon E2ES SGM atmospheric scene'
    nc.createDimension('bins_along_track', n_alt)
    nc.createDimension('bins_across_track', n_act)
    nc.createDimension('layers', n_lay)
    nc.createDimension('levels', n_lev)
    nc.createDimension('vector', 3)  # For source coordinates

    _dims2d = ('bins_along_track', 'bins_across_track')
    _dims3d = ('bins_along_track', 'bins_across_track', 'layers')

    # Atmosphere levels
    var = nc.createVariable('zlev', 'f8', 'levels', fill_value=default_fill)
    var.long_name = 'level height'
    var.units = 'm'
    var.valid_min = 0.0
    var.valid_max = 1e+5
    var[:] = atm.zlev

    var = nc.createVariable('zlay', 'f8', 'layers', fill_value=default_fill)
    var.long_name = 'central layer height'
    var.units = 'm'
    var.valid_min = 0.0
    var.valid_max = 1e+5
    var[:] = atm.zlay

    for gas in atm.gases:
        if not gas.name.startswith('X'):
            var = nc.createVariable('dcol_'+gas.name.lower(),
                                    'f8',
                                    _dims3d,
                                    fill_value=default_fill)
            var.long_name = gas.name + ' layer column density'
            var.units = 'molec./m2'
            var.valid_min = 0.0
            var.valid_max = 1e+28
            var[:] = gas.concentration
        else:
            var = nc.createVariable(
                gas.name.lower(), 'f8', _dims2d, fill_value=default_fill)
            var.long_name = gas.name + ' dry air mol column mixing ratio'
            var.valid_min = 0.0
            var.valid_max = 1e+5
            if gas.name in ('xch4', 'xo3', 'xno', 'xno2'):
                var.units = 'ppbv'
                var[:] = gas.concentration * 1e9
            else:
                var.units = 'ppmv'
                var[:] = gas.concentration * 1e6

        if gas.source is not None:
            var = nc.createVariable(gas.name.lower()+'_source_location',
                                    'f8',
                                    'vector',
                                    fill_value=default_fill)
            var.long_name = 'CO2 source location [height, latitude, longitude]'
            var.units = '[m, degree, degree]'
            var.valid_min = -180.0
            var.valid_max = 10000.0
            var[:] = gas.source

        if gas.emission_in_kgps is not None:
            var = nc.createVariable(
                gas.name.lower()+'_emission', 'f8', fill_value=default_fill)
            var.long_name = 'CO2 emission'
            var.units = 'kg/s'
            var.valid_min = 0.0
            var.valid_max = 10000.0
            var[:] = gas.emission_in_kgps

    # Air
    var = nc.createVariable('col_air', 'f8', _dims2d, fill_value=default_fill)
    var.long_name = 'air total column density'
    var.units = 'molec./m2'
    var.valid_min = 0.0
    var.valid_max = 1e+32
    var[:] = np.sum(atm.air)

    var = nc.createVariable('dcol_air', 'f8', _dims3d, fill_value=default_fill)
    var.long_name = 'dry air layer column density'
    var.units = 'molec./m2'
    var.valid_min = 0.0
    var.valid_max = 1e+35
    var[:] = atm.air

    # Albedos
    for s2_albedo in s2_albedos:
        var_name = 'albedo_' + s2_albedo.band_label.lower()
        var = nc.createVariable(
            var_name, 'f8', _dims2d, fill_value=default_fill)
        var.long_name = f'S2 albedo band {s2_albedo.band_label}'
        var.valid_min = 0.0
        var.valid_max = 1.0
        var.central_wavelength = s2_albedo.central_wavelength
        var.band_width = s2_albedo.bandwidth
        var[:] = s2_albedo.values

    # Coordinates
    if meteo.x.size > 0:
        var = nc.createVariable(
            'x', 'f8', 'bins_across_track', fill_value=default_fill)
        var.long_name = 'x-coordinate in UTM coordinate system'
        var.crs = str(meteo.crs)
        var.valid_min = 0.0
        var.valid_max = 1e7
        var[:] = meteo.x[0, :]

        var = nc.createVariable(
            'y', 'f8', 'bins_along_track', fill_value=default_fill)
        var.long_name = 'y-coordinate in UTM coordinate system'
        var.crs = str(meteo.crs)
        var.valid_min = 0.0
        var.valid_max = 1e7
        var[:] = meteo.y[:, 0]

    var = nc.createVariable('latitude', 'f8', _dims2d)
    var.long_name = 'latitudes'
    var.units = 'degrees'
    var.valid_min = -90.0
    var.valid_max = 90.0
    var[:] = geometry.lat

    var = nc.createVariable('longitude', 'f8', _dims2d)
    var.long_name = 'longitudes'
    var.units = 'degrees'
    var.valid_min = -180.0
    var.valid_max = 180.0
    var[:] = geometry.lon

    nc.close()


def write_atmosphere_ref(filename: str,
                         atm: Atmosphere,
                         albedo: npt.NDArray[np.floating],
                         geometry: Geometry,
                         bin_alt: int,
                         bin_act: int) -> None:
    """Bin and write geophysical reference scene data to output.

    Data is binned before saving to file. Binning factors should be
    chosen such that the final dimensions makes sense for L2
    comparison.

    Parameters
    ----------
    filename
        Reference atmosphere file name
    atm
        Atmosphere
    albedo
        Albedo
    geometry
        Extended geometry
    bin_alt
        Bin data in ALT dimension before saving to file
    bin_act
        Bin data in ACT dimension before saving to file

    """
    n_alt, n_act, n_lay = atm.get_gas('co2').concentration.shape
    n_alt_binned = int(n_alt // bin_alt)
    n_act_binned = int(n_act // bin_act)

    def bin_data(arr: npt.NDArray[np.floating]) -> npt.NDArray[np.floating]:
        """Bin array in both ALT and ACT dimensions."""
        shape = np.array(arr.shape, dtype=int)
        shape[1] = n_act_binned
        arr_binned_act = np.empty(shape)
        for i in range(n_act_binned):
            beg, end = i*bin_act, min(n_act, (i+1)*bin_act)
            arr_binned_act[:, i, ...] = np.mean(arr[:, beg:end, ...], axis=1)
        shape[0] = n_alt_binned
        arr_binned = np.empty(shape)
        for i in range(n_alt_binned):
            beg, end = min(n_alt, i*bin_alt), (i+1)*bin_alt
            arr_binned[i, ...] = np.mean(arr_binned_act[beg:end, ...], axis=0)
        return arr_binned

    # Geometry
    geometry_binned = Geometry.from_shape((n_alt_binned, n_act_binned))
    geometry_binned.lat[:] = bin_data(geometry.lat)
    geometry_binned.lon[:] = bin_data(geometry.lon)
    geometry_binned.height[:] = bin_data(geometry.height)
    geometry_binned.sza[:] = bin_data(geometry.sza)
    geometry_binned.saa[:] = bin_data(geometry.saa)
    geometry_binned.vza[:] = bin_data(geometry.vza)
    geometry_binned.vaa[:] = bin_data(geometry.vaa)

    # Atmosphere
    atm_binned = Atmosphere.from_empty()
    atm_binned.zlay = atm.zlay
    atm_binned.zlev = atm.zlev
    for gas in atm.gases:
        atm_binned.gases.append(Gas(gas.name,
                                    gas.source,
                                    gas.emission_in_kgps,
                                    bin_data(gas.concentration)))
    atm_binned.air = bin_data(atm.air)

    # Meteorological data (nothing to do here)
    meteo = Meteo.from_empty()

    # Albedo
    albedo_binned = DataArray(bin_data(albedo), dims=('y', 'x'))
    albedo_binned.attrs['band_label'] = 'B11'
    albedo_binned.attrs['central_wavelength'] = 1620
    albedo_binned.attrs['bandwidth'] = '60'
    albedo_binned.rio.write_crs('EPSG:4326', inplace=True)
    albedos_binned = [albedo_binned]

    write_atmosphere(filename,
                     atm_binned,
                     meteo,
                     albedos_binned,
                     geometry_binned)


def write_radiance(filename: str,
                   config: dict,
                   rad: L1,
                   hetero_isrf: DataArray,
                   stokes_Q: npt.NDArray[np.floating],
                   stokes_U: npt.NDArray[np.floating]) -> None:
    """Write SGM radiances to NetCDF file."""
    default_fill = -32767
    n_alt, n_act, n_wav = rad.spectra.shape

    nc = Dataset(filename, mode='w')
    nc.title = 'Tango Carbon E2ES SGM radiometric scene'
    nc.product_type = str(rad.proc_level)
    nc.createDimension('along_track_sample', n_alt)
    nc.createDimension('across_track_sample', n_act)
    nc.createDimension('wavelength', n_wav)

    var_config = nc.createVariable('configuration', str)
    config_text = (
        yaml.dump(config).replace(' true', ' yes').replace(' false', ' no'))
    var_config[:] = np.array([config_text], dtype='object')
    var_config.comment = 'configuration parameters used to produce this file'

    var = nc.createVariable(
        'wavelength', 'f8', 'wavelength', fill_value=default_fill)
    var.long_name = 'wavelengths'
    var.units = 'nm'
    var.min_value = 0.0
    var.max_value = 8000.0
    var[:] = rad.wavelengths

    var = nc.createVariable(
        'solar_irradiance', 'f8', 'wavelength', fill_value=default_fill)
    var.long_name = 'solar irradiance line-by-line'
    var.units = 'photons / (nm m2 s)'
    var.min_value = 0.0
    var.max_value = 1e30
    var[:] = rad.solar_irradiance

    var = nc.createVariable(
        'radiance',
        'f8',
        ('along_track_sample', 'across_track_sample', 'wavelength'),
        fill_value=default_fill)
    var.long_name = 'radiance line-by-line'
    var.units = 'photons / (sr nm m2 s)'
    var.min_value = 0.0
    var.max_value = 1e28
    var[:] = rad.spectra

    if stokes_Q.size > 0:
        var = nc.createVariable(
            'q',
            'f8',
            ('along_track_sample', 'across_track_sample', 'wavelength'),
            fill_value=default_fill)
        var.long_name = 'Stokes Q parameter'
        var.units = 'photons / (sr nm m2 s)'
        var.min_value = 0.0
        var.max_value = 1e28
        var[:] = stokes_Q

        var = nc.createVariable(
            'u',
            'f8',
            ('along_track_sample', 'across_track_sample', 'wavelength'),
            fill_value=default_fill)
        var.long_name = 'Stokes U parameter'
        var.units = 'photons / (sr nm m2 s)'
        var.min_value = 0.0
        var.max_value = 1e28
        var[:] = stokes_U

    nc_write_geometry(nc, rad.geometry, False)

    if hetero_isrf.shape:
        grp = nc.createGroup('isrf')
        grp.description = ('ALT and ACT dependent ISRF due to inhomogeneous '
                           'slit illumination')
        grp.createDimension('wavelength', len(hetero_isrf.wavelength.data))
        var = grp.createVariable('wavelength', 'f8', 'wavelength')
        var.min_value = -100.0
        var.max_value = 100.0
        var.units = 'nm'
        var[:] = hetero_isrf.wavelength.data
        var = grp.createVariable('isrf', 'f8', ('along_track_sample',
                                                'across_track_sample',
                                                'wavelength'))
        var.min_value = 0.0
        var.max_value = 1.0
        var[:] = hetero_isrf.data

    nc.close()


def read_atmosphere_and_albedo(filename: str) -> tuple[
        Atmosphere, npt.NDArray[np.floating]]:

    nc = Dataset(filename)

    xco2_source = None
    xch4_source = None
    if 'xco2_source_location' in nc.variables:
        xco2_source = nc['xco2_source_location'][:].data
    if 'xch4_source_location' in nc.variables:
        xch4_source = nc['xch4_source_location'][:].data
    gases = [
        Gas('ch4', None, None, nc['dcol_ch4'][:].data),
        Gas('co2', None, None, nc['dcol_co2'][:].data),
        Gas('h2o', None, None, nc['dcol_h2o'][:].data),
        Gas('XCH4', None, None, nc['xch4'][:].data * 1e-9),
        Gas('XCO2', xco2_source, None, nc['xco2'][:].data * 1e-6),
        Gas('XH2O', xch4_source, None, nc['xh2o'][:].data),
    ]
    atm = Atmosphere.from_empty()
    atm.zlay = nc['zlay'][:].data
    atm.zlev = nc['zlev'][:].data
    atm.gases = gases
    atm.air = nc['dcol_air'][:].data

    albedo = nc['albedo_b11'][:].data

    return atm, albedo
