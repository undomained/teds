# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions for reading and writing L2 data."""
from netCDF4 import Dataset
import numpy as np

from .types import L2
from .types import RefProfiles
from teds.gm.types import Geometry
from teds.gm.io import nc_write_height
from teds.gm.io import nc_write_lat
from teds.gm.io import nc_write_lon
from teds.l1al1b.types import L1
from teds.sgm.atmosphere import Atmosphere


def write_l2(filename: str,
             atm: Atmosphere,
             l2: L2,
             ref_profiles: RefProfiles,
             geometry: Geometry) -> None:
    """Write L2 data product to NetCDF file."""
    default_fill = -32767
    n_alt, n_act = geometry.lat.shape
    n_lay = l2.col_avg_kernels['CO2'].shape[2]

    nc = Dataset(filename, mode='w')
    nc.title = 'Tango Carbon level 2 data'
    nc.createDimension('along_track_sample', n_alt)
    nc.createDimension('across_track_sample', n_act)
    nc.createDimension('layers', n_lay)
    _dims = ('along_track_sample', 'across_track_sample')
    _dims3d = ('along_track_sample', 'across_track_sample', 'layers')

    var = nc.createVariable('zlay', 'f8', 'layers', fill_value=default_fill)
    var.long_name = 'central layer height'
    var.units = 'm'
    var.valid_min = 0.0
    var.valid_max = 1e+5
    var[:] = atm.zlay

    grp = nc.createGroup('geolocation_data')
    geometry.rad2deg()
    nc_write_lat(grp, geometry.lat)
    nc_write_lon(grp, geometry.lon)
    nc_write_height(grp, geometry.height)

    var = nc.createVariable(
        'aquisition_time', 'f8', 'along_track_sample', fill_value=default_fill)
    var.long_name = 'aquisition time of sensing in UTC'
    var.units = 's'
    var.valid_min = 0.0
    var.valid_max = 1e+25
    var[:] = np.full(n_alt, 99999)

    var = nc.createVariable('albedo', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'wavelength-independent component of albedo'
    var.valid_min = 0
    var.valid_max = 1
    var[:] = l2.albedo0

    var = nc.createVariable(
        'spectral_shift', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'constant shift of spectral calibration'
    var.units = 'nm'
    var.valid_min = 0
    var.valid_max = 10
    var[:] = l2.spec_shift

    var = nc.createVariable(
        'spectral_squeeze', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'squeeze of spectral calibration'
    var.valid_min = 0
    var.valid_max = 10
    var[:] = l2.spec_squeeze

    grp_prior = nc.createGroup('prior')
    grp_ns = nc.createGroup('non_scattering_retrieval')

    scale = {'CO2': 1e6, 'CH4': 1e9, 'H2O': 1e6}

    gases = ['CO2', 'CH4', 'H2O']
    for gas in gases:
        xgas = 'X' + gas

        var = grp_ns.createVariable(
            xgas.lower(), 'f8', _dims, fill_value=default_fill)
        var.long_name = f'{gas} dry air column mixing ratio {xgas}'
        var.units = ('ppbv' if gas == 'CH4' else 'ppmv')
        var.valid_min = 0
        var.valid_max = 1e4
        var[:] = l2.mixing_ratios[gas] * scale[gas]

        var = grp_ns.createVariable(
            'precision_'+xgas.lower(), 'f8', _dims, fill_value=default_fill)
        var.long_name = (
            f'{gas} precision of dry air column mixing ratio {xgas}')
        var.units = ('ppbv' if gas == 'CH4' else 'ppmv')
        var.valid_min = 0
        var.valid_max = 1e4
        var[:] = l2.precisions[gas] * scale[gas]

        var = nc.createVariable('col_avg_kernel_'+xgas.lower(),
                                'f8',
                                _dims3d,
                                fill_value=default_fill)
        var.long_name = f'{gas} column averaging kernel of {xgas}'
        var.valid_min = 0
        var.valid_max = 10
        var[:] = l2.col_avg_kernels[gas]

        var = grp_prior.createVariable(
            'apri_prof_'+gas.lower(), 'f8', 'layers', fill_value=default_fill)
        var.long_name = f'a priori profile {gas} in layer column density'
        var.units = 'molecules / cm2'
        var.valid_min = 0
        var.valid_max = 1e28
        var[:] = ref_profiles.gases[gas]

    # Next the main product, proxy
    gases = ['CO2', 'CH4']
    for gas in gases:
        xgas = ('x' + gas).lower()
        xgas_proxy = "X" + gas + " proxy"

        var = nc.createVariable(
            xgas+'_proxy', 'f8', _dims, fill_value=default_fill)
        var.long_name = f'{xgas_proxy} dry air column mixing ratio'
        var.units = ('ppbv' if gas == 'CH4' else 'ppmv')
        var.valid_min = 0
        var.valid_max = 2e20
        var[:] = l2.proxys[gas] * scale[gas]

        var = nc.createVariable(
            'precision_'+xgas+'_proxy', 'f8', _dims, fill_value=default_fill)
        var.long_name = f'{xgas_proxy} dry air column mixing ratio precision'
        var.units = ('ppbv' if gas == 'CH4' else 'ppmv')
        var.valid_min = 0
        var.valid_max = 2e20
        var[:] = l2.proxy_precisions[gas] * scale[gas]

        var = nc.createVariable(
            'accuracy_'+xgas+'_proxy', 'f8', _dims, fill_value=default_fill)
        var.long_name = f'{xgas_proxy} dry air column mixing ratio accuracy'
        var.units = ('ppbv' if gas == 'CH4' else 'ppmv')
        var.valid_min = 0
        var.valid_max = 300
        var[:] = np.full((n_alt, n_act), 99999)

        var = nc.createVariable(
            'qa_value_'+xgas+'_proxy', 'f8', _dims, fill_value=default_fill)
        var.long_name = (
            f'{xgas_proxy} dry air column mixing ratio quality value')
        var.units = ('ppbv' if gas == 'CH4' else 'ppmv')
        var.valid_min = 0
        var.valid_max = 300
        var[:] = 100

    var = grp_prior.createVariable(
        'surface_pressure', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'surface pressure'
    var.units = 'hPa'
    var.valid_min = 0
    var.valid_max = 1400
    var[:] = 1013

    grp = nc.createGroup('diagnostics')

    var = grp.createVariable(
        'processing_flag', 'i2', _dims, fill_value=default_fill)
    var.long_name = 'processing flag to indicate processor anomalies'
    var.valid_min = 0
    var.valid_max = 101
    var[:] = np.full((n_alt, n_act), 100, dtype=np.int16)

    var = grp.createVariable(
        'convergence', 'u1', _dims, fill_value=7)
    var.long_name = 'whether pixel converged'
    var.valid_min = 0
    var.valid_max = 1
    var[:] = l2.converged

    var = grp.createVariable(
        'iterations', 'i4', _dims, fill_value=default_fill)
    var.long_name = 'number of iterations'
    var.valid_min = 0
    var.valid_max = 100
    var[:] = l2.iterations

    var = grp.createVariable('chi2', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'spectral chi square of the fit'
    var.valid_min = 0
    var.valid_max = 100
    var[:] = l2.chi2

    nc.close()


def write_l2_diagnostics(filename: str, l1b: L1, l2: L2) -> None:
    """Write diagnostics not included in the level 2 data product."""
    default_fill = -32767
    n_alt, n_act, n_wave = l1b.spectra.shape

    nc = Dataset(filename, mode='w')
    nc.title = 'Tango Carbon E2ES L2 diagnostics'
    nc.createDimension('along_track_sample', n_alt)
    nc.createDimension('along_across_sample', n_act)
    nc.createDimension('wavelength', n_wave)

    _dims = ('along_track_sample', 'along_across_sample', 'wavelength')

    var = nc.createVariable(
        'wavelength', 'f8', 'wavelength', fill_value=default_fill)
    var.long_name = 'wavelengths'
    var.units = 'nm'
    var.min_value = 0.0
    var.max_value = 8000.0
    var[:] = l1b.wavelengths

    var = nc.createVariable(
        'measurement', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'spectral photon radiance'
    var.units = 'nm-1 s-1 sr-1 m-2'
    var.valid_min = 0.0
    var.valid_max = 1e20
    var[:] = l1b.spectra

    var = nc.createVariable(
        'solar_irradiance', 'f8', 'wavelength', fill_value=default_fill)
    var.long_name = 'solar irradiance convolved with ISRF'
    var.units = 'photons / (nm m2 s)'
    var.min_value = 0.0
    var.max_value = 1e30
    var[:] = l1b.solar_irradiance

    var = nc.createVariable(
        'gain_co2', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'CO2 spectral gain vector'
    var.units = 'ppm/(photons/(nm m2 s sr))'
    var.min_value = -1e30
    var.max_value = 1e30
    var[:] = l2.gains['CO2']

    var = nc.createVariable(
        'gain_ch4', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'CH4 spectral gain vector'
    var.units = 'ppm/(photons/(nm m2 s sr))'
    var.min_value = -1e30
    var.max_value = 1e30
    var[:] = l2.gains['CH4']

    var = nc.createVariable(
        'gain_h2o', 'f8', _dims, fill_value=default_fill)
    var.long_name = 'H2O spectral gain vector'
    var.units = 'ppm/(photons/(nm m2 s sr))'
    var.min_value = -1e30
    var.max_value = 1e30
    var[:] = l2.gains['H2O']

    nc.close()
