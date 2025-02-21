# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Simplified instrument model and L1B processor.

Convert the SGM output to a level 1B product. It does not include any
details on the instrument.

"""
from tqdm import tqdm
import math
import netCDF4 as nc
import numpy as np
import numpy.typing as npt

from ..lib import libNumTools
from ..lib.libWrite import writevariablefromname
from teds.gm.types import Geometry
from teds.l1al1b.python.calibration import bin_l1b
from teds.l1al1b.python.types import L1


def sparse_isrf_convolution(
        isrf: npt.NDArray[np.float64],
        mask: npt.NDArray[np.bool_],
        spectrum: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    nwav = isrf[:, 0].size
    spectrum_conv = np.empty(nwav)
    for iwav in range(nwav):
        idx = mask[iwav, :]
        spectrum_conv[iwav] = isrf[iwav, idx].dot(spectrum[idx])
    return spectrum_conv


def get_sgm_rad_data(filename: str, ialt: int) -> dict:
    input = nc.Dataset(filename, mode='r')
    sgm_data = {}
    sgm_data['wavelength line-by-line'] = input['wavelength'][:]
    sgm_data['solar irradiance line-by-line'] = input['solar_irradiance'][:]
    sgm_data['radiance line-by-line'] = input['radiance'][ialt, :, :]
    input.close()
    return sgm_data


def get_gm_data(filename: str) -> Geometry:
    nc_geo = nc.Dataset(filename)
    return Geometry(nc_geo['latitude'][:].data,
                    nc_geo['longitude'][:].data,
                    np.zeros(nc_geo['latitude'][:].shape),
                    nc_geo['solar_zenith'][:].data,
                    nc_geo['solar_azimuth'][:].data,
                    nc_geo['sensor_zenith'][:].data,
                    nc_geo['sensor_azimuth'][:].data)


def sim_output(filename: str, l1b_product: L1) -> None:
    output = nc.Dataset(filename, mode='w')
    output.title = 'Tango Carbon level 1B data'

    nalt, nact, nwav = l1b_product.spectra.shape
    output.createDimension('wavelength', nwav)
    output.createDimension('across_track_sample', nact)
    output.createDimension('along_track_sample', nalt)

    nc_grp = output.createGroup('geolocation_data')
    _dims = ('along_track_sample', 'across_track_sample')

    nc_var = nc_grp.createVariable(
        'latitude', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'latitude at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees_north'
    nc_var[:] = l1b_product.geometry.lat

    nc_var = nc_grp.createVariable(
        'longitude', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'longitude at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees_east'
    nc_var[:] = l1b_product.geometry.lon

    nc_var = nc_grp.createVariable(
        'solar_zenith', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'solar zenith angle at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_product.geometry.sza

    nc_var = nc_grp.createVariable(
        'solar_azimuth', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'solar azimuth angle at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_product.geometry.saa

    nc_var = nc_grp.createVariable(
        'sensor_zenith', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'sensor zenith angle at bin locations'
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_product.geometry.vza

    nc_var = nc_grp.createVariable(
        'sensor_azimuth', 'f4', _dims, fill_value=-32767.0)
    nc_var.long_name = 'sensor azimuth angle at bin locations'
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = 'degrees'
    nc_var[:] = l1b_product.geometry.vaa

    nc_grp = output.createGroup('observation_data')
    _dims = ('across_track_sample', 'wavelength')

    l1b_wave = np.zeros((nact, nwav))
    for iact in range(nact):
        l1b_wave[iact, :] = l1b_product.wavelengths[0, :]
    writevariablefromname(
        nc_grp, 'wavelength', _dims, l1b_product.wavelengths[0, :])

    _dims3 = ('along_track_sample', 'across_track_sample', 'wavelength')

    writevariablefromname(nc_grp, 'radiance', _dims3, l1b_product.spectra)

    nc_var = nc_grp.createVariable(
        'radiance_stdev', 'f8', _dims3, fill_value=-32767.0)
    nc_var.long_name = 'relative radiance noise (1-sigma)'
    nc_var.units = 'photons/(sr nm m2 s)'
    nc_var.valid_min = 0.0
    nc_var.valid_max = 1e24
    nc_var[:] = l1b_product.spectra_noise

    output.close()


def simplified_instrument_model_and_l1b_processor(config: dict) -> None:

    l1b_product = L1.from_empty()
    l1b_product.wavelengths = np.arange(config['spec_settings']['wave_start'],
                                        config['spec_settings']['wave_end'],
                                        config['spec_settings']['dwave'])  # nm
    l1b_product.geometry = get_gm_data(config['io_files']['input_gm'])

    # basic dimensions
    nwav = l1b_product.wavelengths.size
    nalt, nact = l1b_product.geometry.sza.shape

    # measurement array
    ymeas = np.empty([nalt, nact, nwav])

    # get line-by-line spectral grid and define some pointers

    sgm_data = get_sgm_rad_data(config['io_files']['input_sgm'], ialt=0)
    wave_lbl = sgm_data['wavelength line-by-line']
    wave = l1b_product.wavelengths

    # define isrf function
    isrf_convolution = libNumTools.get_isrf(
        wave, wave_lbl, config['isrf_settings'])

    for ialt in tqdm(range(nalt)):

        # get lbl data from sgm file for scan line ialt
        sgm_data = get_sgm_rad_data(config['io_files']['input_sgm'], ialt)

        for iact in range(nact):
            spectrum_lbl = np.array(sgm_data['radiance line-by-line'][iact, :])
            # isrf convolution
            ymeas[ialt, iact, :] = isrf_convolution(spectrum_lbl)

    # Noise model based on SNR = a I / (sqrt (aI +b )) ; a[(e- m2 sr s
    # nm) / phot], and [b] = e-.
    snr = config['snr_model']['a_snr'] * ymeas / \
        (np.sqrt(config['snr_model']['a_snr']*ymeas
                 + config['snr_model']['b_snr']))

    l1b_product.spectra_noise = ymeas/snr  # units [1]

    # Random noise
    # Get nwav random numbers that are normally distributed with
    # standard deviation 1.
    np.random.seed(config['snr_model']['seed'])

    ynoise = np.empty([nalt, nact, nwav])
    for ialt in range(nalt):
        for iact in range(nact):
            noise_dis = np.random.normal(0., 1., nwav)
            # noise contribution
            ynoise[ialt, iact, :] = (
                1 / snr[ialt, iact, :] * noise_dis * ymeas[ialt, iact, :])

    if config['sim_with_noise']:
        l1b_product.spectra = ymeas + ynoise
    else:
        l1b_product.spectra = ymeas

    l1b_product.wavelengths = np.tile(l1b_product.wavelengths,
                                      nact).reshape(nact, -1)
    if 'bin_spectra' in config:
        bin_l1b(l1b_product, config['bin_spectra'])
        l1b_product.spectra_noise /= math.sqrt(config['bin_spectra'])

    sim_output(config['io_files']['output_l1b'], l1b_product)

    print('=>siml1b calculation finished successfully')
    return
