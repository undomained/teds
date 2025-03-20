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

from ..lib.libWrite import writevariablefromname
from teds.gm.io import nc_write_geometry
from teds.gm.io import read_geometry
from teds.l1al1b.python.calibration import bin_l1b
from teds.l1al1b.python.types import L1
from teds.lib.convolution import KernelGauss
from teds.lib.io import print_heading


def get_sgm_rad_data(filename: str, ialt: int) -> dict:
    input = nc.Dataset(filename, mode='r')
    sgm_data = {}
    sgm_data['wavelength line-by-line'] = input['wavelength'][:]
    sgm_data['solar irradiance line-by-line'] = input['solar_irradiance'][:]
    sgm_data['radiance line-by-line'] = input['radiance'][ialt, :, :]
    input.close()
    return sgm_data


def sim_output(filename: str, l1b_product: L1) -> None:
    output = nc.Dataset(filename, mode='w')
    output.title = 'Tango Carbon level 1B data'

    n_alt, n_act, n_wav = l1b_product.spectra.shape
    output.createDimension('wavelength', n_wav)
    output.createDimension('across_track_sample', n_act)
    output.createDimension('along_track_sample', n_alt)

    nc_grp = output.createGroup('observation_data')
    _dims = ('across_track_sample', 'wavelength')

    l1b_wave = np.zeros((n_act, n_wav))
    for iact in range(n_act):
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

    nc_write_geometry(output.createGroup('geolocation_data'),
                      l1b_product.geometry,
                      False)

    output.close()


def simplified_instrument_model_and_l1b_processor(config: dict) -> None:

    l1b_product = L1.from_empty()
    l1b_product.wavelengths = np.arange(config['spec_settings']['wave_start'],
                                        config['spec_settings']['wave_end'],
                                        config['spec_settings']['dwave'])  # nm
    l1b_product.geometry = read_geometry(config['io_files']['input_sgm'])

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
    isrf = KernelGauss(wave_lbl.data,
                       wave,
                       config['isrf_settings']['fwhm'],
                       config['isrf_settings']['shape'])

    for ialt in tqdm(range(nalt)):

        # get lbl data from sgm file for scan line ialt
        sgm_data = get_sgm_rad_data(config['io_files']['input_sgm'], ialt)

        for iact in range(nact):
            spectrum_lbl = np.array(sgm_data['radiance line-by-line'][iact, :])
            # isrf convolution
            ymeas[ialt, iact, :] = isrf.convolve(spectrum_lbl)

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

    print_heading('Success')
