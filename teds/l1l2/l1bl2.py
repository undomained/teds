# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Level-1B to level-2 processor.

The primary objective is to compute the dry air total column mixing
ratios XCO2 and XCH4.

"""
import os
import numpy as np
import time
from tqdm import tqdm
from copy import deepcopy

from .io import write_l2
from .io import write_l2_diagnostics
from .types import L2
from .types import RefProfiles
from teds import log
from teds.gm.io import read_geometry
from teds.l1al1b.io import read_l1
from teds.lib.convolution import Kernel
from teds.lib.gauss_newton import gauss_newton
from teds.lib.io import check_file_presence
from teds.lib.io import merge_config_with_default
from teds.lib.io import print_heading
from teds.lib.io import print_system_info
from teds.lib.radiative_transfer import MolecularData
from teds.lib.radiative_transfer import OpticAbsProp
from teds.lib.radiative_transfer import read_sun_spectrum_TSIS1HSRS
from teds.sgm import atmosphere
from teds.sgm.io import read_atmosphere_and_albedo


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config
        Configuration parameters.

    """
    check_file_presence(config['io_files']['l1b'], 'l1b')
    check_file_presence(config['io_files']['atmosphere'], 'atmosphere')
    check_file_presence(config['io_files']['afgl'], 'afgl')
    check_file_presence(config['io_files']['sun_reference'], 'sun_reference')
    if config['io_files']['isrf'] or config['isrf']['tabulated']:
        check_file_presence(config['io_files']['isrf'], 'isrf')
    if not config['xsec_forced'] and config['io_files']['dump_xsec']:
        check_file_presence(config['io_files']['dump_xsec'], 'dump_xsec')
    config['retrieval']['gases'] = [x.upper()
                                    for x in config['retrieval']['gases']]


def level1b_to_level2_processor(config_user: dict) -> None:
    """Generate L2 product for the Tango-Carbon instrument.

    Parameters
    ----------
    config
        Configuration dictionary

    """
    print_heading('Tango L1-L2 processor', empty_line=False)
    print_system_info()
    print(flush=True)

    config = merge_config_with_default(config_user, 'teds.l1l2')
    check_config(config)

    print_heading('Preparing input')

    log.info('Reading L1B product')
    l1b = read_l1(config['io_files']['l1b'], 0, None)
    geometry = read_geometry(config['io_files']['l1b'])
    # Work with radians from now on. Will be converted back to degrees
    # just before writing to output.
    geometry.deg2rad()

    # Read atmosphere for reference profiles (might not be used
    # depending on the sw_prof_init flag).
    log.info('Reading Atmosphere')
    atm_sgm, _ = read_atmosphere_and_albedo(config['io_files']['atmosphere'])

    # Define line-by-line spectral grid
    wave_start = config['spec_settings']['wave_start']
    wave_end = config['spec_settings']['wave_end']
    wave_extend = config['spec_settings']['wave_extend']
    dwave_lbl = config['spec_settings']['dwave']
    n_lbl = int((wave_end - wave_start + 2 * wave_extend) / dwave_lbl)
    wave_lbl = (wave_start - wave_extend
                + np.linspace(0, n_lbl - 1, n_lbl) * dwave_lbl)

    # Trim L1B data to match the L2 wavelength window
    i_start = np.searchsorted(l1b.wavelengths, wave_start)
    i_end = np.searchsorted(l1b.wavelengths, wave_end)
    _trim_ratio = (i_end - 1 - i_start) / len(l1b.wavelengths)  # for output
    log.info(f'Trimming L1B wavelength grid by {100 * (1-_trim_ratio):.1f}%')
    l1b.wavelengths = l1b.wavelengths[i_start:i_end+1]
    l1b.spectra = l1b.spectra[:, :, i_start:i_end+1]
    l1b.spectra_noise = l1b.spectra_noise[:, :, i_start:i_end+1]

    # Main dimensions
    n_alt, n_act, n_wave = l1b.spectra.shape
    # Vertical layering
    n_lay = config['atmosphere']['n_layers']
    dz_lay = config['atmosphere']['layer_thickness']
    p_surf = config['atmosphere']['surface_pressure']
    n_lev = n_lay + 1  # number of levels
    # Altitude of layer midpoint
    z_lay = (np.arange(n_lay - 1, -1, -1) + 0.5) * dz_lay
    # Altitude of layer interfaces = levels
    z_lev = np.arange(n_lev - 1, -1, -1) * dz_lay

    # Model atmosphere
    atm = atmosphere.Atmosphere.from_file(
        z_lay, z_lev, p_surf, config['io_files']['afgl'])

    # Calculate optical properties
    optics = OpticAbsProp(wave_lbl, z_lay)
    # If a NetCDF dump file exists read from there instead
    if (
            not os.path.exists(config['io_files']['dump_xsec'])
            or config['xsec_forced']):
        # Download molecular absorption parameters. See hapi manual
        # Sec 6.6.
        iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]
        molec = MolecularData(wave_lbl, config['io_files']['hapi'], iso_ids)
        # Molecular absorption optical properties
        optics.calc_molec_xsec(molec, atm)
        optics.xsec_to_file(config['io_files']['dump_xsec'])
    else:
        optics.xsec_from_file(config['io_files']['dump_xsec'])
        optics.check_wavelength_grid(wave_lbl)

    # Define isrf function
    if config['isrf']['tabulated']:
        log.info('Reading ISRF')
        isrf = Kernel.from_file(
            config['io_files']['isrf'], wave_lbl, l1b.wavelengths)
    else:
        log.info('Generating ISRF from generalized Gaussian parameters')
        isrf = Kernel.from_gauss(wave_lbl,
                                 l1b.wavelengths,
                                 config['isrf']['fwhm'],
                                 config['isrf']['shape'])

    # Solar irradiance: line-by-line and convolved
    sun_wavelengths, sun_spectrum = read_sun_spectrum_TSIS1HSRS(
        config['io_files']['sun_reference'])
    sun_lbl = np.interp(wave_lbl, sun_wavelengths, sun_spectrum)
    l1b.solar_irradiance = isrf.convolve(sun_lbl)

    # Need reference column profiles for fitting
    ref_profiles = RefProfiles()
    for i_gas in range(len(config['retrieval']['gases'])):
        gas = config['retrieval']['gases'][i_gas].upper()
        ref_profiles.gases[gas] = atm.get_gas(gas).concentration
        ref_profiles.initial[gas] = float(
            config['retrieval']['initial_concentrations'][i_gas])

    # Initialize L2 product
    l2 = L2(n_alt,
            n_act,
            n_wave,
            atm.zlay.size,
            config['retrieval']['n_albedos'],
            config['retrieval']['gases'])

    # Timings
    timings = {'opt': 0, 'rtm': 0, 'conv': 0, 'kern': 0}
    total_time = time.time()

    print_heading('Retrieval')
    log.info('Generating proxy product along track:')
    for i_alt in tqdm(range(n_alt)):
        for i_act in range(n_act):
            if (config['retrieval']['sw_prof_init'] == 'sgm'):
                for gas.upper in config['retrieval']['gases']:
                    ref_profiles.gases[gas] = (
                        atm_sgm.get_gas(gas).concentration[i_alt, i_act, :])

            # Initialize each retrieval with the same atmosphere
            atm_ret = deepcopy(atm)

            gauss_newton(config,
                         ref_profiles,
                         atm_ret,
                         optics,
                         wave_lbl,
                         sun_lbl,
                         l1b.solar_irradiance,
                         l1b.spectra[i_alt, i_act, :],
                         l1b.spectra_noise[i_alt, i_act, :]**2,
                         np.cos(geometry.sza[i_alt, i_act]),
                         np.cos(geometry.vza[i_alt, i_act]),
                         isrf,
                         timings,
                         l2,
                         i_alt,
                         i_act)

            if (not l2.converged[i_alt, i_act]):
                log.warning(f'Pixel ({i_alt},{i_act}) not converged')

    # Define proxy products
    xch4_model = 1.8e-6
    xco2_model = 410e-6
    l2.proxys['CO2'] = (l2.mixing_ratios['CO2'] / l2.mixing_ratios['CH4']
                        * xch4_model)
    l2.proxys['CH4'] = (l2.mixing_ratios['CH4'] / l2.mixing_ratios['CO2']
                        * xco2_model)
    rel_error = np.sqrt((l2.precisions['CO2'] / l2.mixing_ratios['CO2'])**2
                        + (l2.precisions['CH4'] / l2.mixing_ratios['CH4'])**2)
    l2.proxy_precisions['CO2'] = rel_error * l2.mixing_ratios['CO2']
    l2.proxy_precisions['CH4'] = rel_error * l2.mixing_ratios['CH4']

    log.info('')
    log.info('Writing output')
    write_l2(config['io_files']['l2'], atm, l2, ref_profiles, geometry)

    if config['io_files']['l2_diag']:
        write_l2_diagnostics(config['io_files']['l2_diag'], l1b, l2)

    total_time = time.time() - total_time
    log.info('')
    log.info('Timings:')
    log.info(f"  Optics:             {timings['opt']:7.2f} s")
    log.info(f"  Radiative transfer: {timings['rtm']:7.2f} s")
    log.info(f"  Convolutions:       {timings['conv']:7.2f} s")
    log.info(f"  Derivatives:        {timings['kern']:7.2f} s")
    log.info(f"  Total:              {total_time:7.2f} s")

    print_heading('Success')
