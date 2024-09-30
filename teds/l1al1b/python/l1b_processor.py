# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Level 1A to 1B processor.

Level 1A data are unprocessed, but annotated instrument data. It contains
detector artifacts such as dark current and instrument artifacts such as stray
light. Level 1B data have been processed to sensor units, in this case radiance
units. Data at an intermediate level can also be read and produced.

This is a Python version of the corresponding C++ code in the TANGO
end-to-end simulator (tango_l1b.x).

Input files are:
- YAML configuration file with at least the keys io.l1a, io.sgm, and io.ckd
- netCDF input data (io.sgm)
- netCDF calibration key data io.ckd
- (optionally) netCDF binning-table set io.binning_table
- (optionally) netCDF geometry data io.geometry

"""
from pathlib import Path
import argparse
import os

import numpy as np
import yaml

from . import calibration as cal
from .binning import bin_data
from .io import print_heading
from .io import print_system_info
from .io import read_binning_pattern
from .io import read_ckd
from .io import read_l1
from .io import read_proc_level
from .io import write_l1
from .types import ProcLevel
from .types import L1
from teds import log
from teds.lib.utils import merge_configs


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config
        Configuration dictionary
    towards_l1b
        Whether the processing direction is to L1B or L1A

    """
    for key in ('l1a', 'ckd'):
        input_file = Path(config['io'][key])
        if not input_file.is_file():
            raise SystemExit(f"ERROR: {input_file} not found")
    if not os.access(config['io']['l1b'], os.W_OK):
        raise SystemExit(f"ERROR: {config['io']['l1b']} not writable")
    proc_level = read_proc_level(config['io']['l1a'])
    log.info(f"Processing from {proc_level} to {config['cal_level']}")


def step_needed(proc_step: ProcLevel,
                in_level: ProcLevel,
                out_level: ProcLevel) -> bool:
    """Test whether a processing step is needed.

    Needed processing steps are between in_level (inclusive) and
    out_level (inclusive).

    Parameters
    ----------
    proc_step
        Processing step to be tested.
    in_level
        Processing level of the input data.
    out_level
        Processing level of the output data.

    Returns
    -------
        Whether the step is needed.

    """
    return in_level < proc_step and out_level >= proc_step


def process_l1b(config_user: dict) -> None:
    """Run a simulation of the L1A-L1B process.

    Read input, process data, and write output.

    Parameters
    ----------
    config_user
        Configuration dictionary directly from file, as given by the user, to
        be expanded with default values for parameters not specified by the
        user.

    """
    print_heading('Tango L1B processor', empty_line=False)
    print_system_info()

    print_heading('Reading CKD and input data')
    # Start with the full default config and then merge in those
    # settings given by the user.
    defaults_filename = Path(__file__).parent / 'default_config.yaml'
    config: dict = yaml.safe_load(open(defaults_filename))
    merge_configs(config, config_user)
    check_config(config)
    # Read input data
    ckd = read_ckd(config['io']['ckd'])

    # Read input data
    log.info('Reading input data')
    l1_products: L1 = read_l1(
        config['io']['l1a'], config['image_start'], config['image_end'])

    # Read binning table corresponding to input data
    log.info('Reading binning table')
    data_binning_table_id = int(l1_products['binning_table_ids'][0])
    bin_indices, count_table = read_binning_pattern(
        config['io']['binning_table'],
        data_binning_table_id,
        ckd['n_detrows'],
        ckd['n_detcols'])
    # Binning CKD instead of unbinning data
    if config['unbinning'] == 'none' and not (
            bin_indices.ravel()
            == np.arange(ckd['n_detrows'] * ckd['n_detcols'])).all():
        ckd['pixel_mask'] = bin_data(
            np.array([ckd['pixel_mask']]), bin_indices, count_table)[0]
        if np.ndim(ckd['noise']['conversion_gain']) == 0:
            (ckd['dark']['offset'],
             ckd['dark']['current'],
             squared_read_noise,
             ckd['prnu']['prnu_qe']) = bin_data(
                 np.array([ckd['dark']['offset'],
                           ckd['dark']['current'],
                           ckd['noise']['read_noise']**2,
                           ckd['prnu']['prnu_qe']]),
                 bin_indices,
                 count_table)/count_table
        else:
            (ckd['dark']['offset'],
             ckd['dark']['current'],
             ckd['noise']['conversion_gain'],
             squared_read_noise,
             ckd['prnu']['prnu_qe']) = bin_data(
                 np.array([ckd['dark']['offset'],
                           ckd['dark']['current'],
                           ckd['noise']['conversion_gain'],
                           ckd['noise']['read_noise']**2,
                           ckd['prnu']['prnu_qe']]),
                 bin_indices,
                 count_table)/count_table
        ckd['noise']['read_noise'] = np.sqrt(squared_read_noise)
    # Processing steps
    print_heading('Retrieval')
    # Output processing level
    cal_level = ProcLevel[config['cal_level'].lower()]
    if step_needed(ProcLevel.raw, l1_products['proc_level'], cal_level):
        log.info('Removing binning')
        cal.remove_coadding_and_binning(
            l1_products, bin_indices, count_table, config['unbinning'])
        # Apply pixel mask just before the first pixel-dependent
        # step towards L1B. C++ code keeps bad signals except when
        # saving L1B data and does not process them.
        log.info('Including pixel mask')
        l1_products['signal'][:, ckd['pixel_mask'].ravel()] = np.nan
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark, l1_products['proc_level'], cal_level):
        log.info('Removing offset')
        cal.remove_offset(l1_products, ckd['dark']['offset'])
    # C++ code checks 'noise'
    if step_needed(ProcLevel.dark, l1_products['proc_level'], cal_level):
        log.info('Determining noise')
        # C++ code does not have first check
        if config['dark']['enabled'] and config['noise']['enabled']:
            cal.determine_noise(
                l1_products, ckd['noise'], ckd['dark']['current'])
        else:
            l1_products['noise'] = np.full_like(l1_products['signal'], np.nan)
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark, l1_products['proc_level'], cal_level):
        log.info('Removing dark signal')
        cal.remove_dark_signal(l1_products, ckd['dark']['current'])
    if config['nonlin']['enabled'] and step_needed(
            ProcLevel.nonlin, l1_products['proc_level'], cal_level):
        log.info('Removing non-linearity')
        cal.remove_nonlinearity(l1_products, ckd['nonlin'])
    if config['prnu']['enabled'] and step_needed(
            ProcLevel.prnu, l1_products['proc_level'], cal_level):
        log.info('Removing PRNU')
        cal.remove_prnu(l1_products, ckd['prnu']['prnu_qe'])
    if step_needed(
            ProcLevel.stray, l1_products['proc_level'], cal_level):
        log.info('Removing stray light')
        cal.stray_light(l1_products, ckd['stray'])
    if step_needed(
            ProcLevel.swath, l1_products['proc_level'], cal_level):
        log.info('Mapping from detector')
        cal.map_from_detector(l1_products,
                              ckd['spectral']['wavelengths'],
                              ckd['swath']['spectrum_rows'],
                              config['swath']['spectrum_width'],
                              ckd['pixel_mask'],
                              bin_indices)
    if step_needed(ProcLevel.l1b, l1_products['proc_level'], cal_level):
        # C++ code keeps this step optional with config['l1b']['enabled']
        log.info('Converting to radiance')
        # exposure time should never be needed explicitly
        cal.convert_to_radiance(l1_products,
                                ckd['radiometric']['rad_corr'],
                                l1_products['exptimes'])
    # C++ code does not have second test
    if (
            config['reverse_wavelength']
            and l1_products['proc_level'] > ProcLevel.stray):
        l1_products['wavelengths'] = np.flip(l1_products['wavelengths'],
                                             axis=-1)
        l1_products['spectra'] = np.flip(l1_products['spectra'], axis=-1)
        l1_products['spectra_noise'] = np.flip(l1_products['spectra_noise'],
                                               axis=-1)

    # Write output data
    log.info('Writing output data')
    write_l1(config['io']['l1b'], config, l1_products, geometry=True)

    # If this is shown then the simulation ran successfully
    print_heading('Success')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Level 1 processing: conversion between detector (L1A) "
        "and radiance (L1B) data.")
    parser.add_argument(
        'config_file', help="path of the YAML file with the configuration")
    args = parser.parse_args()
    config = yaml.safe_load(open(args.config_file))
    process_l1b(config)
