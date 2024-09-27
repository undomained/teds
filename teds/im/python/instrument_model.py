# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Instrument model (level 1B to 1A processor).

Level 1A data are unprocessed, but annotated instrument data. It contains
detector artifacts such as dark current and instrument artifacts such as stray
light. Level 1B data have been processed to sensor units, in this case radiance
units. Data at an intermediate level can also be read and produced.

This is a Python version of the corresponding C++ code in the TANGO
end-to-end simulator (tango_im.x).

Input files are:
- YAML configuration file with at least keys io.l1a, io.l1b, and io.ckd
- netCDF input data (io.l1a)
- netCDF calibration key data io.ckd
- (optionally) netCDF binning-table set io.binning_table
- (optionally) netCDF geometry data io.geometry

"""
from pathlib import Path
import argparse
import os

import numpy as np
import yaml

from . import uncalibration as cal
from teds import log
from teds.l1al1b.io import print_heading
from teds.l1al1b.io import print_system_info
from teds.l1al1b.io import read_binning_pattern
from teds.l1al1b.io import read_ckd
from teds.l1al1b.io import read_l1
from teds.l1al1b.io import read_proc_level
from teds.l1al1b.io import write_l1
from teds.l1al1b.types import ProcLevel
from teds.l1al1b.types import L1
from teds.lib.utils import merge_configs


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Args:
      config_file:
        Path of YAML configuration file.
      towards_l1b:
        Whether the processing direction is to L1B or L1A

    """
    for key in ('sgm', 'ckd'):
        input_file = Path(config['io'][key])
        if not input_file.is_file():
            raise SystemExit(f"ERROR: {input_file} not found")
    if not os.access(config['io']['l1a'], os.W_OK):
        raise SystemExit(f"ERROR: {config['io']['l1a']} not writable")
    proc_level = read_proc_level(config['io']['sgm'])
    log.info(f"Processing from {proc_level} to {config['cal_level']}")


def step_needed(proc_step: ProcLevel,
                in_level: ProcLevel,
                out_level: ProcLevel) -> bool:
    """Test whether a processing step is needed.

    Needed processing steps are between in_level (inclusive) and
    out_level (inclusive).

    This is similar to the function in l1b_processor.py but the
    process chain is traversed in the opposite direction.

    Args:
      proc_step:
        Processing step to be tested.
      in_level:
        Processing level of the input data.
      out_level:
        Processing level of the output data.

    Returns:
      Whether the step is needed.

    """
    return in_level >= proc_step and out_level < proc_step


def process_im(config_user: dict) -> None:
    """Run a simulations of the instrument model.

    Read input, process data, and write output.

    Args:
      config_user:
        Configuration dictionary directly from file, as given by the user, to
        be expanded with default values for parameters not specified by the
        user.

    """
    print_heading('Tango instrument model', empty_line=False)
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
        config['io']['sgm'], config['image_start'], config['image_end'])

    # Read binning table corresponding to input data
    log.info('Reading binning table')
    data_binning_table_id = int(l1_products['binning_table_ids'][0])
    # Last test skips when data is from C++ code
    if (
            l1_products['proc_level'] != 'l1a'
            and data_binning_table_id > 0
            and l1_products['original_file'] != ''):
        raise SystemExit("ERROR: binned data cannot be processed towards L1A")
    # C++ code reads the config binning table instead of that
    # specified in the data the config binning table is only needed
    # for the conversion from RAW to L1A data.
    bin_indices, count_table = read_binning_pattern(
        config['io']['binning_table'],
        data_binning_table_id,
        ckd['n_detrows'],
        ckd['n_detcols'])
    # Processing steps
    print_heading('Forward model')
    # Output processing level
    cal_level = ProcLevel[config['cal_level'].lower()]
    if step_needed(ProcLevel.sgm, l1_products['proc_level'], cal_level):
        # in C++ code also performed if cal_level = 'l1b' or 'sgm'
        log.info('Changing wavelength grid')
        cal.change_wavelength_grid(l1_products,
                                   ckd['spectral']['wavelengths'],
                                   config['isrf']['enabled'],
                                   config['isrf']['fwhm_gauss'],
                                   config['isrf']['shape'])
    if config['l1b']['enabled'] and step_needed(
            ProcLevel.l1b, l1_products['proc_level'], cal_level):
        log.info('Converting from radiance')
        cal.convert_from_radiance(l1_products,
                                  ckd['radiometric']['rad_corr'],
                                  config['detector']['exposure_time'])
    if step_needed(ProcLevel.swath, l1_products['proc_level'], cal_level):
        log.info('Mapping to detector')
        cal.map_to_detector(
            l1_products, ckd['swath']['spectrum_rows'], ckd['n_detrows'])
    if config['stray']['enabled'] and step_needed(
            ProcLevel.stray, l1_products['proc_level'], cal_level):
        log.info('Including stray light')
        cal.stray_light(l1_products, ckd['stray'])
    if step_needed(
            ProcLevel.prnu, l1_products['proc_level'], cal_level):
        # Apply pixel mask just before the first pixel-dependent
        # step towards L1A. C++ code keeps bad signals and does
        # not process them (mostly). C++ code pixel mask has a
        # different format than in CKD but is not written, so not
        # helpful.
        log.info('Including pixel mask')
        l1_products['signal'][:, ckd['pixel_mask'].ravel()] = np.nan
        l1_products['noise'][:, ckd['pixel_mask'].ravel()] = np.nan
    if config['prnu']['enabled'] and step_needed(
            ProcLevel.prnu, l1_products['proc_level'], cal_level):
        log.info('Including PRNU')
        cal.include_prnu(l1_products, ckd['prnu']['prnu_qe'])
    if config['nonlin']['enabled'] and step_needed(
            ProcLevel.nonlin, l1_products['proc_level'], cal_level):
        log.info('Including non-linearity')
        cal.include_nonlinearity(l1_products, ckd['nonlin'])
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark, l1_products['proc_level'], cal_level):
        log.info('Including dark signal')
        cal.include_darksignal(l1_products, ckd['dark']['current'])
    if config['noise']['enabled'] and step_needed(
            ProcLevel.noise, l1_products['proc_level'], cal_level):
        log.info('Including noise')
        cal.include_noise(l1_products,
                          ckd['noise'],
                          ckd['dark']['current'],
                          config['noise']['seed'])
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark, l1_products['proc_level'], cal_level):
        log.info('Including offset')
        cal.include_offset(l1_products, ckd['dark']['offset'])
    if step_needed(ProcLevel.raw, l1_products['proc_level'], cal_level):
        log.info('Including binning')
        cal.include_coadding_and_binning(
            l1_products,
            config['detector']['nr_coadditions'],
            config['io']['binning_table'],
            config['detector']['binning_table_id'],
            ckd['n_detrows'],
            ckd['n_detcols'])

    # Write output data
    log.info('Writing output data')
    write_l1(config['io']['l1a'], config, l1_products, geometry=True)

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
    process_im(config)
