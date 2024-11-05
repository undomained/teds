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
from importlib.resources import files
from pathlib import Path

import numpy as np
import yaml

from . import forward_models as fw
from teds import log
from teds.l1al1b.python.io import print_heading
from teds.l1al1b.python.io import print_system_info
from teds.l1al1b.python.io import read_binning_table
from teds.l1al1b.python.io import read_ckd
from teds.l1al1b.python.io import read_l1
from teds.l1al1b.python.io import read_proc_level
from teds.l1al1b.python.io import write_l1
from teds.l1al1b.python.types import ProcLevel
from teds.l1al1b.python.types import L1
from teds.lib.io import merge_configs


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config_file
        Path of YAML configuration file.
    towards_l1b
        Whether the processing direction is to L1B or L1A

    """
    for key in ('sgm', 'ckd'):
        input_file = Path(config['io'][key])
        if not input_file.is_file():
            raise SystemExit(f"ERROR: {input_file} not found")
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
    return in_level >= proc_step and out_level < proc_step


def set_l1_meta(config: dict, l1_product: L1) -> None:
    """Set L1 meta data from settings.

    Meta data such as the exposure time are read from the L1A input
    file by the L1B processor but ignored by the instrument
    model. Instead, they are set by the user.

    Parameters
    ----------
    config
        Configuration dictionary
    l1_product
        L1 product (signal and detector settings).

    """
    l1_product['binning_table_ids'][:] = config['detector']['binning_table_id']
    l1_product['coad_factors'][:] = config['detector']['nr_coadditions']
    l1_product['exptimes'][:] = config['detector']['exposure_time']


def run_instrument_model(config_user: dict | None = None) -> None:
    """Run a simulations of the instrument model.

    Read input, process data, and write output. If run without an
    argument, print the list of all settings and exit.

    Parameters
    ----------
    config_user
        Configuration dictionary directly from file, as given by the
        user, to be expanded with default values for parameters not
        specified by the user.

    """
    defaults_filename = str(files('teds.im.python') / 'default_config.yaml')
    if config_user is None:
        print(open(defaults_filename).read())
        return
    assert isinstance(config_user, dict)

    print_heading('Tango instrument model', empty_line=False)
    print_system_info()

    print_heading('Reading CKD and input data')
    # Start with the full default config and then merge in those
    # settings given by the user.
    config: dict = yaml.safe_load(open(defaults_filename))
    merge_configs(config, config_user)
    check_config(config)
    ckd = read_ckd(config['io']['ckd'])
    log.info('Reading input data')
    l1_product: L1 = read_l1(
        config['io']['sgm'], config['image_start'], config['image_end'])
    set_l1_meta(config, l1_product)
    # Read binning table corresponding to input data
    binning_table = read_binning_table(config['io']['binning_table'],
                                       config['detector']['binning_table_id'],
                                       ckd['n_detector_rows'],
                                       ckd['n_detector_cols'])
    # Processing steps
    print_heading('Forward model')
    # Output processing level
    cal_level = ProcLevel[config['cal_level'].lower()]
    if step_needed(ProcLevel.sgm, l1_product['proc_level'], cal_level):
        log.info('ISRF convolution')
        fw.apply_isrf(l1_product,
                      ckd['swath']['wavelengths'],
                      config['isrf']['enabled'],
                      config['isrf']['fwhm_gauss'],
                      config['isrf']['shape'])
    if config['l1b']['enabled'] and step_needed(
            ProcLevel.l1b, l1_product['proc_level'], cal_level):
        log.info('Radiometric')
        fw.radiometric(l1_product, ckd['radiometric']['rad_corr'])
    if step_needed(ProcLevel.swath, l1_product['proc_level'], cal_level):
        log.info('Detector mapping')
        fw.map_to_detector(l1_product, ckd['swath'])
    if config['stray']['enabled'] and step_needed(
            ProcLevel.stray, l1_product['proc_level'], cal_level):
        log.info('Stray light')
        fw.stray_light(l1_product, ckd['stray'])
    if config['prnu']['enabled'] and step_needed(
            ProcLevel.prnu, l1_product['proc_level'], cal_level):
        log.info('PRNU')
        fw.prnu(l1_product, ckd['prnu']['prnu_qe'])
    if config['nonlin']['enabled'] and step_needed(
            ProcLevel.nonlin, l1_product['proc_level'], cal_level):
        log.info('Nonlinearity')
        fw.nonlinearity(l1_product, ckd['nonlin'])
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark_current, l1_product['proc_level'], cal_level):
        log.info('Dark signal')
        fw.dark_current(l1_product, ckd['dark']['current'])
    if config['noise']['enabled'] and step_needed(
            ProcLevel.noise, l1_product['proc_level'], cal_level):
        log.info('Noise')
        fw.noise(l1_product,
                 ckd['noise'],
                 ckd['dark']['current'],
                 config['noise']['seed'])
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark_offset, l1_product['proc_level'], cal_level):
        log.info('Dark offset')
        fw.dark_offset(l1_product, ckd['dark']['offset'])
    if step_needed(ProcLevel.raw, l1_product['proc_level'], cal_level):
        log.info('Coadding and binning')
        fw.coadding_and_binning(l1_product, binning_table)

    # Write output data
    log.info('Writing output data')
    write_l1(config['io']['l1a'], config, l1_product)

    # If this is shown then the simulation ran successfully
    print_heading('Success')
