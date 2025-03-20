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
import math

from . import forward_models as fw
from teds import log
from teds.l1al1b.python.io import copy_navigation_data
from teds.l1al1b.python.io import read_binning_table
from teds.l1al1b.python.io import read_ckd
from teds.l1al1b.python.io import read_l1
from teds.l1al1b.python.io import read_proc_level
from teds.l1al1b.python.io import write_l1
from teds.l1al1b.python.types import L1
from teds.l1al1b.python.types import ProcLevel
from teds.lib.io import merge_config_with_default
from teds.lib.io import print_heading
from teds.lib.io import print_system_info


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config
        Configuration parameters.

    """
    # If swath.exact_drawing is true then do not bin the detector
    # image but instead artifically scale noise. Binning table ID is
    # always set to 1 in that case. The value is determined by binning
    # and the detector mapping algorithm choice and is not a user
    # parameter.
    if config['swath']['exact_drawing']:
        config['noise']['artificial_scaling'] = 1 / math.sqrt(
            config['detector']['binning_table_id'])
        config['detector']['binning_table_id'] = 1
    else:
        config['noise']['artificial_scaling'] = 1
    for key in ('sgm', 'ckd'):
        input_file = Path(config['io_files'][key])
        if not input_file.is_file():
            raise SystemExit(f"ERROR: {input_file} not found")
    proc_level = read_proc_level(config['io_files']['sgm'])
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
    l1_product.binning_table_id = config['detector']['binning_table_id']
    l1_product.exposure_time = config['detector']['exposure_time']


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
    print_heading('Tango instrument model', empty_line=False)
    print_system_info()

    print_heading('Reading CKD and input data')
    # Start with the full default config and then merge in those
    # settings given by the user.
    config = merge_config_with_default(config_user, 'teds.im.python')
    check_config(config)
    ckd = read_ckd(config['io_files']['ckd'])
    log.info('Reading input data')
    l1_product: L1 = read_l1(config['io_files']['sgm'],
                             config['alt_beg'],
                             config['alt_end'],
                             config['isrf']['in_memory'],
                             config['io_files']['geometry'])
    set_l1_meta(config, l1_product)
    if len(l1_product.spectra.shape) >= 2 and (
            l1_product.spectra.shape[1] != len(ckd.swath.act_angles)):
        log.error(
            f'L1 product has {l1_product.spectra.shape[1]} ACT positions '
            f'whereas the CKD expects {len(ckd.swath.act_angles)}')
        exit(1)
    # Read binning table corresponding to input data
    binning_table = read_binning_table(config['io_files']['binning_table'],
                                       config['detector']['binning_table_id'],
                                       ckd.n_detector_rows,
                                       ckd.n_detector_cols)
    # Processing steps
    print_heading('Forward model')
    # Output processing level
    cal_level = ProcLevel[config['cal_level'].lower()]
    if step_needed(ProcLevel.sgm, l1_product.proc_level, cal_level):
        log.info('ISRF convolution')
        fw.apply_isrf(l1_product,
                      ckd.swath.wavelengths,
                      config['isrf']['enabled'],
                      config['isrf']['fwhm_gauss'],
                      config['isrf']['shape'],
                      config['io_files']['sgm'],
                      config['alt_beg'])
    if config['l1b']['enabled'] and step_needed(
            ProcLevel.l1b, l1_product.proc_level, cal_level):
        log.info('Radiometric')
        fw.radiometric(l1_product, ckd.radiometric.rad_corr)
    if step_needed(ProcLevel.swath, l1_product.proc_level, cal_level):
        log.info('Detector mapping')
        fw.map_to_detector(l1_product,
                           ckd.swath,
                           ckd.spectral.wavelengths,
                           config['swath']['exact_drawing'])
    if config['stray']['enabled'] and step_needed(
            ProcLevel.stray, l1_product.proc_level, cal_level):
        log.info('Stray light')
        fw.stray_light(l1_product, ckd.stray)
    if config['prnu']['enabled'] and step_needed(
            ProcLevel.prnu, l1_product.proc_level, cal_level):
        log.info('PRNU')
        fw.prnu(l1_product, ckd.prnu.prnu_qe)
    if config['nonlin']['enabled'] and step_needed(
            ProcLevel.nonlin, l1_product.proc_level, cal_level):
        log.info('Nonlinearity')
        fw.nonlinearity(l1_product, ckd.nonlin)
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark_current, l1_product.proc_level, cal_level):
        log.info('Dark signal')
        fw.dark_current(l1_product, ckd.dark.current)
    if config['noise']['enabled'] and step_needed(
            ProcLevel.noise, l1_product.proc_level, cal_level):
        log.info('Noise')
        # Only scale noise by coaddition factor if target level is
        # L1A. Otherwise the signal will not be scaled by coadditions
        # and neither should the noise be scaled.
        if cal_level == ProcLevel.l1a:
            n_coadditions = config['detector']['nr_coadditions']
        else:
            n_coadditions = 1
        fw.noise(l1_product,
                 ckd.noise,
                 ckd.dark.current,
                 n_coadditions,
                 config['noise']['artificial_scaling'],
                 config['noise']['seed'])
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark_offset, l1_product.proc_level, cal_level):
        log.info('Dark offset')
        fw.dark_offset(l1_product, ckd.dark.offset)
    if l1_product.signal.size > 0:
        log.info('Detector image binning '
                 f'({config["detector"]["binning_table_id"]}x1)')
        fw.bin_detector_images(
            l1_product, config['detector']['binning_table_id'], binning_table)
    if step_needed(ProcLevel.raw, l1_product.proc_level, cal_level):
        log.info('Coadding and analog-to-digital conversion')
        fw.coadd_and_adc(l1_product, config['detector']['nr_coadditions'])

    # Write output data
    log.info('Writing output')
    write_l1(config['io_files']['l1a'], config, l1_product)
    if config['io_files']['navigation']:
        copy_navigation_data(config['io_files']['navigation'],
                             config['io_files']['l1a'])

    # If this is shown then the simulation ran successfully
    print_heading('Success')
