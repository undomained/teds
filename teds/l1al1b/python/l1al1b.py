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
import math
import numpy as np

from . import calibration as cal
from .binning import bin_data
from .io import copy_geometry
from .io import read_binning_table
from .io import read_ckd
from .io import read_l1
from .io import read_proc_level
from .io import write_l1
from .types import L1
from .types import ProcLevel
from teds import log
from teds.lib.io import merge_config_with_default
from teds.lib.io import print_heading
from teds.lib.io import print_system_info


def check_config(config: dict) -> None:
    """Check consistency of some of the configuration settings.

    Parameters
    ----------
    config
        Configuration dictionary

    """
    # If swath.exact_drawing is true then do not bin the L1B data but
    # instead artifically scale noise. bin_spectra is always set to 1
    # in that case. The value is determined by binning and the
    # detector mapping algorithm choice and is not a user parameter.
    config['noise']['artificial_scaling'] = (
        1 / math.sqrt(config['bin_spectra']))
    if config['swath']['exact_drawing']:
        config['bin_spectra'] = 1
    for key in ('l1a', 'ckd'):
        input_file = Path(config['io'][key])
        if not input_file.is_file():
            raise SystemExit(f"ERROR: {input_file} not found")
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


def run_l1al1b(config_user: dict | None = None) -> None:
    """Run a simulation of the L1A-L1B process.

    Read input, process data, and write output. If run without an
    argument, print the list of all settings and exit.

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
    config = merge_config_with_default(config_user, 'teds.l1al1b.python')
    check_config(config)
    ckd = read_ckd(config['io']['ckd'])
    log.info('Reading input data')
    l1_product: L1 = read_l1(
        config['io']['l1a'], config['alt_beg'], config['alt_end'])

    # Read binning table corresponding to input data
    binning_table = read_binning_table(config['io']['binning_table'],
                                       l1_product.binning_table_id,
                                       ckd.n_detector_rows,
                                       ckd.n_detector_cols)
    # Bin the CKD
    ckd.pixel_mask = bin_data(binning_table, ckd.pixel_mask)
    ckd.dark.offset = bin_data(binning_table, ckd.dark.offset)
    ckd.dark.current = bin_data(binning_table, ckd.dark.current)
    ckd.noise.conversion_gain = bin_data(
        binning_table, ckd.noise.conversion_gain)
    ckd.noise.read_noise = bin_data(binning_table, ckd.noise.read_noise)
    ckd.prnu.prnu_qe = bin_data(binning_table, ckd.prnu.prnu_qe)
    # Processing steps
    print_heading('Retrieval')
    # Output processing level
    cal_level = ProcLevel[config['cal_level'].lower()]
    if step_needed(ProcLevel.raw, l1_product.proc_level, cal_level):
        log.info('Coadding and binning')
        cal.coadding_and_binning(l1_product, binning_table.count_table)
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark_offset, l1_product.proc_level, cal_level):
        log.info('Dark offset')
        cal.dark_offset(l1_product, ckd.dark.offset)
    if step_needed(ProcLevel.noise, l1_product.proc_level, cal_level):
        log.info('Noise')
        if config['dark']['enabled'] and config['noise']['enabled']:
            cal.noise(l1_product,
                      binning_table.count_table,
                      ckd.noise,
                      ckd.dark.current,
                      config['noise']['artificial_scaling'])
        else:
            l1_product.noise = np.full_like(l1_product.signal, 1.0)
    if config['dark']['enabled'] and step_needed(
            ProcLevel.dark_current, l1_product.proc_level, cal_level):
        log.info('Dark current')
        cal.dark_current(l1_product, ckd.dark.current)
    if config['nonlin']['enabled'] and step_needed(
            ProcLevel.nonlin, l1_product.proc_level, cal_level):
        log.info('Nonlinearity')
        cal.nonlinearity(l1_product, ckd.pixel_mask, ckd.nonlin)
    if config['prnu']['enabled'] and step_needed(
            ProcLevel.prnu, l1_product.proc_level, cal_level):
        log.info('PRNU')
        cal.prnu(l1_product, ckd.pixel_mask, ckd.prnu.prnu_qe)
    if l1_product.signal.size > 0:
        log.info('Smoothing over bad values')
        cal.remove_bad_values(
            ckd.n_detector_cols, ckd.pixel_mask, l1_product.signal)
        cal.remove_bad_values(
            ckd.n_detector_cols, ckd.pixel_mask, l1_product.noise)
    if (
            step_needed(ProcLevel.stray, l1_product.proc_level, cal_level)
            and config['stray']['van_cittert_steps'] > 0):
        log.info('Stray light')
        cal.stray_light(l1_product,
                        binning_table,
                        ckd.stray,
                        config['stray']['van_cittert_steps'])
    if step_needed(ProcLevel.swath, l1_product.proc_level, cal_level):
        log.info('Mapping from detector')
        cal.map_from_detector(l1_product,
                              ckd.swath,
                              binning_table.count_table,
                              ckd.spectral.wavelengths,
                              config['swath']['exact_drawing'])
    if l1_product.spectra.shape[2] != ckd.spectral.wavelengths.shape[1]:
        log.info(
            'Interpolating from intermediate to main CKD wavelength grids')
        cal.change_wavelength_grid(l1_product, ckd.spectral.wavelengths)
    if step_needed(ProcLevel.l1b, l1_product.proc_level, cal_level):
        log.info('Radiometric')
        cal.radiometric(l1_product, ckd.radiometric.rad_corr)
    copy_geometry(config['io']['l1a'],
                  config['io']['geometry'],
                  config['alt_beg'],
                  l1_product)
    if l1_product.proc_level >= ProcLevel.swath:
        cal.bin_l1b(l1_product, config['bin_spectra'])
    log.info('Writing output data')
    write_l1(config['io']['l1b'], config, l1_product)

    # If this is shown then the simulation ran successfully
    print_heading('Success')
