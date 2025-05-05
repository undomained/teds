# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Integration tests for the instrument model."""
from pytest import approx

from teds.im import run_instrument_model
from teds.l1al1b.io import read_l1


def test_full_chain(ckd_file, binningtable_file, sgm_file, tmp_path):
    config = {
        'detector': {
            'exposure_time': 0.01724385,
            'nr_coadditions': 10,
        },
        'isrf': {
            'tabulated': False,
            'in_memory': True,
            'fwhm_gauss': 0.5,
        },
        'io_files': {
            'ckd': ckd_file,
            'binning_table': binningtable_file,
            'sgm': sgm_file,
            'l1a': tmp_path / 'tango_l1a.nc',
        },
    }
    run_instrument_model(config)
    l1 = read_l1(tmp_path / 'tango_l1a.nc', 0, None, True)
    assert abs(l1.signal).sum() == approx(82179075.0)


def test_full_chain_no_adc_binning(
        ckd_file, binningtable_file, sgm_file, tmp_path):
    config = {
        'cal_level': 'raw',
        'detector': {
            'exposure_time': 0.01724385,
            'nr_coadditions': 10,
        },
        'isrf': {
            'enabled': False,
            'tabulated': False,
            'in_memory': True,
        },
        'io_files': {
            'ckd': ckd_file,
            'binning_table': binningtable_file,
            'sgm': sgm_file,
            'l1a': tmp_path / 'tango_l1a.nc',
        },
    }
    run_instrument_model(config)
    l1 = read_l1(tmp_path / 'tango_l1a.nc', 0, None, True)
    assert abs(l1.signal).sum() == approx(8217417.1351802)


def test_full_chain_binning_4(
        ckd_file, binningtable_file, sgm_file, tmp_path):
    config = {
        'detector': {
            'exposure_time': 0.01724385,
            'nr_coadditions': 10,
            'binning_table_id': 4,
        },
        'isrf': {
            'enabled': False,
            'tabulated': False,
            'in_memory': True,
        },
        'io_files': {
            'sgm': sgm_file,
            'ckd': ckd_file,
            'binning_table': binningtable_file,
            'l1a': tmp_path / 'tango_l1a.nc',
        },
    }
    run_instrument_model(config)
    l1 = read_l1(tmp_path / 'tango_l1a.nc', 0, None, True)
    assert abs(l1.signal).sum() == approx(82167568.0)


def test_full_chain_binning_4_no_adc(
        ckd_file, binningtable_file, sgm_file, tmp_path):
    config = {
        'cal_level': 'raw',
        'detector': {
            'exposure_time': 0.01724385,
            'nr_coadditions': 10,
            'binning_table_id': 4,
        },
        'isrf': {
            'enabled': False,
            'tabulated': False,
            'in_memory': True,
        },
        'io_files': {
            'ckd': ckd_file,
            'sgm': sgm_file,
            'binning_table': binningtable_file,
            'l1a': tmp_path / 'tango_l1a.nc',
        },
    }
    run_instrument_model(config)
    l1 = read_l1(tmp_path / 'tango_l1a.nc', 0, None, True)
    assert abs(l1.signal).sum() == approx(8216287.8303390)
