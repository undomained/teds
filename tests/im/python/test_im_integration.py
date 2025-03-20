# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Integration tests for the instrument model."""
from pytest import approx

from teds.im import run_instrument_model
from teds.l1al1b.python.io import read_l1


def test_full_chain(ckd_file, binningtable_file, sgm_file, tmp_path):
    config = {
        'detector': {
            'exposure_time': 0.01724385,
            'nr_coadditions': 10,
        },
        'isrf': {
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
    assert abs(l1.signal).sum() == approx(82386986.0)


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
    assert abs(l1.signal).sum() == approx(8238286.4908721)


def test_full_chain_exact_drawing(
        ckd_file, binningtable_file, sgm_file, tmp_path):
    config = {
        'cal_level': 'raw',
        'detector': {
            'exposure_time': 0.01724385,
            'nr_coadditions': 10,
            'binning_table_id': 1,
        },
        'isrf': {
            'enabled': False,
            'in_memory': True,
        },
        'swath': {
            'exact_drawing': True,
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
    assert abs(l1.signal).sum() == approx(28508076.4154292)


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
    assert abs(l1.signal).sum() == approx(82371824.0)


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
    assert abs(l1.signal).sum() == approx(8236740.085948976)
