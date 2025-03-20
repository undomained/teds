# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Integration tests for the L1A-L1B processor."""
from pytest import approx
import pytest

from teds.gm.io import read_geometry
from teds.l1al1b import run_l1al1b
from teds.l1al1b.io import read_l1


def test_full_chain_no_binning(
        l1a_file, ckd_file, binningtable_file, tmp_path):
    config = {
        'io_files': {
            'binning_table': binningtable_file,
            'l1a': l1a_file,
            'ckd': ckd_file,
            'l1b': tmp_path / 'tango_l1b.nc',
        },
    }
    run_l1al1b(config)
    l1 = read_l1(tmp_path / 'tango_l1b.nc', 0, None, True)
    assert abs(l1.spectra).sum() == approx(1.0760962e+21)
    assert abs(l1.spectra_noise).sum() == approx(4.4683595e+18)


@pytest.mark.parametrize('l1a_file', [{'bin': 2}], indirect=True)
def test_full_chain_l1a_binning(
        l1a_file, ckd_file, binningtable_file, tmp_path):
    config = {
        'io_files': {
            'binning_table': binningtable_file,
            'l1a': l1a_file,
            'ckd': ckd_file,
            'l1b': tmp_path / 'tango_l1b.nc',
        },
    }
    run_l1al1b(config)
    l1 = read_l1(tmp_path / 'tango_l1b.nc', 0, None, True)
    assert abs(l1.spectra).sum() == approx(1.0867994e+21)
    assert abs(l1.spectra_noise).sum() == approx(3.1926969e+18)


def test_full_chain_l1b_binning(
        l1a_file, ckd_file, binningtable_file, tmp_path):
    config = {
        'bin_spectra': 5,
        'io_files': {
            'binning_table': binningtable_file,
            'l1a': l1a_file,
            'ckd': ckd_file,
            'l1b': tmp_path / 'tango_l1b.nc',
        },
    }
    run_l1al1b(config)
    l1 = read_l1(tmp_path / 'tango_l1b.nc', 0, None, True)
    assert abs(l1.spectra).sum() == approx(2.1521924e+20)
    assert abs(l1.spectra_noise).sum() == approx(3.9966223e+17)


def test_full_chain_exact_mapping(
        l1a_file, ckd_file, binningtable_file, tmp_path):
    config = {
        'swath': {
            'exact_drawing': True,
        },
        'io_files': {
            'binning_table': binningtable_file,
            'l1a': l1a_file,
            'ckd': ckd_file,
            'l1b': tmp_path / 'tango_l1b.nc',
        },
    }
    run_l1al1b(config)
    l1 = read_l1(tmp_path / 'tango_l1b.nc', 0, None, True)
    assert abs(l1.spectra).sum() == approx(1.0795310e+21)
    assert abs(l1.spectra_noise).sum() == approx(4.1243568e+18)


def test_geolocation(
        l1a_file, ckd_file, binningtable_file, tmp_path):
    config = {
        'cal_level': 'raw',
        'io_files': {
            'binning_table': binningtable_file,
            'l1a': l1a_file,
            'ckd': ckd_file,
            'l1b': tmp_path / 'tango_l1b.nc',
        },
    }
    run_l1al1b(config)
    geo = read_geometry(tmp_path / 'tango_l1b.nc')
    assert abs(geo.lat).sum() == approx(2590.0187941)
    assert abs(geo.lon).sum() == approx(722.7365456)
    assert abs(geo.height).sum() == approx(0.0)
    assert abs(geo.vza).sum() == approx(46.5414420)
    assert abs(geo.vaa).sum() == approx(4495.8432605)
    assert abs(geo.sza).sum() == approx(1959.7164234)
    assert abs(geo.saa).sum() == approx(7716.1682456)
