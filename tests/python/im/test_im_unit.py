# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Unit tests for the instrument model."""
from pytest import approx

import teds.im.forward_models as fw


def test_dark_offset(l1, ckd_dark):
    fw.dark_offset(l1, ckd_dark.offset)
    assert abs(l1.signal).sum() == approx(8922319.6454672)


def test_noise(l1, binning_table, ckd_dark, ckd_noise):
    fw.noise(l1, ckd_noise, ckd_dark.current, 1, 0)
    assert abs(l1.signal).sum() == approx(8229957.1059201)


def test_dark_current(l1, ckd_dark):
    fw.dark_current(l1, ckd_dark.offset)
    assert abs(l1.signal).sum() == approx(8291669.4367081)


def test_nonlin(l1, ckd_nonlin):
    fw.nonlinearity(l1, ckd_nonlin)
    assert abs(l1.signal).sum() == approx(8261203.0)


def test_prnu(l1, ckd_prnu):
    fw.prnu(l1, ckd_prnu.prnu_qe)
    assert abs(l1.signal).sum() == approx(6175154.6307808)


def test_stray(l1, ckd_stray):
    fw.stray_light(l1, ckd_stray)
    assert abs(l1.signal).sum() == approx(7435369.4285431)


def test_detector_mapping(l1, ckd_swath, ckd_spectral):
    fw.map_to_detector(l1, ckd_swath, ckd_spectral.wavelengths, 5)
    assert abs(l1.spectra).sum() == approx(4.3375163e+20)


def test_detector_radiometric(l1, ckd_radiometric):
    fw.radiometric(l1, ckd_radiometric.rad_corr)
    assert abs(l1.spectra).sum() == approx(4733988.3035614)
