# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Tests for the L1A-L1B processor."""
from pytest import approx

import teds.l1al1b.calibration as cal


def test_dark_offset(l1, ckd_dark):
    cal.remove_offset(l1, ckd_dark['offset'])
    assert abs(l1['signal']).sum() == approx(7600086.354532824)


def test_noise(l1, ckd_dark, ckd_noise):
    cal.determine_noise(l1, ckd_noise, ckd_dark['current'])
    assert abs(l1['noise']).sum() == approx(3441708450.447448)


def test_dark_current(l1, ckd_dark):
    cal.remove_dark_signal(l1, ckd_dark['current'])
    assert abs(l1['signal']).sum() == approx(26064709.448407017)


def test_nonlinearity(l1, ckd_nonlin):
    cal.remove_nonlinearity(l1, ckd_nonlin)
    assert abs(l1['signal']).sum() == approx(8261203.0)
    assert abs(l1['noise']).sum() == approx(3198.5)


def test_prnu(l1, ckd_prnu):
    cal.remove_prnu(l1, ckd_prnu['prnu_qe'])
    assert abs(l1['signal']).sum() == approx(11443847.197698373)
    assert abs(l1['noise']).sum() == approx(4624.803962459934)
