# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Tests for the L1A-L1B processor."""
from pytest import approx

import teds.l1al1b.calibration as cal


def test_dark_offset(l1, ckd_dark):
    cal.dark_offset(l1, ckd_dark.offset)
    assert abs(l1.signal).sum() == approx(7600086.354532824)


def test_noise(l1, binning_table, ckd_dark, ckd_noise):
    cal.noise(l1, binning_table.count_table, ckd_noise, ckd_dark.current, 1.0)
    assert abs(l1.noise).sum() == approx(1054444.9909459099)


def test_dark_current(l1, ckd_dark):
    cal.dark_current(l1, ckd_dark.current)
    assert abs(l1.signal).sum() == approx(26064709.448407017)


def test_nonlinearity(l1, ckd_pixel_mask, ckd_nonlin):
    cal.nonlinearity(l1, ckd_pixel_mask, ckd_nonlin)
    assert abs(l1.signal).sum() == approx(8261203.0)
    assert abs(l1.noise).sum() == approx(3264.0)


def test_prnu(l1, ckd_pixel_mask, ckd_prnu):
    cal.prnu(l1, ckd_pixel_mask, ckd_prnu.prnu_qe)
    assert abs(l1.signal).sum() == approx(10834706.965660967)
    assert abs(l1.noise).sum() == approx(4563.8899392561925)
