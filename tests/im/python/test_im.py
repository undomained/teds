# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Tests for the instrument model."""
from pytest import approx

import teds.im.forward_models as cal


def test_dark_offset(l1, ckd_dark):
    cal.dark_offset(l1, ckd_dark.offset)
    assert abs(l1.signal).sum() == approx(8922319.645467177)


def test_noise(l1, binning_table, ckd_dark, ckd_noise):
    cal.noise(l1, ckd_noise, ckd_dark.current, 1.0, 0)
    assert abs(l1.signal).sum() == approx(8229957.105920106)
