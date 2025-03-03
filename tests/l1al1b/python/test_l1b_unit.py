# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Unit tests for the L1A-L1B processor."""
from pytest import approx

import teds.l1al1b.calibration as cal


def test_dark_offset(l1, ckd_dark):
    cal.dark_offset(l1, ckd_dark.offset)
    assert abs(l1.signal).sum() == approx(7600086.3545328)


def test_noise(l1, binning_table, ckd_dark, ckd_noise):
    cal.noise(l1, binning_table.count_table, ckd_noise, ckd_dark.current, 1.0)
    assert abs(l1.noise).sum() == approx(1054444.9909459)


def test_dark_current(l1, ckd_dark):
    cal.dark_current(l1, ckd_dark.current)
    assert abs(l1.signal).sum() == approx(26064709.4484070)


def test_nonlinearity(l1, ckd_pixel_mask, ckd_nonlin):
    cal.nonlinearity(l1, ckd_pixel_mask, ckd_nonlin)
    assert abs(l1.signal).sum() == approx(8261203.0)
    assert abs(l1.noise).sum() == approx(3264.0)


def test_prnu(l1, ckd_pixel_mask, ckd_prnu):
    cal.prnu(l1, ckd_pixel_mask, ckd_prnu.prnu_qe)
    assert abs(l1.signal).sum() == approx(10834706.9656610)
    assert abs(l1.noise).sum() == approx(4563.8899393)


# Stray light with 1 kernel
def test_stray_1(l1, binning_table, ckd_stray):
    ckd_stray.kernels_fft = ckd_stray.kernels_fft[:1, :]
    ckd_stray.weights = ckd_stray.weights[:1, :]
    ckd_stray.edges = ckd_stray.edges[:1, :]
    cal.stray_light(l1, binning_table, ckd_stray, 3)
    assert abs(l1.signal).sum() == approx(9178937.4549845)
    assert abs(l1.noise).sum() == approx(3264.0)


# Stray light with 2 kernels
def test_stray_2(l1, binning_table, ckd_stray):
    cal.stray_light(l1, binning_table, ckd_stray, 3)
    assert abs(l1.signal).sum() == approx(9178760.4724600)
    assert abs(l1.noise).sum() == approx(3264.0)


def test_swath_baseline_mapping(l1, binning_table, ckd_swath, ckd_spectral):
    cal.map_from_detector(l1,
                          ckd_swath,
                          binning_table.count_table,
                          ckd_spectral.wavelengths,
                          False)
    assert abs(l1.spectra).sum() == approx(14858166.2984863)
    assert abs(l1.spectra_noise).sum() == approx(5000.0)


def test_swath_exact_mapping(l1, binning_table, ckd_swath, ckd_spectral):
    cal.map_from_detector(l1,
                          ckd_swath,
                          binning_table.count_table,
                          ckd_spectral.wavelengths,
                          True)
    assert abs(l1.spectra).sum() == approx(8136262.5185032)
    assert abs(l1.spectra_noise).sum() == approx(2936.414486)


def test_radiometric(l1, ckd_radiometric):
    l1.spectra /= 1e12
    l1.spectra_noise /= 1e12
    cal.radiometric(l1, ckd_radiometric.rad_corr)
    assert abs(l1.spectra).sum() == approx(2.6684607e+22)
    assert abs(l1.spectra_noise).sum() == approx(1.0147225e+20)
