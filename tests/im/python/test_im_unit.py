# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Unit tests for the instrument model."""
from pytest import approx
from scipy.interpolate import CubicSpline
import numpy as np

import teds.im.forward_models as fw


def test_dark_offset(l1, ckd_dark):
    fw.dark_offset(l1, ckd_dark.offset)
    assert abs(l1.signal).sum() == approx(8922319.6454672)


def test_noise(l1, binning_table, ckd_dark, ckd_noise):
    fw.noise(l1, ckd_noise, ckd_dark.current, 1, 1.0, 0)
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
    n_alt, n_act, _ = l1.spectra.shape
    spectra_new = np.empty((n_alt, n_act, len(ckd_swath.wavelengths)))
    for i_alt in range(n_alt):
        for i_act in range(n_act):
            spline = CubicSpline(ckd_spectral.wavelengths[i_act, ::-1],
                                 l1.spectra[i_alt, i_act, ::-1])
            spectra_new[i_alt, i_act, :] = spline(ckd_swath.wavelengths)
    l1.spectra = spectra_new
    l1.wavelengths = np.tile(ckd_swath.wavelengths, n_act).reshape(n_act, -1)
    fw.map_to_detector(l1, ckd_swath, ckd_spectral.wavelengths, False)
    assert abs(l1.spectra).sum() == approx(4.5481734e+20)


def test_detector_radiometric(l1, ckd_radiometric):
    fw.radiometric(l1, ckd_radiometric.rad_corr)
    assert abs(l1.spectra).sum() == approx(3178578.5298749)
