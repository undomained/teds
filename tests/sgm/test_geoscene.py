# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Tests for the scene generation module."""
from pytest import approx
import numpy as np
import numpy.typing as npt

from teds.sgm.atmosphere import Atmosphere
from teds.sgm.atmosphere import get_AFGL_atm_homogenous_distribution


def get_gas(atm: Atmosphere, gas_name: str) -> npt.NDArray[np.float64]:
    return atm.get_gas(gas_name).concentration


def test_atmosphere(atmosphere_path):
    nlay = 20
    dzlay = 1000
    psurf = 101300
    nlev = nlay + 1
    zlay = (np.arange(nlay - 1, -1, -1) + 0.5) * dzlay
    zlev = np.arange(nlev - 1, -1, -1) * dzlay
    atmosphere = Atmosphere.from_file(zlay, zlev, psurf, atmosphere_path)
    assert abs(get_gas(atmosphere, 'O3')).sum() == approx(5.8943662e+22)
    assert abs(get_gas(atmosphere, 'H2O')).sum() == approx(4.7657731e+26)
    assert abs(get_gas(atmosphere, 'CO2')).sum() == approx(7.0316294e+25)
    assert abs(get_gas(atmosphere, 'NO2')).sum() == approx(2.6392520e+19)
    assert abs(get_gas(atmosphere, 'O2')).sum() == approx(4.4608026e+28)
    assert abs(get_gas(atmosphere, 'CH4')).sum() == approx(3.8326705e+23)
    psurf = 1e5
    atmosphere = Atmosphere.from_file(zlay, zlev, psurf, atmosphere_path)
    assert abs(get_gas(atmosphere, 'O3')).sum() == approx(5.6443625e+21)
    assert abs(get_gas(atmosphere, 'H2O')).sum() == approx(1.643443e+27)
    assert abs(get_gas(atmosphere, 'CO2')).sum() == approx(7.002405e+25)
    assert abs(get_gas(atmosphere, 'NO2')).sum() == approx(4.880464e+18)
    assert abs(get_gas(atmosphere, 'O2')).sum() == approx(4.441659e+28)
    assert abs(get_gas(atmosphere, 'CH4')).sum() == approx(3.816223e+23)
    psurf = 2e5
    atmosphere = Atmosphere.from_file(zlay, zlev, psurf, atmosphere_path)
    assert abs(get_gas(atmosphere, 'O3')).sum() == approx(1.1288725e+22)
    assert abs(get_gas(atmosphere, 'H2O')).sum() == approx(3.286886e+27)
    assert abs(get_gas(atmosphere, 'CO2')).sum() == approx(1.400481e+26)
    assert abs(get_gas(atmosphere, 'NO2')).sum() == approx(9.760928e+18)
    assert abs(get_gas(atmosphere, 'O2')).sum() == approx(8.883318e+28)
    assert abs(get_gas(atmosphere, 'CH4')).sum() == approx(7.632445e+23)


def test_get_AFGL_atm_homogenous_distribution(atmosphere_path):
    atmosphere = get_AFGL_atm_homogenous_distribution(
        atmosphere_path, 20, 1000)
    assert abs(get_gas(atmosphere, 'O3')).sum() == approx(5.8943662e+22)
    assert abs(get_gas(atmosphere, 'H2O')).sum() == approx(2.129261e+27)
    assert abs(get_gas(atmosphere, 'CO2')).sum() == approx(8.623509e+25)
    assert abs(get_gas(atmosphere, 'NO2')).sum() == approx(2.639252e+19)
    assert abs(get_gas(atmosphere, 'O2')).sum() == approx(4.460803e+28)
    assert abs(get_gas(atmosphere, 'CH4')).sum() == approx(3.832671e+23)
