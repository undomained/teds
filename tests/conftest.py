# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Fixtures and settings used by more than one modules."""
from pathlib import Path
import pytest
import numpy as np

from teds.l1al1b.types import CKDDark
from teds.l1al1b.types import CKDNoise
from teds.l1al1b.types import CKDNonlin
from teds.l1al1b.types import CKDPRNU
from teds.l1al1b.types import L1

_cur_dir = Path(__file__).parent
_fix_dir = _cur_dir / 'common'


# Fixtures used by the instrument model and L1A-L1B processor

@pytest.fixture
def signal_i4(scope='session'):
    signal = np.loadtxt(_fix_dir / 'detector_image.txt', 'i4')
    return signal.reshape((1, -1))


@pytest.fixture
def signal(signal_i4, scope='session'):
    return signal_i4.astype('f8')


@pytest.fixture
def ckd_dark(scope='session'):
    ckd: CKDDark = {
        'offset': np.loadtxt(_fix_dir / 'ckd_dark_offset.txt', 'f8'),
        'current': np.loadtxt(_fix_dir / 'ckd_dark_current.txt', 'f8')
    }
    return ckd


@pytest.fixture
def ckd_noise(scope='session'):
    ckd: CKDNoise = {
        'conversion_gain':
        np.loadtxt(_fix_dir / 'ckd_conversion_gain.txt', 'f8'),
        'read_noise': np.loadtxt(_fix_dir / 'ckd_read_noise.txt', 'f8')
    }
    return ckd


@pytest.fixture
def ckd_nonlin(scope='session'):
    data = np.loadtxt(_fix_dir / 'ckd_nonlin.txt', 'f8')
    ckd: CKDNonlin = {'expected': data[:, 0], 'observed': data[:, 1]}
    return ckd


@pytest.fixture
def ckd_prnu(scope='session'):
    ckd: CKDPRNU = {'prnu_qe': np.loadtxt(_fix_dir / 'ckd_prnu.txt', 'f8')}
    return ckd


@pytest.fixture
def l1(signal, scope='session'):
    l1: L1 = {
        'signal': signal,
        'noise': np.ones(signal.shape),
        'binning_table_ids': np.array([0]),
        'coad_factors': np.ones(signal.shape),
        'exptimes': np.array([0.0460833]),
    }
    return l1
