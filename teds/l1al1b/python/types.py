# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Types used by various TEDS modules."""
from dataclasses import dataclass
from enum import IntEnum
from enum import auto
from typing import Self
import numpy as np
import numpy.typing as npt

from teds.gm.types import Geometry
from teds.gm.types import Navigation


# Process ladder, an ordered list of possible states of data. Data
# levels run from L1A to L1B.
class ProcLevel(IntEnum):
    l1a = 0
    raw = auto()
    dark_offset = auto()
    noise = auto()
    dark_current = auto()
    nonlin = auto()
    prnu = auto()
    stray = auto()
    swath = auto()
    l1b = auto()
    sgm = auto()

    def __str__(self) -> str:
        _names = (
            'L1A',
            'raw',
            'dark offset',
            'noise',
            'dark current',
            'nonlinearity',
            'PRNU',
            'stray light',
            'swath',
            'L1B',
            'SGM')
        return _names[self.value]


@dataclass
class BinningTable:
    """Binning table indices and counts"""
    # Index of each detector pixel on the binned array
    bin_indices: npt.NDArray[np.int32]
    # Multiplicity of each binned pixel
    count_table: npt.NDArray[np.int32]


@dataclass
class CKDDark:
    """Dark offset and dark current CKD."""
    # Dark offset (independent of integration time)
    offset: npt.NDArray[np.float64]
    # Dark current per second of integration time
    current: npt.NDArray[np.float64]


@dataclass
class CKDNoise:
    """Noise CKD."""
    # Conversion gain (signal dependent noise term)
    conversion_gain: npt.NDArray[np.float64]
    # Read noise (signal independent noise term)
    read_noise: npt.NDArray[np.float64]


@dataclass
class CKDNonlin:
    """Signal nonlinearity CKD."""
    # Observed, nonlinear signal in counts
    observed: npt.NDArray[np.float64]
    # Expected, linear signal in counts
    expected: npt.NDArray[np.float64]


@dataclass
class CKDPRNU:
    """Photoresponse non-uniformity (PRNU) CKD."""
    # PRNU including quantum efficiency
    prnu_qe: npt.NDArray[np.float64]


@dataclass
class CKDStray:
    """Stray light CKD."""
    # Fourier transforms of stray light kernels
    kernels_fft: list[npt.NDArray[np.complex128]]
    # Total internal scattering factor
    eta: npt.NDArray[np.float64]
    # Kernel weights. For a given kernel, its weight is 0 outside its
    # domain of influence.
    weights: npt.NDArray[np.float64]
    # Boundaries of subimages that must be extracted for the
    # convolutions. The order of coefficients is 'bottom', 'top',
    # 'left', 'right'.
    edges: npt.NDArray[np.int32]


@dataclass
class CKDSwath:
    """CKD related to the satellite swath."""
    # Across track angles
    act_angles: npt.NDArray[np.float64]
    # ACT angle of each detector pixel
    act_map: npt.NDArray[np.float64]
    # Wavelength of each detector pixel
    wavelength_map: npt.NDArray[np.float64]
    # Row index of each L1B spectral element
    row_map: npt.NDArray[np.float64]
    # Column index of each L1B spectral element
    col_map: npt.NDArray[np.float64]
    # Line of sight vectors
    line_of_sights: npt.NDArray[np.float64]


@dataclass
class CKDSpectral:
    """Spectral CKD."""
    # L1B wavelength grid, with uniform spacing. The same grid is
    # assigned to all L1B spectra.
    wavelengths: npt.NDArray[np.float64]


@dataclass
class CKDRadiometric:
    """Radiometric CKD."""
    # Radiometric calibration (correction) constants
    rad_corr: npt.NDArray[np.float64]


@dataclass
class CKD:
    """Dictionary describing all calibration parameter for Tango Carbon."""
    # Number of detector pixels in the spatial direction
    n_detector_rows: int
    # Number of detector pixels in the spectral direction
    n_detector_cols: int
    # Bad pixel mask
    pixel_mask: npt.NDArray[np.bool_]

    dark: CKDDark
    noise: CKDNoise
    nonlin: CKDNonlin
    prnu: CKDPRNU
    stray: CKDStray
    swath: CKDSwath
    spectral: CKDSpectral
    radiometric: CKDRadiometric


@dataclass
class L1:
    """Partially or fully calibrated level 1 data.

    The data level can range from L1A to L1B depending on the
    calibration level.

    """
    # Current calibration or processing level
    proc_level: ProcLevel

    # Signal and noise of each detector pixel
    signal: npt.NDArray[np.float64]
    noise: npt.NDArray[np.float64]
    # Once spectra have been extracted from the detector, the signal
    # and noise arrays are discarded and we work with wavelengths and
    # spectra and their noise levels.
    wavelengths: npt.NDArray[np.float64]
    spectra: npt.NDArray[np.float64]
    spectra_noise: npt.NDArray[np.float64]
    solar_irradiance: npt.NDArray[np.float64]

    # Detector image attributes
    time_units: str
    tai_seconds: npt.NDArray[np.uint]
    tai_subsec: npt.NDArray[np.float64]
    binning_table_id: int
    coad_factor: int
    exposure_time: float

    # Navigation data interpolated to detector timestamps
    navigation: Navigation

    # Geolocation result from navigation data
    geometry: Geometry

    @classmethod
    def from_empty(cls) -> Self:
        return cls(proc_level=ProcLevel.l1a,
                   signal=np.empty(0),
                   noise=np.empty(0),
                   wavelengths=np.empty(0),
                   spectra=np.empty(0),
                   spectra_noise=np.empty(0),
                   solar_irradiance=np.empty(0),
                   time_units='seconds since 2025-01-01',
                   tai_seconds=np.empty(0, dtype=np.uint),
                   tai_subsec=np.empty(0),
                   binning_table_id=0,
                   coad_factor=1,
                   exposure_time=0.0,
                   navigation=Navigation.from_shape((0,)),
                   geometry=Geometry.from_shape((0,)))
