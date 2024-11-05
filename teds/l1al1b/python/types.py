# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Types and other constants used in the L1B processor.

Much of this is also used by the instrument model (IM).

"""
from enum import auto
from enum import IntEnum
from typing import TypedDict

import numpy as np
import numpy.typing as npt


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


class BinningTable(TypedDict):
    """Binning table indices and counts"""
    # Index of each detector pixel on the binned array
    bin_indices: npt.NDArray[np.int32]
    # Multiplicity of each binned pixel
    count_table: npt.NDArray[np.int32]


class CKDDark(TypedDict):
    """Dark offset and dark current CKD."""
    # Dark offset (independent of integration time)
    offset: npt.NDArray[np.float64]
    # Dark current per second of integration time
    current: npt.NDArray[np.float64]


class CKDNoise(TypedDict):
    """Noise CKD."""
    # Conversion gain (signal dependent noise term)
    conversion_gain: npt.NDArray[np.float64]
    # Read noise (signal independent noise term)
    read_noise: npt.NDArray[np.float64]


class CKDNonlin(TypedDict):
    """Signal nonlinearity CKD."""
    # Expected, linear signal in counts
    expected: npt.NDArray[np.float64]
    # Observed, nonlinear signal in counts
    observed: npt.NDArray[np.float64]


class CKDPRNU(TypedDict):
    """Photoresponse non-uniformity (PRNU) CKD."""
    # PRNU including quantum efficiency
    prnu_qe: npt.NDArray[np.float64]


class CKDStray(TypedDict, total=False):
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


class CKDSwath(TypedDict):
    """CKD related to the satellite swath."""
    # Across track angles
    act_angles: npt.NDArray[np.float64]
    # Intermediate wavelengths after ISRF convolution
    wavelengths: npt.NDArray[np.float64]
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


class CKDSpectral(TypedDict):
    """Spectral CKD."""
    # Wavelengths assigned to each detector column of each L1B
    # spectrum
    wavelengths: npt.NDArray[np.float64]


class CKDRadiometric(TypedDict):
    """Radiometric CKD."""
    # Radiometric calibration (correction) constants
    rad_corr: npt.NDArray[np.float64]


class CKD(TypedDict, total=False):
    """Dictionary describing all calibration parameter for Tango Carbon."""
    # Number of detector pixels in the spatial direction
    n_detector_rows: int
    # Number of detector pixels in the spectral direction
    n_detector_cols: int
    # Bad pixel mask
    pixel_mask: npt.NDArray[bool]

    dark: CKDDark
    noise: CKDNoise
    nonlin: CKDNonlin
    prnu: CKDPRNU
    stray: CKDStray
    swath: CKDSwath
    spectral: CKDSpectral
    radiometric: CKDRadiometric


class Geometry(TypedDict):
    """Viewing and solar geometry describing a slice or full orbit."""
    latitude: npt.NDArray[np.float64]
    longitude: npt.NDArray[np.float64]
    height: npt.NDArray[np.float64]
    # Solar/viewing azimuth/zenith angles
    saa: npt.NDArray[np.float64]
    sza: npt.NDArray[np.float64]
    vaa: npt.NDArray[np.float64]
    vza: npt.NDArray[np.float64]


class L1(TypedDict, total=False):
    """Partially or fully calibrated level 1 data.

    The data level can range from L1A to L1B depending on the
    calibration level.

    """
    # Current calibration or processing level
    proc_level: ProcLevel

    # Signal and noise of each detector pixel
    image_i32: npt.NDArray[np.int32]
    image: npt.NDArray[np.float64]
    noise: npt.NDArray[np.float64]
    # Once spectra have been extracted from the detector, the signal
    # and noise arrays are discarded and we work with wavelengths and
    # spectra and their noise levels.
    wavelengths: npt.NDArray[np.float64]
    spectra: npt.NDArray[np.float64]
    spectra_noise: npt.NDArray[np.float64]
    solar_irradiance: npt.NDArray[np.float64]

    # Detector and other settings
    timestamps: npt.NDArray[np.float64]
    binning_table_ids: npt.NDArray[np.int32]
    coad_factors: npt.NDArray[np.int32]
    exptimes: npt.NDArray[np.float64]

    geometry: Geometry
