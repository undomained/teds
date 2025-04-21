# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Types used by SGM."""
from dataclasses import dataclass
from netCDF4 import Dataset
from teds import log
from typing import Self
import numpy as np
import numpy.typing as npt
import sys


@dataclass
class Gas:
    """Concentration and other attributes of a greenhouse gas"""
    name: str
    source: npt.NDArray[np.floating] | None
    emission_in_kgps: float | None
    concentration: npt.NDArray[np.floating]  # (x, y, z)


@dataclass
class Meteo:
    """Class to hold meteorological data for one of more gases"""
    crs: str
    dx: float
    dy: float
    dz: float
    x: npt.NDArray[np.floating]
    y: npt.NDArray[np.floating]
    z: npt.NDArray[np.floating]
    lat: npt.NDArray[np.floating]
    lon: npt.NDArray[np.floating]
    znodes: npt.NDArray[np.floating]
    gases: list[Gas]

    @classmethod
    def from_empty(cls) -> Self:
        return cls('',
                   0,
                   0,
                   0,
                   np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   np.empty(0),
                   [])

    @classmethod
    def from_files(cls, filenames: list[str]) -> Self:
        """Initialize from one or more plume files."""
        gases = []
        # Grid needs to be read only once
        read_grid = True
        for filename in filenames:
            nc = Dataset(filename)
            gas_names = list(filter(
                lambda x: x not in
                ('time', 'z', 'x', 'y', 'longitude', 'latitude', 'znodes',
                 'U', 'V', 'W'),
                nc.variables.keys()))
            for gas_name in gas_names:
                if read_grid:
                    x, y, z = nc['x'][:].data, nc['y'][:].data, nc['z'][:].data
                    dz, dy, dx = z[1] - z[0], y[1] - y[0], x[1] - x[0]
                    znodes = nc['znodes'][:].data
                    lat = nc['latitude'][:].data
                    lon = nc['longitude'][:].data
                    read_grid = False
                # Input data is in the order (t, z, y, x)
                concentration = np.swapaxes(nc[gas_name][0].data, 0, 2)
                gases.append(Gas(gas_name,
                                 nc[gas_name].source,
                                 nc[gas_name].emission_in_kgps,
                                 concentration))
        return cls('', dx, dy, dz, x, y, z, lat, lon, znodes, gases)

    def get_gas(self, name: str) -> Gas:
        res = list(filter(lambda x: x.name == name, self.gases))
        if not res:
            log.error(f'object does not contain {name}')
            sys.exit(1)
        return res[0]
