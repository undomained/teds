# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.


import numpy as np
from scipy.interpolate import RegularGridInterpolator
from xarray import DataArray
from typing import List
from teds import log


def interp_sentinel2_albedo(s2_albedos: list[DataArray],
                            lat,
                            lon) -> List[DataArray]:

    s2_albedos_regridded = []
    for s2_albedo in s2_albedos:
        log.info(f'Sentinel 2 band {s2_albedo.band_label}:')

        #   Change coordinate system to WGS84
        log.info('  Projecting to WSG84')
        s2_albedo = s2_albedo.rio.reproject('EPSG:4326')
        s2_albedo = s2_albedo.rename({'x': 'lon', 'y': 'lat'})

        # Extract data on target grid. Define an interpolating
        # function interp such that interp(lat,lon) is an interpolated
        # value.
        log.info('  Interpolating to MicroHH grid')
        interp = RegularGridInterpolator(
            (s2_albedo.lat, s2_albedo.lon), s2_albedo.values, method='linear')

        target_points = np.array(list(zip(lat.ravel(), lon.ravel())))
        res = interp(target_points).reshape(lat.shape)

        crs = s2_albedo.rio.crs
        s2_albedo = DataArray(res,
                              dims=('y', 'x'),
                              coords={
                                  'lat': (['y', 'x'], lat),
                                  'lon': (['y', 'x'], lon)
                              },
                              attrs={
                                  'gsd': s2_albedo.gsd,
                                  'band_label': s2_albedo.band_label,
                              })
        s2_albedo.rio.write_crs(crs, inplace=True)

        # Add additional metadata
        central_wavelengths = {
            'B01': 442.1, 'B02': 492.4, 'B03': 559.8, 'B04': 664.6,
            'B05': 704.1, 'B06': 740.5, 'B07': 782.8, 'B08': 832.8,
            'B8A': 864.7, 'B09': 945.1, 'B10': 1372.5, 'B11': 1613.7,
            'B12': 2202.4,
        }
        bandwidths = {
            'B01': 21.0, 'B02': 66.0, 'B03': 36.0, 'B04': 31.0, 'B05': 15.0,
            'B06': 15.0, 'B07': 20.0, 'B08': 106.0, 'B8A': 21.0, 'B09': 20.0,
            'B10': 31.0, 'B11': 91.0, 'B12': 175.0
        }
        s2_albedo.attrs['central_wavelength'] = (
            central_wavelengths[s2_albedo.band_label])
        s2_albedo.attrs['bandwidth'] = (
            bandwidths[s2_albedo.band_label])
        s2_albedos_regridded.append(s2_albedo)

    return s2_albedos_regridded