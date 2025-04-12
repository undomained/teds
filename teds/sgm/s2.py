# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Operations for fetching and manipulating Sentinel 2 granules."""
from geopandas import GeoDataFrame
from io import BytesIO
from itertools import groupby
from netCDF4 import Dataset
from pystac.item import Item as PystacItem
from pystac_client import Client
from requests import get
from rioxarray.merge import merge_arrays
from shapely import Point
from shapely import Polygon
from tqdm import tqdm
from typing import cast
from xarray import DataArray
import datetime
import numpy as np
import numpy.typing as npt
import os
import rioxarray
import shapely

from teds import log
from teds.gm.io import read_geometry


def fetch_granules(
        lat: npt.NDArray[np.floating],
        lon: npt.NDArray[np.floating],
        date_range: tuple[datetime.date, datetime.date],
        max_cloud_cover: float = 0.1) -> list[PystacItem]:
    """Return S2 granules corresponding to a target box.

    Parameters
    ----------
    lat
        Latitudes of the target box
    lon
        Longitudes of the target box
    date_range
        Date range to narrow the search results
    max_cloud_cover
        Maximum cloud cover for each granule [0-100]

    Returns
    -------
        List of pystac objects representing the granule metadata. The
        actual albedos are not downloaded.

    """
    # Generate geometry from latitude-longitude points
    latlon = [Point(xy) for xy in zip(lon.ravel(), lat.ravel())]
    # Link points to the WGS84 reference ellipsiod
    gdf = GeoDataFrame(geometry=latlon, crs='EPSG:4326')
    # Include a small margin because later when we convolve the scene
    # the areas near the granule borders become zero and we want to
    # cut those parts out.
    target_box = gdf.unary_union.convex_hull.buffer(0.1)
    # Specify a STAC API. Sentinel 2 and Landsat collections are stored at AWS.
    api = Client.open('https://earth-search.aws.element84.com/v0')
    # Search for granules that intersect fully or partially with the
    # bounding box. Collection will contain just the meta data. It
    # would be too much to download the full contents of all matching
    # granules. Filter the list first and then download only one or
    # more granules.
    datetime_range = (
        datetime.datetime.combine(date_range[0], datetime.time()),
        datetime.datetime.combine(date_range[1], datetime.time()))
    search = api.search(
        collections=['sentinel-s2-l2a-cogs'],
        datetime=datetime_range,
        query={
            'eo:cloud_cover': {'lt': max_cloud_cover},
            'sentinel:valid_cloud_cover': {'eq': True}, },
        intersects=target_box,)
    log.info(f'Number of S2 granules found: {search.matched()}')
    collection = search.get_all_items()
    # Group S2 granules by date
    collection_grouped = []
    for _, group in groupby(collection, key=lambda x:
                            cast(datetime.datetime, x.get_datetime()).date()):
        granules = list(group)
        # Find the area of each group that is within the target box
        geometries = [Polygon(
            cast(dict[str, dict], granule.geometry)['coordinates'][0])
                      for granule in granules]
        group_polygon = (
            shapely.ops.unary_union(geometries).intersection(target_box))
        collection_grouped.append({
            'granules': granules,
            'area': group_polygon.area,
        })
    collection_sorted = sorted(
        collection_grouped, key=lambda x: x['area'], reverse=True)
    print('|-------------------------------------------------------------|')
    print('| Platform    | UTM zone | Cloud cover | Datetime             |')
    for item in collection_sorted:
        print(
            '|-------------+----------+-------------+----------------------|\n'
            f'|                 overlap area: {item["area"]:.5f} deg^2     '
            '            |')
        for granule in item['granules']:
            print(f'| {granule.properties["platform"]} |'
                  f'       {granule.properties["sentinel:utm_zone"]} |'
                  f'        {granule.properties["eo:cloud_cover"]:.2f} |'
                  f' {granule.properties["datetime"]} |'
                  )
    print('|-------------------------------------------------------------|')
    # Return the collection of granules that have the largest overlap
    # with the target box.
    return collection_sorted[0]['granules']


def write_albedo(
        albedo_file: str, albedos: list[DataArray], metadata: dict) -> None:
    nc = Dataset(albedo_file, 'w')
    nc.title = 'Sentinel 2 albedos for different wavelength bands'
    for albedo in albedos:
        nc_grp = nc.createGroup(albedo.band_label)
        nc_grp.createDimension('x', len(albedo.x))
        nc_grp.createDimension('y', len(albedo.y))

        nc_var = nc_grp.createVariable('gsd', 'f8')
        nc_var.long_name = 'ground sampling distance'
        nc_var.units = 'm'
        nc_var[:] = albedo.gsd

        nc_var = nc_grp.createVariable('x', 'f8', ['x'], fill_value=-32767)
        nc_var.long_name = 'easting'
        nc_var.units = 'm'
        nc_var.valid_min = 0.0
        nc_var.valid_max = 900e3
        nc_var[:] = albedo.x.values

        nc_var = nc_grp.createVariable('y', 'f8', ['y'], fill_value=-32767)
        nc_var.long_name = 'northing'
        nc_var.units = 'm'
        nc_var.valid_min = 0.0
        nc_var.valid_max = 10e6
        nc_var[:] = albedo.y.values

        if albedo.band_label == 'SCL':
            nc_var = nc_grp.createVariable(
                'scl', 'u1', ['y', 'x'], fill_value=0)
            nc_var.long_name = 'scene classification layer'
            nc_var.flag_values = [np.uint8(i) for i in range(12)]
            nc_var.flag_meanings = ('No Data (Missing data)',
                                    'Saturated or defective pixel',
                                    'Topographic casted shadows',
                                    'Cloud shadows',
                                    'Vegetation',
                                    'Not-vegetated',
                                    'Water',
                                    'Unclassified',
                                    'Cloud medium probability',
                                    'Cloud high probability',
                                    'Thin cirrus',
                                    'Snow or ice')
            nc_var.valid_min = 1
            nc_var.valid_max = 11
        else:
            nc_var = nc_grp.createVariable(
                'albedo', 'u2', ['y', 'x'], fill_value=32767)
            nc_var.long_name = 'albedo'
            nc_var.valid_min = 0
            nc_var.valid_max = 10_000
            nc_var.set_auto_scale(False)
            nc_var.scale_factor = 1e-4
        nc_var[:] = albedo.values[0]
    for k, v in metadata.items():
        setattr(nc, k, v)
    nc.close()


def download_albedo(config: dict) -> list[DataArray]:
    """Download Sentinel 2 albedo using a target area from GM output.

    The geometry file is given by config['io_files']['input_gm'].

    Parameters
    ----------
    config
        Configuration dictionary

    Returns
    -------
        Sentinel 2 albedos, each corresponding to a different
        wavelength band.

    """
    geometry = read_geometry(config['io_files']['input_gm'])
    return download_albedo_for_coords(config, geometry.lat, geometry.lon)


def download_albedo_for_coords(
        config: dict,
        lat: npt.NDArray[np.floating],
        lon: npt.NDArray[np.floating]) -> list[DataArray]:
    """Download Sentinel 2 albedo covering a target area.

    Parameters
    ----------
    config
        Configuration dictionary
    lat
        Latitudes of the target area
    lon
        Longitudes of the target area

    Returns
    -------
        Sentinel 2 albedos, each corresponding to a different
        wavelength band.

    """
    mcc = 0.1
    if "max_cloud_cover" in config["sentinel2"]:
        mcc = config["sentinel2"]["max_cloud_cover"]
    granules = fetch_granules(
        lat, lon, config['sentinel2']['date_range'], mcc)
    # Extract the high resolution albedo map of selected wavelength bands
    albedos = []
    first_crs = None
    for band_label in config['sentinel2']['band_label'] + ['SCL']:
        log.info(f'Downloading Sentinel 2 albedo for band {band_label}')
        # Fetch one or more S2 albedo maps at a given wavelength.
        # These will be merged into one so there will always be one
        # albedo per wavelength.
        albedos_partial = []
        for granule in granules:
            tiff_url = granule.assets[band_label].href
            short_name = f'{tiff_url.split("/")[-2]}/{tiff_url.split("/")[-1]}'
            resp = get(tiff_url, stream=True)
            in_memory_object = BytesIO()
            with tqdm(desc=short_name,
                      total=int(resp.headers.get('content-length', 0)),
                      unit='B',
                      unit_scale=True,
                      unit_divisor=1000) as bar:
                for data in resp.iter_content(chunk_size=1024):
                    size = in_memory_object.write(data)
                    bar.update(size)
            # Option to read the tiff file from memory or first write
            # to hard drive and read from there.
            if 'in_memory' in config['sentinel2'] and (
                    config['sentinel2']['in_memory']):
                albedo = rioxarray.open_rasterio(in_memory_object)
            else:
                _filename = 'albedo.tif'
                with open(_filename, 'bw') as tif_file:
                    in_memory_object.seek(0)
                    tif_file.write(in_memory_object.read())
                albedo = rioxarray.open_rasterio(_filename)
                os.remove(_filename)
            assert isinstance(albedo, DataArray)
            albedos_partial.append(albedo)
        # Merge albedos from different granules into one
        first_crs = albedos_partial[0].rio.crs
        for i in range(1, len(albedos_partial)):
            albedos_partial[i] = albedos_partial[i].rio.reproject(first_crs)
        albedo_combined = merge_arrays(albedos_partial)
        albedo_combined = albedo_combined.assign_attrs({
            'band_label': band_label,
            'gsd': abs(albedos_partial[0].x.values[1]
                       - albedos_partial[0].x.values[0])
        })
        albedos.append(albedo_combined)
    metadata = {
        'crs': str(first_crs),
        'granule_epsg': [
            granule.properties['proj:epsg'] for granule in granules
        ],
        'granule_datetime': [
            str(granule.properties['datetime']) for granule in granules
        ],
        'granule_product_id': [
            granule.properties['sentinel:product_id'] for granule in granules
        ],
        'granule_cloud_cover': [
            granule.properties['eo:cloud_cover'] for granule in granules
        ]
    }
    write_albedo(config['sentinel2']['albedo_file'], albedos, metadata)
    return albedos


def read_albedo(filename: str) -> list[DataArray]:
    """Read a list of Sentinel 2 albedos from a NetCDF file."""
    nc = Dataset(filename)
    albedos = []
    for group in [x for x in nc.groups if x != 'SCL']:
        albedo = DataArray(nc[group]['albedo'][:],
                           dims=('y', 'x'),
                           coords={
                               'y': nc[group]['y'][:],
                               'x': nc[group]['x'][:]
                               })
        albedo.attrs['gsd'] = nc[group]['gsd'][:]
        albedo.attrs['band_label'] = group
        if hasattr(nc, 'crs'):
            albedo.rio.write_crs(nc.crs, inplace=True)
        else:
            albedo.rio.write_crs(nc[group].crs, inplace=True)
        # Replace nan with closest non-nan value
        mask = np.isnan(albedo.values)
        idx = np.where(~mask, np.arange(mask.shape[1]), 0)
        np.maximum.accumulate(idx, axis=1, out=idx)
        albedo.values[mask] = albedo.values[np.nonzero(mask)[0], idx[mask]]
        albedos.append(albedo)
    return albedos


def read_scl(filename: str) -> DataArray:
    """Read Sentinel 2 surface classification layer from NetCDF file."""
    nc = Dataset(filename)
    grp = nc['SCL']
    scl = DataArray(grp['scl'][:],
                    dims=('y', 'x'),
                    coords={'y': grp['y'][:], 'x': grp['x'][:]})
    scl.attrs['gsd'] = grp['gsd'][:]
    scl.rio.write_crs(grp.crs, inplace=True)
    return scl
