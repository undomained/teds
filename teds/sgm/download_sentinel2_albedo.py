# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

from geopandas import GeoDataFrame
from io import BytesIO
from shapely import Point
from shapely import Polygon
from pystac_client import Client
from pystac.item import Item as PystacItem
from netCDF4 import Dataset
from requests import get
from tqdm import tqdm
from xarray import DataArray
import numpy as np
import numpy.typing as npt
import os
import rioxarray
import shapely
import typing as tp

from teds import log


def generate_geometry(lat: DataArray,
                      lon: npt.NDArray[np.float64]) -> tp.List[PystacItem]:
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
    search = api.search(
        collections=['sentinel-s2-l2a-cogs'],
        query={
            'eo:cloud_cover': {'lt': 0.1},
            'sentinel:valid_cloud_cover': {'eq': True}, },
        intersects=target_box,)
    collection = search.get_all_items()
    # Narrow the list of granules to the minimal number that touch the
    # bounding box. Consider the intersection of each granule with the
    # target area and then sum up the intersections. Start with an
    # empty polygon for the combined intersects.
    all_boxes = Polygon()
    # List of granules (metadata) being considered in the end
    collection_filtered: list[PystacItem] = []
    for i_granule, granule in enumerate(collection):
        # Crop the bounding box of this granule with the target bounding box
        assert granule.geometry is not None
        coords = granule.geometry['coordinates']
        intersect = Polygon(coords[0]).intersection(
            target_box)
        # If adding this box to the list of previously accepted box increases
        # the total area by a certain margin then accept this box.
        all_boxes_cur = shapely.ops.unary_union([all_boxes, intersect])
        if abs(all_boxes_cur.area - all_boxes.area) > 0.02:
            all_boxes = all_boxes_cur
            collection_filtered.append(granule)
    return collection_filtered


def write_albedo_to_netcdf(albedo_file: str,
                           albedos: tp.List[DataArray]) -> None:
    nc = Dataset(albedo_file, 'w')
    nc.title = 'Sentinel 2 albedos for different wavelength bands'
    for albedo in albedos:
        nc_grp = nc.createGroup(albedo.band_label)
        nc_grp.createDimension('x', len(albedo.x))
        nc_grp.createDimension('y', len(albedo.y))

        nc_grp.crs = int(albedo.rio.crs.to_epsg())

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
                'scl', 'u1', ['x', 'y'], fill_value=0)
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
                'albedo', 'u2', ['x', 'y'], fill_value=32767)
            nc_var.long_name = 'albedo'
            nc_var.valid_min = 0
            nc_var.valid_max = 10_000
            nc_var.set_auto_scale(False)
            nc_var.scale_factor = 1e-4
        nc_var[:] = albedo.values[0]
    nc.close()


def download_sentinel2_albedo(config: dict) -> None:
    """Download Sentinel 2 albedo.

    Args:
      config: configuration settings

    """
    nc = Dataset(config['gm_input'])
    lat = nc['lat'][:]
    lon = nc['lon'][:]
    collection = generate_geometry(lat, lon)
    # Extract the high resolution albedo map of selected wavelength bands
    albedos = []
    for band_label in config['sentinel2']['band_section'] + ['SCL']:
        log.info(f'Downloading Sentinel 2 albedo for band {band_label}')
        tiff_url = collection[-1].assets[band_label].href
        short_name = tiff_url.split('sentinel-s2-l2a-cogs/')[1]
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
        # Option to read the tiff file from memory or first write to
        # hard drive and read from there.
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
        albedo = albedo.assign_attrs({
            'band_label': band_label,
            'gsd': abs(albedo.x.values[1] - albedo.x.values[0])
        })
        albedos.append(albedo)
    write_albedo_to_netcdf(config['sentinel2']['albedo_file'], albedos)
