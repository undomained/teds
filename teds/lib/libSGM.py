# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.


import numpy as np
from scipy.interpolate import RegularGridInterpolator
import logging
from xarray import DataArray
from .libNumTools import convolution_2d
from typing import List
from shapely import Point, Polygon
from geopandas import GeoDataFrame
from pystac_client import Client
from pystac.item import Item as PystacItem
from netCDF4 import Dataset
import rioxarray
import shapely
from requests import get
from tqdm import tqdm
from io import BytesIO
from teds import log


def generate_geometry(lat: List[float], lon: List[float]) -> List[PystacItem]:
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
        collections=('sentinel-s2-l2a-cogs',),
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
    collection_filtered = []
    for i_granule, granule in enumerate(collection):
        # Crop the bounding box of this granule with the target bounding box
        intersect = Polygon(granule.geometry['coordinates'][0]).intersection(
            target_box)
        # If adding this box to the list of previously accepted box increases
        # the total area by a certain margin then accept this box.
        all_boxes_cur = shapely.ops.unary_union([all_boxes, intersect])
        if abs(all_boxes_cur.area - all_boxes.area) > 0.02:
            all_boxes = all_boxes_cur
            collection_filtered.append(granule)
    return collection_filtered


def write_albedo_to_netcdf(albedo_file: str, albedos: List[DataArray]) -> None:
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

        nc_var = nc_grp.createVariable(
            'albedo', 'u2', ['x', 'y'], fill_value=32767)
        nc_var.long_name = 'albedo'
        nc_var.valid_min = 0
        nc_var.valid_max = 10_000
        nc_var.set_auto_scale(False)
        nc_var.scale_factor = 1e-4
        nc_var[:] = albedo.values[0]
    nc.close()


def download_sentinel2_albedo(config) -> None:
    nc = Dataset(config['io']['gm'])
    lat = nc['lat'][:]
    lon = nc['lon'][:]
    collection = generate_geometry(lat, lon)
    # Extract the high resolution albedo map of selected wavelength bands
    albedos = []
    for band_label in config['sentinel2']['band_section']:
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
        albedo = rioxarray.open_rasterio(in_memory_object)
        albedo = albedo.assign_attrs({
            'band_label': band_label,
            'gsd': abs(albedo.x.values[1] - albedo.x.values[0])
        })
        albedo = albedo.clip(min=1e-5*10_000, max=1.0*10_000)
        albedos.append(albedo)
    write_albedo_to_netcdf(config['sentinel2']['albedo_file'], albedos)


def interp_sentinel2_albedo(s2_albedos: List[DataArray],
                            lat,
                            lon) -> List[DataArray]:

    s2_albedos_regridded = []
    for s2_albedo in s2_albedos:
        log.info(f'Sentinel 2 band {s2_albedo.band_label}:')
        # # Define the settings for the convolution
        # conv_settings = {}
        # if (kernel_para['type'] == '2D Gaussian'):
        #     fwhm_x = kernel_para['fwhm_x']
        #     fwhm_y = kernel_para['fwhm_y']
        #     fsize = kernel_para['size_factor']

        #     conv_settings['type'] = kernel_para['type']
        #     conv_settings['1D kernel extension'] = int(
        #         fsize * np.max([fwhm_x, fwhm_y]) / s2_albedo.gsd)
        #     # Convert all kernel parameter to units of sampling distance
        #     conv_settings['fwhm x'] = int(fwhm_x / s2_albedo.gsd)
        #     conv_settings['fwhm y'] = int(fwhm_y / s2_albedo.gsd)
        #     log.info('  Convolving with Gaussian')
        #     s2_albedo.values = convolution_2d(s2_albedo.values, conv_settings)

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
            'B08A': 864.7, 'B09': 945.1, 'B10': 1372.5, 'B11': 1613.7,
            'B12': 2202.4,
        }
        bandwidths = {
            'B01': 21.0, 'B02': 66.0, 'B03': 36.0, 'B04': 31.0, 'B05': 15.0,
            'B06': 15.0, 'B07': 20.0, 'B08': 106.0, 'B08A': 21.0, 'B09': 20.0,
            'B10': 31.0, 'B11': 91.0, 'B12': 175.0
        }
        s2_albedo.attrs['central_wavelength'] = (
            central_wavelengths[s2_albedo.band_label])
        s2_albedo.attrs['bandwidth'] = (
            bandwidths[s2_albedo.band_label])
        s2_albedos_regridded.append(s2_albedo)

    return s2_albedos_regridded