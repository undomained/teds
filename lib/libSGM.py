# library SGM functions

import sys
import numpy as np
import scipy
from copy import deepcopy
from shapely import Point
from shapely import Polygon
import geopandas
from geopandas import GeoDataFrame
from pystac_client import Client
from scipy.interpolate import RegularGridInterpolator
import rioxarray
import shapely
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from scipy.signal import fftconvolve

from end_to_end.lib import libRT
from end_to_end.lib import libNumTools
from end_to_end.lib import hapi

def get_raw_sentinel2_data(lat, lon, S2_reading_log):

    # generate a geometry object with the latitude-longitude points
    latlon = [Point(xy) for xy in zip(lon.flatten(), lat.flatten())]
    # link the data point to the WGS84 reference ellipsiod
    gdf = GeoDataFrame(geometry=latlon, crs="EPSG:4326")

#    fig=gdf.explore(width=600, height=600)

    # Include a small margin because later when we convolve the scene the areas
    # near the granule borders become zero and we want to cut those parts out.
    target_box = gdf.unary_union.convex_hull.buffer(0.01)

    # Specify a STAC API. Sentinel 2 and Landsat collections are stored at AWS.
    api = Client.open('https://earth-search.aws.element84.com/v0')
    if(S2_reading_log):
        print('This API provides access to the following STAC catalogs')
        for link in api.get_links('child'):
            print(link)

    # Search for granules that intersect fully or partially with the bounding box.
    # The variable collection is a collection of just the meta data. It would be
    # too much to download the full data of all matching granules. We'll filter the
    # list first and then download only one of the granules.

    search = api.search(
        # max_items=10,
        collections=('sentinel-s2-l2a-cogs',),
        query={
            'eo:cloud_cover': {'lt': 0.1},  # filter by cloud coverage
            'sentinel:valid_cloud_cover': {'eq': True}, },
        intersects=target_box,)
    collection = search.get_all_items()

    # Narrow the list of granules down to the minimal number that touch the bounding box.
    # Consider the intersect bof ounding data granule with the target area.
    # Sum up intersecttions of all granules to optimize coverage.

    # We start with an empty polygon for the combined intersects
    all_boxes = Polygon()
    # List of granules (metadata) being considered in the end
    collection_filtered = []
    for i_granule, granule in enumerate(collection):  # loop over all extracted granules
        # Crop the bounding box of this granule with the target bounding box
        intersect = Polygon(granule.geometry['coordinates'][0]).intersection(target_box)
        # If adding this box to the list of previously accepted box increases
        # the total area by a certain margin then accept this box.
        all_boxes_cur = shapely.ops.unary_union([all_boxes, intersect])
        if abs(all_boxes_cur.area - all_boxes.area) > 0.02:
            all_boxes = all_boxes_cur
            collection_filtered.append(granule)
        if(S2_reading_log):
            print(f'Number matched: {search.matched()}')
            print(f'Target area: {traget_box.area}')
            print(f'Current area: {all_boxes.area}')
            print(f'Number of granules: {len(collection_filtered)}')

    # Have a quick look at the accepted granule(s)
    #    overview = rioxarray.open_rasterio(collection_filtered[0].assets['overview'].get_absolute_href())
    #    overview = overview.rio.reproject('EPSG:4326', shape=(overview.shape[1], overview.shape[2]))
    #   overview.plot.imshow()

    # Extract the high resolution albedo map of a selected wavelength (B07 is 783 nm)
    S2_albedo = rioxarray.open_rasterio(collection_filtered[0].assets['B11'].get_absolute_href())
    S2_ssd = collection_filtered[0].assets['B11'].extra_fields['gsd']
    return(S2_albedo, S2_ssd)


def get_sentinel2_albedo(gm_data, conf):

    lon = gm_data['lon']
    lat = gm_data['lat']
    S2_reading_log = False
    S2_albedo_raw, S2_ssd = get_raw_sentinel2_data(lat, lon, S2_reading_log)

    # Note that the S2 data are scaled by a factor 1.E4
    # store data

    # Note that the albedo values need to be divided by 10,000.
    if(S2_reading_log):
        S2_albedo_raw.plot(robust=True)

    # Define the settings for the convolution
    conv_settings = {}
    if(conf['kernel_parameter']['type'] == '2D Gaussian'):
        fwhm_x = conf['kernel_parameter']['fwhm_x']
        fwhm_y = conf['kernel_parameter']['fwhm_y']
        fsize = conf['kernel_parameter']['size_factor']

        conv_settings['type'] = conf['kernel_parameter']['type']
        conv_settings['1D kernel extension'] = np.int0(fsize*np.max([fwhm_x, fwhm_y])/S2_ssd)
        # convert all kernel parameter in units of sampling distance
        conv_settings['fwhm x'] = np.int0(fwhm_x/S2_ssd)
        conv_settings['fwhm y'] = np.int0(fwhm_y/S2_ssd)

    # copy data type
    S2_albedo_conv = deepcopy(S2_albedo_raw)
    S2_albedo_conv.data[0, :, :] = libNumTools.convolution_2d(S2_albedo_raw.data[0, :, :], conv_settings)

    if(S2_reading_log):
        S2_albedo_conv.plot(robust=True)

#   Change coordinate system to WGS84
    S2_albedo_resampled = S2_albedo_conv.rio.reproject('EPSG:4326')
    if(S2_reading_log):
        S2_albedo_resampled[:, :].plot(robust=True)

    # Extract data on target grid
    # Define an interpolating function interp such that interp(lat,lon) is an
    # interpolated value. Note that data.y and data.x are lat, long coordinates.
    interp = RegularGridInterpolator((S2_albedo_resampled.y, S2_albedo_resampled.x),
                                     S2_albedo_resampled.values[0], method='cubic')

    # The 2D interpolator only works with 1D lat/lon grids.
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    target_points = np.zeros((len(lat_flat), 2))
    for i_point in range(len(target_points)):
        target_points[i_point] = [lat_flat[i_point], lon_flat[i_point]]

    albedo = np.reshape(interp(target_points), lat.shape)

    # Note that the albedo values need to be divided by 10,000.

    albedo = albedo / 1e4

    return(albedo)
