#====================================================================
# simple routine to get convolved S2 data for GM scene area
#====================================================================
import sys
from copy import deepcopy
import netCDF4 as nc
import numpy as np
from teds.lib import libSGM, libNumTools
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt


def get_gm_data(filename):
    input = nc.Dataset(filename, mode='r')
    gm_data = {}
    gm_data['sza'] = deepcopy(input['sza'][:, :])
    gm_data['saa'] = deepcopy(input['saa'][:, :])
    gm_data['vza'] = deepcopy(input['vza'][:, :])
    gm_data['vaa'] = deepcopy(input['vaa'][:, :])
    gm_data['lat'] = deepcopy(input['lat'][:, :])
    gm_data['lon'] = deepcopy(input['lon'][:, :])
    input.close()
    return gm_data


def get_sentinel2_albedo(conf):

    gm_data = get_gm_data(conf['gm_input'])

    lon = gm_data['lon']
    lat = gm_data['lat']
    S2_reading_log = False
    print('Getting raw S2 data...')
    S2_albedo_raw, S2_ssd = libSGM.get_raw_sentinel2_data(lat, lon, S2_reading_log)
    print('                   ...done')
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

#    plt.contourf(lon.T, lat.T, albedo.T)
#    plt.show()

    return(albedo)

