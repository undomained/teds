# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
# =============================================================================
#     geophysical scene generation module for different E2E simulator profiles
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================
import numpy as np
import scipy as sc  # .interpolate as interp
import time

import logging
_logger = logging.getLogger(__name__)


def interpolate_single_parameter_to_orbit(lat_orbit,
                                          lon_orbit,
                                          n_pixels,
                                          lat_in,
                                          lon_in,
                                          parameters_in,
                                          lonlat_gridded=False,
                                          method='nearest_linear'):
    parameter_orbit = np.zeros(n_pixels)
    parameter_orbit_tmp = np.zeros(n_pixels)

    if lonlat_gridded:
        x_sparse = lon_in[:, :]
        y_sparse = lat_in[:, :]
    else:
        x_sparse, y_sparse = np.meshgrid(lon_in, lat_in)

    # mark invalid values
    parameters_in = np.ma.masked_invalid(parameters_in)

    # get only the valid values
    x1 = x_sparse[~parameters_in.mask]
    y1 = y_sparse[~parameters_in.mask]

    t_start = time.perf_counter()

    if (len(x1) > 0):
        newarr = parameters_in[~parameters_in.mask]

        parameter_orbit_tmp[:] = sc.interpolate.griddata(
            (x1, y1),
            newarr,
            (lon_orbit.flatten(), lat_orbit.flatten()),
            method='nearest')
        t_inter = time.perf_counter()
        _logger.debug("Nearest interpolation done. Linear following.")
        parameter_orbit[:] = sc.interpolate.griddata(
            (x1, y1),
            newarr,
            (lon_orbit.flatten(), lat_orbit.flatten()),
            method='linear')

        if method == 'nearest':
            parameter_orbit[:] = sc.interpolate.griddata(
                (x1, y1),
                newarr,
                (lon_orbit.flatten(), lat_orbit.flatten()),
                method='nearest')
        if method == "nearest_linear":

            parameter_orbit_tmp[:] = sc.interpolate.griddata(
                (x1, y1),
                newarr,
                (lon_orbit.flatten(), lat_orbit.flatten()),
                method='nearest')

            t_inter = time.perf_counter()

            parameter_orbit[:] = sc.interpolate.griddata(
                (x1, y1),
                newarr,
                (lon_orbit.flatten(), lat_orbit.flatten()),
                method='linear')

            idx_nan = np.where((np.isnan(parameter_orbit)))
            parameter_orbit[idx_nan] = parameter_orbit_tmp[idx_nan]

        if method == 'linear':
            parameter_orbit[:] = sc.interpolate.griddata(
                (x1, y1),
                newarr,
                (lon_orbit.flatten(), lat_orbit.flatten()),
                method='linear')

        _logger.debug("interpolation succes.")

        tmp_arr = parameter_orbit[:]
        threshold_min = np.nanmin(newarr)
        threshold_max = np.nanmax(newarr)
        tmp_arr[tmp_arr < threshold_min] = threshold_min
        tmp_arr[tmp_arr > threshold_max] = threshold_max
        parameter_orbit[:] = tmp_arr[:]

    t_end = time.perf_counter()
    if method == 'nearest_linear':
        _logger.debug(
            "Time for interpolations: nearest {}, linear {}, total {}".format(
                t_inter - t_start, t_end - t_inter, t_end - t_start))
    else:
        _logger.debug("Time for interpolations: nearest {}".format(
            t_end - t_start))

    return parameter_orbit


def interpolate_to_orbit(lat_orbit,
                         lon_orbit,
                         n_pixels,
                         lat_in,
                         lon_in,
                         parameters_in,
                         n_parameters,
                         transpose=False,
                         lonlat_gridded=False,
                         keep_minmax=True,
                         method='nearest_linear'):

    _logger.debug('lat orbit {}'.format(lat_orbit.shape))
    _logger.debug('lon orbit {}'.format(lon_orbit.shape))
    _logger.debug('n_parameters {}'.format(n_parameters))
    _logger.debug('lat_in {}'.format(lat_in.shape))
    _logger.debug('lon_in {}'.format(lon_in.shape))

    orbit_data = np.zeros((n_parameters, n_pixels))
    orbit_data_tmp = np.zeros((n_parameters, n_pixels))
    if lonlat_gridded:
        x_sparse = lon_in[:, :]
        y_sparse = lat_in[:, :]
    else:
        x_sparse, y_sparse = np.meshgrid(lon_in, lat_in)

    for ipar in range(0, n_parameters):
        t_start = time.perf_counter()

        _logger.debug("collocation for i_parameter: {}".format(ipar))

        if transpose:
            parameters_in_i = np.transpose(parameters_in[ipar, :, :])
        else:
            parameters_in_i = parameters_in[ipar, :, :]

        # mark invalid values
        parameters_in_i = np.ma.masked_invalid(parameters_in_i)

        # get only the valid values
        x1 = x_sparse[~parameters_in_i.mask]
        y1 = y_sparse[~parameters_in_i.mask]

        if len(x1) > 0:
            newarr = parameters_in_i[~parameters_in_i.mask]

            if method == 'nearest':
                orbit_data[ipar, :] = sc.interpolate.griddata(
                    (x1, y1),
                    newarr,
                    (lon_orbit.flatten(), lat_orbit.flatten()),
                    method='nearest')
            if method == "nearest_linear":

                orbit_data_tmp[ipar, :] = sc.interpolate.griddata(
                    (x1, y1),
                    newarr,
                    (lon_orbit.flatten(), lat_orbit.flatten()),
                    method='nearest')

                t_inter = time.perf_counter()

                orbit_data[ipar, :] = sc.interpolate.griddata(
                    (x1, y1),
                    newarr,
                    (lon_orbit.flatten(), lat_orbit.flatten()),
                    method='linear')

                idx_nan = np.where((np.isnan(orbit_data[ipar, :])))
                orbit_data[ipar, idx_nan] = orbit_data_tmp[ipar, idx_nan]

            if method == 'linear':
                orbit_data[ipar, :] = sc.interpolate.griddata(
                    (x1, y1),
                    newarr,
                    (lon_orbit.flatten(), lat_orbit.flatten()),
                    method='linear')

            _logger.debug("interpolation succes.")

            if keep_minmax:

                tmp_arr = orbit_data[ipar, :]
                threshold_min = np.nanmin(newarr)
                threshold_max = np.nanmax(newarr)
                tmp_arr[tmp_arr < threshold_min] = threshold_min
                tmp_arr[tmp_arr > threshold_max] = threshold_max
                orbit_data[ipar, :] = tmp_arr[:]

            t_end = time.perf_counter()
            if method == 'nearest_linear':
                _logger.debug(
                    "Time for interpolations: nearest {}, linear {}, "
                    "total {}".format(t_inter - t_start,
                                      t_end - t_inter,
                                      t_end - t_start))
            else:
                _logger.debug("Time for interpolations: {} {}".format(
                    method, t_end - t_start))

    return orbit_data


def collocate_cloudmask(lat_orbit,
                        lon_orbit,
                        n_pixels,
                        lat_in,
                        lon_in,
                        cloudmask_in,
                        fraction_minimum=0.5,
                        lonlat_gridded=True):

    if lonlat_gridded:
        x_sparse = lon_in[:, :]
        y_sparse = lat_in[:, :]
    else:
        x_sparse, y_sparse = np.meshgrid(lon_in, lat_in)

    # Interpolation to orbit data for Cloud parameters. Outside range,
    # set to one, cloudy.
    cloud_mask_orbit = np.ones((n_pixels))
    cloud_fraction_orbit = np.ones((n_pixels))

    idx_cloud = np.where((cloud_mask_orbit >= fraction_minimum))
    idx_clear = np.where((cloud_mask_orbit < fraction_minimum))

    # mark invalid values
    cloud_mask_day_now = np.ma.masked_invalid(cloudmask_in)

    # get only the valid values
    x1 = x_sparse[~cloud_mask_day_now.mask]
    y1 = y_sparse[~cloud_mask_day_now.mask]
    t_start = time.perf_counter()

    if (len(x1) > 0):
        newarr = cloud_mask_day_now[~cloud_mask_day_now.mask]

        cloud_mask_orbit[:] = sc.interpolate.griddata(
            (x1, y1), newarr, (lon_orbit, lat_orbit), method='nearest')
        indices = np.where(
            (lat_orbit < np.min(lat_in)) | (lat_orbit > np.max(lat_in)) |
            (lon_orbit < np.min(lon_in)) | (lon_orbit > np.max(lon_in)))
        cloud_mask_orbit[indices] = 1
        cloud_fraction_orbit[:] = cloud_mask_orbit[:]
        # Note opposite def as in original cloud...idx_cloud !=1,
        # dx_no_cloud == 1?!)
        idx_cloud = np.where((cloud_mask_orbit >= fraction_minimum))
        idx_clear = np.where((cloud_mask_orbit < fraction_minimum))

        cloud_mask_orbit[idx_clear] = 0
        cloud_mask_orbit[idx_cloud] = 1
        _logger.debug('cloud percentage in orbit = {}'.format(
            np.sum(cloud_mask_orbit[idx_cloud]) / float(n_pixels)))

    t_end = time.perf_counter()
    _logger.debug("Time for interpolations: nearest total {}".format(
        t_end - t_start))

    ###############################################
    cloudmask_orbit = np.zeros((n_pixels)).astype(int)

    cloudmask_orbit[:] = cloud_mask_orbit[:].astype(int)

    return cloudmask_orbit, idx_cloud, idx_clear, cloud_fraction_orbit
