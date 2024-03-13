# =============================================================================
# geometry module for different E2E simulator profiles
#
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================
import numpy as np
import sys
import os
import netCDF4 as nc
import yaml
from teds.lib.libWrite import writevariablefromname
from teds.lib.libOrbSim import Sensor, Satellite
import teds.lib.data_netcdf.data_netcdf as dn
import teds.lib.lib_utils as Utils
import datetime

def check_input(logger, nact, check_list, place):
    """
        Check the input for simple profiles (single swath and individual spectra).
        Length should be equal to nact.
    """

    for view in check_list:

        if nact != len(config['scene_spec'][view]):
#            sys.exit(f"input error in gm for {view} nact ({nact}) not equal to {view} length ({len(config['scene_spec'][view])}), {place}")
            error_string = f"input error in gm for {view} nact ({nact}) not equal to {view} length ({len(config['scene_spec'][view])}), {place}.\nCan not continue."
            logger.error(error_string)
            sys.exit(255)

def get_individual_spectra(logger, config):
    """
        Generate the gm output for individual spectra
        first check consistencies of 'indivual_spectra' input.
        Here we use the 2-dimensional (nalt, nact) data structure in an artificial way
    
        return vza : viewing zentih angle; numpy array
               vaa : viewing azimuth angle; numpy array
               sza : solar zenith angle; numpy array
               saa : solar azimuth angle; numpy array
               lat_grid : pixel latitude grid; numpy array
               lon_grid: pixel_longitude grid, numpy array
    """

    nn = len(config['scene_spec']['sza'])

    check_input(logger, nn, ['sza','saa','vza','vaa','albedo'], "code 1")

    # here we use the 2-dimensional data structure in an artificial way
    nact = nn
    nalt = 1

    lon_grid = np.empty([nalt, nact])
    lat_grid = np.empty([nalt, nact])
    sza = np.empty([nalt, nact])
    saa = np.empty([nalt, nact])
    vza = np.empty([nalt, nact])
    vaa = np.empty([nalt, nact])

    lat_grid[0, :] = np.nan
    lon_grid[0, :] = np.nan

    nact = len(lat_grid[0])
    nalt = len(lat_grid)

    sza[0, :] = config['scene_spec']['sza'][:]
    saa[0, :] = config['scene_spec']['saa'][:]
    vza[0, :] = config['scene_spec']['vza'][:]
    vaa[0, :] = config['scene_spec']['vaa'][:]

    return vza, vaa, sza, saa, lat_grid, lon_grid

# TODO get_single_swath and get_indivudual_spectra very much the same, Can code be merged?
def get_single_swath(logger, config):
    """
        Generate the gm output for single swath

        return vza : viewing zentih angle; numpy array
               vaa : viewing azimuth angle; numpy array
               sza : solar zenith angle; numpy array
               saa : solar azimuth angle; numpy array
               lat_grid : pixel latitude grid; numpy array
               lon_grid: pixel_longitude grid, numpy array
    """

    nact = config['field_of_regard']['nact']       
    nalt = 1

    check_input(logger, nact, ['sza','saa','vza','vaa','albedo'], "code 2")

    lon_grid = np.empty([nalt, nact])
    lat_grid = np.empty([nalt, nact])
    sza = np.empty([nalt, nact])
    saa = np.empty([nalt, nact])
    vza = np.empty([nalt, nact])
    vaa = np.empty([nalt, nact])

    lat_grid[0][:] = np.nan
    lon_grid[0][:] = np.nan

    sza[0, :] = config['scene_spec']['sza'][:]
    saa[0, :] = config['scene_spec']['saa'][:]
    vza[0, :] = config['scene_spec']['vza'][:]
    vaa[0, :] = config['scene_spec']['vaa'][:]

    return vza, vaa, sza, saa, lat_grid, lon_grid

def get_orbit(logger, config):
    """
        Generate the gm output for an orbit
        configure satellite and propagate orbit

        return vza : viewing zentih angle; numpy array
               vaa : viewing azimuth angle; numpy array
               sza : solar zenith angle; numpy array
               saa : solar azimuth angle; numpy array
               lat_grid : pixel latitude grid; numpy array
               lon_grid: pixel_longitude grid, numpy array
    """

    sat = orbit_simulation(logger, config)

    # configure sensors and compute ground pixel information
    sensors = sensor_simulation(logger, config, sat)

    vza = sensors['gpxs'][0]['vza']
    vaa = sensors['gpxs'][0]['vaa']
    sza = sensors['gpxs'][0]['sza']
    saa = sensors['gpxs'][0]['saa']
    lat_grid = sensors['gpxs'][0]['lat']
    lon_grid = sensors['gpxs'][0]['lon']

    return vza, vaa, sza, saa, lat_grid, lon_grid

def get_S2_microHH(logger, config):
    """
        Generate the gm output for S2 microHH

        return vza : viewing zentih angle; numpy array
               vaa : viewing azimuth angle; numpy array
               sza : solar zenith angle; numpy array
               saa : solar azimuth angle; numpy array
               lat_grid : pixel latitude grid; numpy array
               lon_grid: pixel_longitude grid, numpy array
    """
    from lib import libGM
    from lib import constants
    nact = config["field_of_regard"]["nact"]
    nalt = config["field_of_regard"]["nalt"]

    vza = np.empty([nalt, nact])
    vaa = np.empty([nalt, nact])

    # across track angles, assume equi-distant sampling
    alpha_act_min = config["field_of_regard"]["alpha_act_min"]
    alpha_act_max = config["field_of_regard"]["alpha_act_max"]
    delta_alpha = (alpha_act_max - alpha_act_min) / (nact - 1)

    alpha_act = (delta_alpha * np.arange(nact) + alpha_act_min) / 180.0 * np.pi

    # extract time and location data from YAML settings
    when = [
        config["time"]["year"],
        config["time"]["month"],
        config["time"]["day"],
        config["time"]["hour"],
        config["time"]["minute"],
        config["time"]["timezone"],
    ]
    lat_ref = config["geometry"]["lat_initial"]
    lon_ref = config["geometry"]["lon_initial"]

    # geocentric radius as a function of latitude for WGS84, see https://en.wikipedia.org/wiki/Earth_radius
    coslat = np.cos(lat_ref / 180.0 * np.pi)
    sinlat = np.sin(lat_ref / 180.0 * np.pi)

    Rearth = np.sqrt(
        ((constants.a_axis_earth**2 * coslat) ** 2 + (constants.b_axis_earth**2 * sinlat) ** 2)
        / ((constants.a_axis_earth * coslat) ** 2 + (constants.b_axis_earth * sinlat) ** 2)
    )

    # spatial sampling in ACT direction, see atbd for defintion of phi
    bcoeff = (
        -2
        * (Rearth + config["satellite"]["sat_height"])
        * np.cos(alpha_act + config["satellite"]["alpha_roll"] / 180.0 * np.pi)
    )
    ccoeff = (Rearth + config["satellite"]["sat_height"]) ** 2 - Rearth**2
    Lalpha = -0.5 * (bcoeff + np.sqrt(bcoeff**2 - 4 * ccoeff))
    cosphi = (Rearth**2 + (Rearth + config["satellite"]["sat_height"]) ** 2 - Lalpha**2) / (
        2.0 * Rearth * (Rearth + config["satellite"]["sat_height"])
    )
    phi = np.arccos(cosphi)
    sact = phi * Rearth
    # convention: negative distances to the west of the subsatellite point, positive distances to the east
    idx = (alpha_act + config["satellite"]["alpha_roll"] / 180.0 * np.pi) < 0.0
    sact[idx] = -sact[idx]

    # calculate viewing zenith and viewing azimuth angle. Here we assume, that these angles are the same for each scanline
    vza_tmp = (alpha_act) / np.pi * 180.0 + config["satellite"]["alpha_roll"]
    idx = vza_tmp < 0.0
    vza_tmp[idx] = -vza_tmp[idx]
    for ialt in range(nalt):
        vza[ialt, :] = vza_tmp
    # the longitude coordinate line and a swath perpendicular to this. So all points to the west of
    # the subsatellite have an azimuth angle of -90 degree, for points all to Earst it is 90 degree.
    vaa_tmp = np.zeros(alpha_act.size) + 90.0  # true for all LOS with sact > 0
    vaa_tmp[idx] = 270.0
    for ialt in range(nalt):
        vaa[ialt, :] = vaa_tmp

    # Spatial sampling in ALT direction
    # angular velocity
    radicand = constants.grav * constants.mearth / (Rearth + config["satellite"]["sat_height"]) ** 3
    omega_sat = np.sqrt(radicand)
    # satellite ground speed, term cos(phi) need as ground speed varies with swath position
    v_ground = omega_sat * Rearth * cosphi

    # relative spatial sampling distance
    salt = v_ground * config["time"]["time_incremental"]
    # spatial grid
    spat_grid = [salt, sact]
    # transformation to a (lat,lon) grid
    sza, saa, lat_grid, lon_grid = libGM.trans_lat_lon(lat_ref, lon_ref, when, spat_grid, nalt, nact)

    return vza, vaa, sza, saa, lat_grid, lon_grid

def gm_output(logger, config, vza, vaa, sza, saa, lat_grid, lon_grid):
    """
       Write gm oputput to filename (set in config file) as nc file.
    """

    filename = config['gm_file']
    # Check if directory exists, otherwise create:
    out_dir = os.path.split(filename)[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    title = config['gm_title']

    nact = len(lat_grid[0,:])
    nalt = len(lat_grid[:,0])
    output = nc.Dataset(filename, mode='w')
#    output.title = 'Tango Carbon E2ES GM output'
    output.title = title
    output.createDimension('bins_across_track', nact)    # across track axis
    output.createDimension('bins_along_track', nalt)     # along track axis
    # dimensions
    dims = ('bins_along_track', 'bins_across_track')
    _ = writevariablefromname(output, "solarzenithangle", dims, sza)
    _ = writevariablefromname(output, "solarazimuthangle", dims, saa)
    _ = writevariablefromname(output, "viewingzenithangle", dims, vza)
    _ = writevariablefromname(output, "viewingazimuthangle", dims, vaa)
    _ = writevariablefromname(output, "latitude", dims, lat_grid)
    _ = writevariablefromname(output, "longitude", dims, lon_grid)
    output.close()

def gm_output_via_object(logger, config, vza, vaa, sza, saa, lat_grid, lon_grid):
    """
       Write gm oputput to filename (set in config file) as nc file.
       Using the data_netCDF class
    """

    filename = config['gm_file']
    # Check if directory exists, otherwise create:
    out_dir = os.path.split(filename)[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    title = config['gm_title']

    nact = len(lat_grid[0,:])
    nalt = len(lat_grid[:,0])

    gm_data = dn.DataNetCDF(filename, title=title)
    gm_data.add('E2E_configuration', value=str(config), kind='attribute')

    dims = ('bins_along_track', 'bins_across_track')
    gm_data.add(name=dims[0], value=nalt, kind='dimension')    # along track axis
    gm_data.add(name=dims[1], value=nact, kind='dimension')     # across track axis
    gm_data.add(name="solarzenithangle", dimensions=dims, value=sza)
    gm_data.add(name="solarazimuthangle", dimensions=dims, value=saa)
    gm_data.add(name="viewingzenithangle", dimensions=dims, value=vza)
    gm_data.add(name="viewingazimuthangle", dimensions=dims, value=vaa)
    gm_data.add(name="latitude", dimensions=dims, value=lat_grid)
    gm_data.add(name="longitude", dimensions=dims, value=lon_grid)
    gm_data.write()


def orbit_simulation(logger, config):
    """
        define the orbit
    """
    sat = Satellite(logger, config['orbit'])

    # propagate the orbit
    logger.info('propagate orbit...')
    dt_start = config['orbit']['epoch']
    dt_end = dt_start + datetime.timedelta(hours=config['orbit']['propagation_duration'])
    # compute the satellite position for every 10 seconds
    satpos = sat.compute_position(dt_start, dt_end, dt_interval=10.0)

    return {'sat': sat, 'sat_pos': satpos, 'dt_start': dt_start, 'dt_end': dt_end}


def sensor_simulation(logger, config, sat):
    """
        propogate sensors
    """
    gpxs = []
    sensor_config = []
    pitch = []
    roll = []
    yaw = []

    for key in config['sensors'].keys():

        logger.info('defining sensor ' + key)

        # sensor config wrt satellite reference frame
        sensor_half_swath_deg = np.rad2deg(
            np.arctan(0.5 * config['sensors'][key]['swath_width'] / config['orbit']['sat_height']))
        sensor_interval_deg = 2 * sensor_half_swath_deg / config['sensors'][key]['n_ground_pixels']

        # across track angle range and interval (homogenous sampling assumed)
        act_angle_range = [-sensor_half_swath_deg, sensor_half_swath_deg, sensor_interval_deg]
        thetas = np.arange(act_angle_range[0]+0.5*sensor_interval_deg, act_angle_range[1], act_angle_range[2])

        # along track angle range - commented out because unused
        # alt_angle_range = [-np.rad2deg(np.arctan(0.5 * config['sensors'][key]['alt_sampling']/config['orbit']['sat_height'])),
        #                   np.rad2deg(np.arctan(0.5 * config['sensors'][key]['alt_sampling']/config['orbit']['sat_height']))]

        # time range of observations and interval
        dt_range = [sat['dt_start'] + datetime.timedelta(minutes=config['sensors'][key]['start_time']),
                    sat['dt_start'] + datetime.timedelta(minutes=config['sensors'][key]['end_time']),
                    config['sensors'][key]['integration_time']]

        # make and propage the sensor
        sensor = Sensor(logger, act_angle_range[0], act_angle_range[1], act_angle_range[2])
        sensor.compute_groundpoints(sat['sat_pos'], pitch=config['sensors'][key]['pitch'],
                                    yaw=config['sensors'][key]['yaw'], roll=config['sensors'][key]['roll'])

        logger.info('compute the ground pixels (gpx)')
        gpx = sensor.get_groundpoints(dt_range[0], dt_range[1], dt_range[2], thetas)

        # collect the output data
        sensor_config.append(config['sensors'][key])
        sensor_config[-1]['name'] = key

        gpxs.append(gpx)

        pitch.append(config['sensors'][key]['pitch'])
        yaw.append(config['sensors'][key]['yaw'])
        roll.append(config['sensors'][key]['roll'])

    return {'gpxs': gpxs, 'pitch': np.array(pitch), 'roll': np.array(roll), 'yaw': np.array(yaw)}


def interpolate_pitch(logger, sensors, pitch, i_time):
    """
        What does this do exactly?
    """
    # assumes that ordering of sensors['pitch'] is either monotonically increasing or monotonically decreasing
    sign = np.sign(sensors['pitch'][1] - sensors['pitch'][0])
    # note that it will extrapolate by using the edge values
    s_f = np.interp(sign * pitch, sign * sensors['pitch'][:], np.arange(sensors['pitch'].shape[0], dtype=float))

    idx = [np.floor(s_f).astype(int), np.ceil(s_f).astype(int)]
    w = [1 - s_f + idx[0], s_f - idx[0]]

    gpx = {}

    # interpolate 3D array
    keys = ['p']
    for key in keys:
        gpx[key] = w[0] * sensors['gpxs'][idx[0]][key][i_time, :, :] + w[1] * sensors['gpxs'][idx[1]][key][i_time, :, :]

    # interpolate 2D array
    keys = ['lat', 'lon', 'height',  'vza', 'vaa', 'sza', 'saa']
    for key in keys:
        gpx[key] = w[0] * sensors['gpxs'][idx[0]][key][i_time, :] + w[1] * sensors['gpxs'][idx[1]][key][i_time, :]

    # time sample
    keys = ['sat_p', 'sat_lat', 'sat_lon', 'sat_height', 'seconds_from_epoch']
    for key in keys:
        gpx[key] = sensors['gpxs'][idx[0]][key][i_time]

    # copy
    keys = ['thetas', 'start_time', 'end_time', 'epoch']
    for key in keys:
        gpx[key] = sensors['gpxs'][idx[0]][key]

    return gpx

def gm_output_old(logger, config, vza, vaa, sza, saa, lat_grid, lon_grid,):
    """
       Write  gm oputput to filename (set in config file) as nc file.
       why is this labeled old?
    """

    filename = config['output']
    # Check if directory exists, otherwise create:
    out_dir = os.path.split(filename)[0]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    title = config['output_title']

    nact = len(lat_grid[0])
    nalt = len(lat_grid)

    output = nc.Dataset(filename, mode='w')

#    output.title = 'Tango Carbon E2ES GM output'
    output.title = title

    output.createDimension('bins_across_track', nact)    # across track axis
    output.createDimension('bins_along_track', nalt)     # along track axis

    gm_sza = output.createVariable('sza', np.float64, ('bins_along_track', 'bins_across_track',))
    gm_sza.units = 'degree'
    gm_sza.long_name = 'solar zenith angle'
    gm_sza.valid_min = 0.
    gm_sza.valid_max = 90.
    gm_sza.FillValue = -32767
    gm_sza[:] = sza[:]

    gm_saa = output.createVariable('saa', np.float64, ('bins_along_track', 'bins_across_track',))
    gm_saa.units = 'degree'
    gm_saa.long_name = 'solar azimuth angle'
    gm_saa.valid_min = 0.
    gm_saa.valid_max = +360.
    gm_saa.FillValue = -32767
    gm_saa[:] = saa[:]

    gm_vza = output.createVariable('vza', np.float64, ('bins_along_track', 'bins_across_track',))
    gm_vza.units = 'degree'
    gm_vza.long_name = 'viewing zenith angle'
    gm_vza.valid_min = 0.
    gm_vza.valid_max = 90.
    gm_vza.FillValue = -32767
    gm_vza[:] = vza[:]

    gm_vaa = output.createVariable('vaa', np.float64, ('bins_along_track', 'bins_across_track',))
    gm_vaa.units = 'degree'
    gm_vaa.long_name = 'viewing azimuth angle'
    gm_vaa.valid_min = 0.
    gm_vaa.valid_max = +360.
    gm_vaa.FillValue = -32767
    gm_vaa[:] = vaa[:]

    gm_lat = output.createVariable('lat', np.float64, ('bins_along_track', 'bins_across_track',))
    gm_lat.units = 'degree'
    gm_lat.long_name = 'latitude'
    gm_lat.valid_min = -90.
    gm_lat.valid_max = +90.
    gm_lat.FillValue = -32767
    gm_lat[:] = lat_grid[:]

    gm_lon = output.createVariable('lon', np.float64, ('bins_along_track', 'bins_across_track',))
    gm_lon.units = 'degree'
    gm_lon.long_name = 'longitude'
    gm_lon.valid_min = -180.
    gm_lon.valid_max = +180.
    gm_lon.FillValue = -32767
    gm_lon[:] = lon_grid[:]

    output.close()

    return


def geometry_module(logger, config):
    """
    Geometry module to specify geometry.

    """

    # TODO: dict? they are numpy arays. Toch?
    # the gm output is orginazed in dictionaries of the format dic[nalt, nact]

    if config['profile'] == "individual_spectra":
        vza, vaa, sza, saa, lat_grid, lon_grid = get_individual_spectra(logger, config)

    elif (config['profile'] == "single_swath"):

        vza, vaa, sza, saa, lat_grid, lon_grid = get_single_swath(logger, config)
    
    elif (config['profile'] == "orbit"):
        # configure satellite and propagate orbit
        vza, vaa, sza, saa, lat_grid, lon_grid = get_orbit(logger, config)

    elif (config['profile'] == "S2_microHH"):
        vza, vaa, sza, saa, lat_grid, lon_grid = get_S2_microHH(logger, config)

    else:
        error_string = f"something went wrong in gm, code=3, unrecognized profile choise: {config['profile']}\nCan not continue. Stopping now!"
        logger.error(error_string)
        sys.exit(255)

    # write data to output file
    gm_output(logger, config, vza, vaa, sza, saa, lat_grid, lon_grid)

    logger.info("=>gm calculation finished successfully. ")
    return


if __name__ == '__main__' :

    # Get logger for GM
    gm_logger = Utils.get_logger()
    # Get configuration info
    cfgFile = sys.argv[1]
    config = Utils.getConfig(gm_logger, cfgFile)
    # Get information (like git hash and config file name and version (if available) 
    # that will be added to the output file as attributes
    main_attribute_dict = Utils.get_main_attributes(config, config_attributes_name='GM_configuration')

    geometry_module(gm_logger, config)

    # add attributes to the output file
    Utils.add_attributes_to_output(gm_logger, config['gm_file'], main_attribute_dict)

