# =============================================================================
# geometry module for different E2E simulator profiles
#
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================
import numpy as np
import sys
import netCDF4 as nc
import yaml


def gm_output(filename, vza, vaa, sza, saa, lat_grid, lon_grid,):

    nact = len(lat_grid[0])
    nalt = len(lat_grid)

    output = nc.Dataset(filename, mode='w')

    output.title = 'Tango Carbon E2ES GM output'

    output.createDimension('bins_across_track', nact)     # across track axis
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


def geometry_module(config):
    """
    Geometry module to specify geometry.

    """

    # sys.path.append(global_config['path']['e2es_path']+'lib')
    from end_to_end.lib import libGM
    from end_to_end.lib import constants

    # the gm output is orginazed in dictionaries of the format dic[nalt, nact]

    ninit = 0
    if config['profile'] == "individual_spectra":
        # Generate the gm output to calculate E2E performance for individual spectra
        # first check consistencies of 'indivual_spectra' nput.

        nn = len(config['sza'])

        ns = (
            nn
            + len(config['saa']) + len(config['vza']) \
            + len(config['vaa'])+ len(config['albedo'])) / 5

        if nn != ns:
            sys.exit("input error in gm, code 1")

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

        sza[0, :] = config['sza'][:]
        saa[0, :] = config['saa'][:]
        vza[0, :] = config['vza'][:]
        vaa[0, :] = config['vaa'][:]
        ninit = ninit + 1

    if config['profile'] == "single_swath":

        ncheck = len(config['sza']) + len(config['saa']) + \
                 len(config['vza']) + len(config['vaa'])
       
        if (ncheck != 4*config['numb_atm_scenes']):
            sys.exit("input error in gm, code 2")

        for iscen in range(config['numb_atm_scenes']+1):
            outofrange = (config['scene_trans_index'][iscen] > 99) & \
                (config['scene_trans_index'][iscen] < 0) 
            if(outofrange):
                sys.exit('config parameter scene_trans_index out of range')

            
        nact = config["field_of_regard"]["nact"]
        nalt = 1

        lon_grid = np.empty([nalt, nact])
        lat_grid = np.empty([nalt, nact])
        sza = np.empty([nalt, nact])
        saa = np.empty([nalt, nact])
        vza = np.empty([nalt, nact])
        vaa = np.empty([nalt, nact])

        lat_grid[0][:] = np.nan
        lon_grid[0][:] = np.nan

        for iscen in range(config['numb_atm_scenes']):
            ind_start = config['scene_trans_index'][iscen]
            ind_end   = config['scene_trans_index'][iscen+1]
            sza[0, ind_start:ind_end] = config['sza'][iscen]
            saa[0, ind_start:ind_end] = config['saa'][iscen]
            vza[0, ind_start:ind_end] = config['vza'][iscen]
            vaa[0, ind_start:ind_end] = config['vaa'][iscen]

        ninit = ninit + 1
    
    if (config['profile'] == "S2_microHH"):
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

        ninit = ninit + 1

    if ninit != 1:
        sys.exit("something went wrong in gm, code=3, ninit = " +  str(ninit))

    # write data to output file
    gm_output(config['output'], vza, vaa, sza, saa, lat_grid, lon_grid)

    print(
        "=>gm calcultion finished successfully. ")
    return


if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    geometry_module(config)
