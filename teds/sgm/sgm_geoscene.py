# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
# =============================================================================
#     geophysical scene generation module for different E2E simulator profiles
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================

import sys
import netCDF4 as nc
import numpy as np
import yaml
import logging
from ..lib import constants
from ..lib import libATM, libSGM
from ..lib.libWrite import writevariablefromname

from netCDF4 import Dataset
from xarray import DataArray
from typing import List

from ..lib.remotap_preproc.exceptions import ProcessError, chain_exception_message
from ..lib.remotap_preproc.aerosol_echam import AerosolEcham
from ..lib.remotap_preproc.ocean import Ocean
from ..lib.remotap_preproc.refractive_index import RefractiveIndex
from ..lib.remotap_preproc.xbpdf_polder import XbpdfPolder
from ..lib.remotap_preproc.xbrdf_gome2_modis import XbrdfGomeModis
from ..lib.remotap_preproc.cloudmask_s2 import CloudmaskS2
from ..lib.remotap_preproc.albedo_s2 import AlbedoS2
from ..lib.remotap_preproc.tool_algorithm import convert2julian7

# Copied verbatim from CO2M preprocessing script by Maud van den Broek,
# https://bitbucket.org/sron_earth/co2m_sgm/src/master/sources/co2mpreproc/processing_modules/aerosol_echam.py
# DO NOT REARRANGE THIS LIST; THE ORDER MATTERS!
varname_aerosol_out = [
    'sphere_mode4', 'sphere_mode5', 'sphere_mode6', 'sphere_mode7', 'sphere_mode1', 'sphere_mode2',
    'sphere_mode3',
    'num_mode4', 'num_mode5', 'num_mode6', 'num_mode7', 'num_mode1', 'num_mode2', 'num_mode3',
    'reff_mode4', 'reff_mode5', 'reff_mode6', 'reff_mode7', 'reff_mode1', 'reff_mode2', 'reff_mode3',
    'vfrac_species1_mode4', 'vfrac_species1_mode5', 'vfrac_species1_mode6', 'vfrac_species1_mode7',
    'vfrac_species1_mode1', 'vfrac_species1_mode2', 'vfrac_species1_mode3',
    'vfrac_species2_mode4', 'vfrac_species2_mode5', 'vfrac_species2_mode6', 'vfrac_species2_mode7',
    'vfrac_species2_mode1', 'vfrac_species2_mode2', 'vfrac_species2_mode3',
    'vfrac_species3_mode4', 'vfrac_species3_mode5', 'vfrac_species3_mode6', 'vfrac_species3_mode7',
    'vfrac_species3_mode1', 'vfrac_species3_mode2', 'vfrac_species3_mode3',
    'vfrac_species4_mode4', 'vfrac_species4_mode5', 'vfrac_species4_mode6', 'vfrac_species4_mode7',
    'vfrac_species4_mode1', 'vfrac_species4_mode2', 'vfrac_species4_mode3',
    'vfrac_species5_mode4', 'vfrac_species5_mode5', 'vfrac_species5_mode6', 'vfrac_species5_mode7',
    'vfrac_species5_mode1', 'vfrac_species5_mode2', 'vfrac_species5_mode3',
    'vfrac_species6_mode4', 'vfrac_species6_mode5', 'vfrac_species6_mode6', 'vfrac_species6_mode7',
    'vfrac_species6_mode1', 'vfrac_species6_mode2', 'vfrac_species6_mode3',
    'ALH_mode4', 'ALH_mode5', 'ALH_mode6', 'ALH_mode7', 'ALH_mode1', 'ALH_mode2', 'ALH_mode3',
    'veff_mode4', 'veff_mode5', 'veff_mode6', 'veff_mode7', 'veff_mode1', 'veff_mode2', 'veff_mode3',
    'pressure_surface_aerosol', 'elevation_surface_aerosol']


def get_sentinel2_albedo(filename: str) -> List[DataArray]:
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
        albedo.rio.write_crs(nc[group].crs, inplace=True)
        albedos.append(albedo)
    return albedos


def get_sentinel2_scl(filename: str) -> DataArray:
    """Read Sentinel 2 surface classification layer from NetCDF file."""
    nc = Dataset(filename)
    grp = nc['SCL']
    scl = DataArray(grp['scl'][:],
                    dims=('y', 'x'),
                    coords={'y': grp['y'][:], 'x': grp['x'][:]})
    scl.attrs['gsd'] = grp['gsd'][:]
    scl.rio.write_crs(grp.crs, inplace=True)
    return scl


class Emptyclass:
    """Empty class. Data container."""

    pass


def get_gm_data(filename):

    names = ['sza', 'saa', 'vza', 'vaa', 'lat', 'lon']

    input = nc.Dataset(filename, mode='r')

    gm_data = Emptyclass()

    for name in names:
        gm_data.__setattr__(name, input[name][:])

    input.close()

    return gm_data


def geosgm_output(
    filename,
    atm,
    atm_std,
    refractive_index = None,
    aerosol_echam = None,
    xbrdf = None,
    xbpdf = None
):
    # write geophysical scene data to output

    nalt, nact, nlay = atm.zlay.shape
    nlev = nlay+1

    output_atm = nc.Dataset(filename, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis
    output_atm.createDimension('location', 3)                 # source coordiantes
    output_atm.createDimension('emission', 1)                 # emission strength

    _dims3dlay = ('bins_along_track', 'bins_across_track', 'number_layers')
    _dims3dlev = ('bins_along_track', 'bins_across_track', 'number_levels')
    _dims2d    = ('bins_along_track', 'bins_across_track')

    out_shape_2d = (nalt, nact)

    gases = [
        x.removeprefix('dcol_') for x in atm.__dict__.keys()
        if 'dcol_' in x and x != 'dcol_air'
    ]

    # level height
    _ = writevariablefromname(output_atm, 'levelheight', _dims3dlev, atm.zlev)
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims3dlay, atm.zlay)
    # level pressure
    _ = writevariablefromname(
        output_atm,
        'pressure_layers',
        _dims3dlay,
        np.einsum("xy,z->xyz", np.ones(out_shape_2d), atm_std.play)
    )

    _ = writevariablefromname(
        output_atm,
        'pressure_levels',
        _dims3dlev,
        np.einsum("xy,z", np.ones(out_shape_2d), atm_std.plev)
    )

    _ = writevariablefromname(
        output_atm,
        'temperature',
        _dims3dlay,
        np.einsum("xy,z", np.ones(out_shape_2d), atm_std.tlay)
    )

    for gas in gases:
        # subcolumn density
        _ = writevariablefromname(output_atm, 'subcol_density_'+gas,
                                  _dims3dlay, atm.__getattribute__('dcol_'+gas))
        # column mixing ratio
        _ = writevariablefromname(output_atm, 'X'+gas, _dims2d,
                                  constants.__getattribute__('scale_X'+gas)*
                                  atm.__getattribute__('X'+gas))
        
#    if(config['profile']=='orbit'):
    _ = writevariablefromname(
        output_atm,
        'subcol_density_air',
        _dims3dlay,
        atm.dcol_air
    )

    # albedo
    for s2_albedo in atm.albedo:
        var_alb = writevariablefromname(
            output_atm,
            'albedo_' + s2_albedo.band_label,
            _dims2d,
            s2_albedo.values)
        var_alb.setncattr("central wavelength", s2_albedo.central_wavelength)
        var_alb.setncattr("band width", s2_albedo.bandwidth)

#    _ = writevariablefromname(output_atm, 'albedo', _dims2d, atm.albedo)
    # xpos and ypos
    _ = writevariablefromname(output_atm, 'xpos', 'bins_across_track', atm.xpos)
    _ = writevariablefromname(output_atm, 'ypos', 'bins_along_track', atm.ypos)

    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims2d, atm.col_air)
    # longitude/latitude coordiantes
    _ = writevariablefromname(output_atm, 'latitude', _dims2d, atm.lat)
    _ = writevariablefromname(output_atm, 'longitude', _dims2d, atm.lon)

    #information on emission source
    substr = 'source'
    attr = atm.__dict__.keys()
    attr_src = [string for string in attr if substr in string]
    for src in attr_src:
        _ = writevariablefromname(output_atm, src+'_location', 'location', atm.__getattribute__(src))

    substr = 'emission'
    attr = atm.__dict__.keys()
    attr_emi = [string for string in attr if substr in string]
    for emi in attr_emi:
        _ = writevariablefromname(output_atm, emi.removesuffix('_in_kgps'), 'emission' , atm.__getattribute__(emi))

    grp_aerosol = None
    # Aerosol
    if refractive_index is not None:
        grp_aerosol = output_atm.createGroup('aerosol')

        data_species = [
            refractive_index.dataBC,
            refractive_index.dataDUST,
            refractive_index.dataH2O,
            refractive_index.dataOC,
            refractive_index.dataIO,
            refractive_index.dataIO,
        ]

        for spi, spd in enumerate(data_species):
            spi1 = spi + 1
            grp_aerosol.createDimension(
                f'nwave_species{spi1}',
                len(spd[:,0])
            )
            _ = writevariablefromname(
                grp_aerosol,
                f'aer_wl_sp{spi1}',
                f'nwave_species{spi1}',
                np.array(spd[:,0])
            )
            _ = writevariablefromname(
                grp_aerosol,
                f'aer_rri_sp{spi1}',
                f'nwave_species{spi1}',
                np.array(spd[:,1])
            )
            _ = writevariablefromname(
                grp_aerosol,
                f'aer_iri_sp{spi1}',
                f'nwave_species{spi1}',
                np.array(spd[:,2])
            )

    if aerosol_echam is not None:
        if grp_aerosol is None:
            grp_aerosol = output_atm.createGroup('aerosol')

        grp_aerosol.createDimension('nmode_aer', 7)
        grp_aerosol.createDimension('nspecies_aer', 6)


        aerosol_echam.parameter_echam_orbit[0:7, :] = (
            1.0 - aerosol_echam.parameter_echam_orbit[28:35, :]
        )

        for i, var in enumerate(varname_aerosol_out):
            if i >= len(aerosol_echam.varname_echam) - 2:
                if (
                    i < (
                        len(aerosol_echam.varname_echam)
                        + len(aerosol_echam.varname_echam_part3) - 2
                    )
                ):
                    _ = writevariablefromname(
                        grp_aerosol,
                        "aer_" + var.lower(),
                        _dims2d,
                        aerosol_echam.parameter_echam_part3_orbit[
                            i-(len(aerosol_echam.varname_echam) - 2)
                        ][:].reshape(out_shape_2d)
                    )
                else:
                    par = (
                        aerosol_echam.parameter_echam_orbit[
                            i-len(aerosol_echam.varname_echam_part3)
                        ][:].reshape(out_shape_2d)
                    )
                    _ = writevariablefromname(
                        grp_aerosol,
                        "aer_" + var.lower(),
                        _dims2d,
                        par
                    )
                    if var == 'pressure_surface_aerosol':
                        surface_pressure_dem = np.array(par, copy=True)
                    elif var == 'elevation_surface_aerosol':
                        surface_elevation_dem = np.array(par, copy=True)
            else:
                _ = writevariablefromname(
                    grp_aerosol,
                    "aer_" + var.lower(),
                    _dims2d,
                    aerosol_echam.parameter_echam_orbit[i][:].reshape(
                        out_shape_2d
                    )
                )

    try:
        surface_pressure_dem
    except:
        surface_pressure_dem = 1008. * np.ones((nalt, nact))

    # Surface pressure and elevation from DEM are currently copied from ECHAM
    # aerosol data in the for-loop above. However, if a proper DEM is added,
    # the values written below should be taken from the DEM, whereas the ones
    # in the aerosol group should remain the same!
    _ = writevariablefromname(
        output_atm,
        'surface_pressure',
        _dims2d,
        surface_pressure_dem
    )
    _ = writevariablefromname(
        output_atm,
        'surface_elevation',
        _dims2d,
        #surface_elevation_dem
        np.zeros((nalt, nact))
    )

    grp_surface = None
    if xbrdf is not None:
        grp_surface = output_atm.createGroup("surface")
        _ = writevariablefromname(
            grp_surface,
            "surf_lisparse",
            _dims2d,
            xbrdf.xbrdf_rli_orbit[0,:].reshape(out_shape_2d)
        )
        _ = writevariablefromname(
            grp_surface,
            "surf_rossthick",
            _dims2d,
            xbrdf.xbrdf_rli_orbit[1,:].reshape(out_shape_2d)
        )
        _ = writevariablefromname(
            grp_surface,
            "surf_snow",
            _dims2d,
            np.zeros(out_shape_2d)
        )
        _ = writevariablefromname(
            grp_surface,
            "surf_water_frac",
            _dims2d,
            np.zeros(out_shape_2d)
        )
        _ = writevariablefromname(
            grp_surface,
            "surf_water_wind",
            _dims2d,
            7.0*np.ones(out_shape_2d)
        )

    if xbpdf is not None:
        if grp_surface is None:
            grp_surface = output_atm.createGroup("surface")

        _ = writevariablefromname(
            grp_surface,
            "surf_bpdf",
            _dims2d,
            xbpdf.xbpdf_orbit.reshape(out_shape_2d)
        )
    output_atm.close()

    return

def geosgm_output_ind_spec(filename, atm):
    # write geophysical scene data to output

    nalt, nact, nlay = atm.zlay.shape
    nlev = nlay+1

    output_atm = nc.Dataset(filename, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
    output_atm.createDimension('bins_along_track', nalt)      # along track axis
    output_atm.createDimension('bins_across_track', nact)     # across track axis
    output_atm.createDimension('number_layers', nlay)         # layer axis
    output_atm.createDimension('number_levels', nlev)         # level axis

    _dims3dlay = ('bins_along_track', 'bins_across_track', 'number_layers')
    _dims3dlev = ('bins_along_track', 'bins_across_track', 'number_levels')
    _dims2d    = ('bins_along_track', 'bins_across_track')

    gases = [x.removeprefix('dcol_') for x in atm.__dict__.keys() if 'dcol_' in x]

    # level height
    _ = writevariablefromname(output_atm, 'levelheight', _dims3dlev, atm.zlev)
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims3dlay, atm.zlay)


    for gas in gases:
        # subcolumn density
        _ = writevariablefromname(output_atm, 'subcol_density_'+gas,
                                  _dims3dlay, atm.__getattribute__('dcol_'+gas))
        # column mixing ratio
        _ = writevariablefromname(output_atm, 'X'+gas, _dims2d,
                                  constants.__getattribute__('scale_X'+gas)*
                                  atm.__getattribute__('X'+gas))
    # albedo
    _ = writevariablefromname(output_atm, 'albedo_B11',_dims2d,atm.albedo)

    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims2d, atm.col_air)
    # longitude/latitude coordiantes
    _ = writevariablefromname(output_atm, 'latitude', _dims2d, atm.lat)
    _ = writevariablefromname(output_atm, 'longitude', _dims2d, atm.lon)
    # xpos and ypos
    _ = writevariablefromname(output_atm, 'xpos', 'bins_across_track', atm.xpos)
    _ = writevariablefromname(output_atm, 'ypos', 'bins_along_track', atm.ypos)

    output_atm.close()

    return


def geoscene_generation(config: dict) -> None:
    """Generate a geophysical scene.

    Args:
      config
        Configuration dictionary

    """
    # first  the geometry data
    gm_data = get_gm_data(config['io_files']['input_gm'])
    nalt, nact = gm_data.sza.shape

    # =============================================================================
    # get a model atmosphere form AFGL files
    # =============================================================================
    nlay = config['atmosphere']['nlay']  # number of layers
    dzlay = config['atmosphere']['dzlay']
    # we assume the same standard atmosphere for all pixels of the granule

    atm_std = libATM.get_AFGL_atm_homogenous_distribution(
        config['io_files']['input_afgl'], nlay, dzlay, config['scale_gas']['xco2'],
        config['scale_gas']['xch4'], config['scale_gas']['xh2o'])

    # individual spectra and single swath
    if (config['profile'] == 'individual_spectra'):
        #use this profile whe you want to study the retrieval dependence as a
        #function of one parameter
        if(nalt!= 1):
            sys.exit("input error in sgm, for profile = indiudual spectra, nalt!=1")
        if (len(config['scene_spec']['albedo']) != nact):
            sys.exit("input error in sgm, albedo dimension not consistent with gm")

        alb = np.zeros([nalt, nact])
        alb[0, :] = config['scene_spec']['albedo'][:]
        albedo = DataArray(alb,
                           dims=('y', 'x'),)
        albedo.attrs['gsd'] = 300  #m 
        albedo.attrs['band_label'] = 'B11'
        albedo.attrs['central_wavelength'] = 1620
        albedo.attrs['bandwidth'] = '60'

        atm = libATM.create_atmosphere_ind_spectra(nalt, nact, atm_std, albedo, gm_data)

        geosgm_output(config['io_files']['output_geo'], atm)
    # Orbit
    if (config['profile'] == 'orbit'):

        # meteorological data

        meteodata = libATM.get_atmosphericdata_new(gm_data.lat, gm_data.lon, config['io_files']['meteo'])

        # Get albedo on the microHH grid
        s2_albedos = get_sentinel2_albedo(config['io_files']['input_s2'])

        #replace nan with closest non-nan value
        for s2_alb in s2_albedos:
            mask = np.isnan(s2_alb.values)
            idx = np.where(~mask,np.arange(mask.shape[1]),0)
            np.maximum.accumulate(idx,axis=1, out=idx)
            s2_alb.values[mask] = s2_alb.values[np.nonzero(mask)[0], idx[mask]]

        s2_albedos = libSGM.interp_sentinel2_albedo(
            s2_albedos,
            meteodata.lat,
            meteodata.lon,
            config['sentinel2']['band_label'])

        meteodata.__setattr__("albedo", s2_albedos)

        atm = libATM.combine_meteo_standard_atm_new(meteodata, atm_std, config)

        atm_lat_flat = atm.lat.reshape(-1)
        atm_lon_flat = atm.lon.reshape(-1)

        refractive_index = None
        aerosol_echam = None
        xbrdf = None
        xbpdf = None

        shape2d = atm.lat.shape
        surface_elevation_dem = atm_std.zlev[-1]*np.ones(shape2d)
        surface_pressure_dem = atm_std.plev[-1]*np.ones(shape2d)

        for module in [
            'refractive_index',
            'aerosol_echam',
            'collocation_algorithm',
            'xbrdf_gome2_modis',
            'xbpdf_polder'
        ]:
            logger = logging.getLogger(
                f'teds.lib.remotap_preproc.{module}'
            )
            logger.setLevel(logging.INFO)

        if (
            'aerosol' in config
            and 'add_aerosol' in config['aerosol']
            and config['aerosol']['add_aerosol']
        ):
            refractive_index = RefractiveIndex(
                config['io_files']['aerosol']['path_refr_index']
            )
            mode_coeff = refractive_index.get_coefficients_for_wavelength(550)
            aerosol_echam = AerosolEcham(
                config['io_files']['aerosol']['path_echam'],
                config['aerosol']['echam_year'],
                config['aerosol']['echam_month'],
                config['aerosol']['echam_day'],
                mode_coeff
            )
            jd = convert2julian7(
                config['aerosol']['echam_year'],
                config['aerosol']['echam_month'],
                config['aerosol']['echam_day'],
                0,
                0,
                0,
                0
            )
            jd_arr = jd*np.ones_like(atm_lat_flat)

            aerosol_echam.collocate(
                jd_arr,
                atm_lat_flat,
                atm_lon_flat,
                len(atm_lat_flat)
            )
            surface_elevation_dem = aerosol_echam.parameter_echam_orbit[
                varname_aerosol_out.index('elevation_surface_aerosol')
                -len(aerosol_echam.varname_echam_part3)
            ][:]
            surface_pressure_dem = aerosol_echam.parameter_echam_orbit[
                varname_aerosol_out.index('pressure_surface_aerosol')
                -len(aerosol_echam.varname_echam_part3)
            ][:]

        if (
            'surface' in config
            and 'add_brdf' in config['surface']
            and config['surface']['add_brdf']
        ):
            jd = convert2julian7(
                config['surface']['brdf_year'],
                config['surface']['brdf_month'],
                config['surface']['brdf_day'],
                0,
                0,
                0,
                0
            )
            jd_arr = jd*np.ones_like(atm_lat_flat)
            if not config['surface']['use_surface_classification']:
                xbrdf = XbrdfGomeModis(
                    config['io_files']['surface']['path_xbrdf'],
                    config['surface']['brdf_month'],
                    combined_with_sentinel=False
                )
                xbrdf.collocate(
                    jd_arr,
                    atm_lat_flat,
                    atm_lon_flat,
                    len(atm_lat_flat)
                )
            else:
                raise NotImplementedError(
                    "Surface classification coming soon!"
                )

            xbpdf = XbpdfPolder(
                config['io_files']['surface']['path_xbpdf'],
                '2006',
                config['surface']['brdf_month']
            )
            xbpdf.collocate(
                jd_arr,
                atm_lat_flat,
                atm_lon_flat,
                len(atm_lat_flat)
            )

        geosgm_output(
            config['io_files']['output_geo'],
            atm,
            atm_std,
            refractive_index,
            aerosol_echam,
            xbrdf,
            xbpdf
        )

    print('=>sgm geoscene calculation finished successfully')

if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    geoscene_generation(config)
