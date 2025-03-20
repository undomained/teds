# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
# =============================================================================
#     geophysical scene generation module for different E2E simulator profiles
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# =============================================================================

import netCDF4 as nc
import numpy as np
import numpy.typing as npt
import os

from ..lib import constants
from ..lib import libNumTools
from ..lib.libWrite import variable_dict
from ..lib.libWrite import writevariablefromname
from .sgm_Carbon_radscene import get_geosgm_data
from .sgm_Carbon_radscene import get_gm_data


def create_remotap_geometry(config: dict,
                            satellite_altitude: npt.NDArray[np.float64],
                            julian_day: npt.NDArray[np.float64],
                            saa: npt.NDArray[np.float64],
                            vaa: npt.NDArray[np.float64],
                            sza: npt.NDArray[np.float64],
                            vza: npt.NDArray[np.float64],
                            lat: npt.NDArray[np.float64],
                            lon: npt.NDArray[np.float64]) -> None:
    with nc.Dataset(config['io_files']['output_geometry'], 'w') as geom:
        nalt, nact = lat.shape
        npix = nalt*nact

        geom.createDimension('Nacross_track', nact)
        geom.createDimension('Nalong_track', nalt)
        geom.createDimension('Npix', npix)
        geom.createDimension('maxvza', 1)

        writevariablefromname(geom,
                              'remotap_altitude',
                              'Npix',
                              satellite_altitude.reshape(-1))
        writevariablefromname(geom,
                              'remotap_julday',
                              'Npix',
                              julian_day.reshape(-1))
        writevariablefromname(geom,
                              'remotap_latitude',
                              'Npix',
                              lat.reshape(-1))
        writevariablefromname(geom,
                              'remotap_longitude',
                              'Npix',
                              lon.reshape(-1))
        writevariablefromname(geom,
                              'remotap_solarazimuthangle',
                              ('Npix', 'maxvza'),
                              saa.reshape(-1, 1))
        writevariablefromname(geom,
                              'remotap_solarzenithangle',
                              ('Npix', 'maxvza'),
                              sza.reshape(-1, 1))
        writevariablefromname(geom,
                              'remotap_viewingazimuthangle',
                              ('Npix', 'maxvza'),
                              vaa.reshape(-1, 1))
        writevariablefromname(geom,
                              'remotap_viewingzenithangle',
                              ('Npix', 'maxvza'),
                              vza.reshape(-1, 1))
        writevariablefromname(geom,
                              'remotap_pixelid',
                              'Npix',
                              np.arange(npix)+1)


def create_remotap_input(config: dict) -> None:
    with nc.Dataset(config['io_files']['input_sgm_geo']) as in_nc, \
         nc.Dataset(config['io_files']['output_input'], 'w') as out_nc:
        npix = (
            in_nc.dimensions['bins_along_track'].size
            * in_nc.dimensions['bins_across_track'].size
        )
        out_nc.createDimension('Npix', npix)

        # Aerosol
        grp_aerosol = out_nc.createGroup('aerosol')
        for k, d in in_nc.groups['aerosol'].dimensions.items():
            grp_aerosol.createDimension(k, d.size)

        for k, v in in_nc.groups['aerosol'].variables.items():
            out_dict_k = None
            for dict_k, dict_v in variable_dict.items():
                if dict_v['name'] == k and dict_k.startswith('aer'):
                    out_dict_k = 'remotap_' + dict_k

            out_dims = v.dimensions
            out_data = v[()]
            if (
                out_dims[0] == 'bins_along_track'
                and out_dims[1] == 'bins_across_track'
            ):
                out_dims = ('Npix',) + v.dimensions[2:]
                out_data = out_data.reshape((npix,) + out_data.shape[2:])

            _ = writevariablefromname(
                grp_aerosol,
                out_dict_k,
                out_dims,
                out_data
            )

        # Surface
        grp_surface = out_nc.createGroup("surface")
        for k, d in in_nc.groups['surface'].dimensions.items():
            grp_surface.createDimension(k, d.size)

        for k, v in in_nc.groups['surface'].variables.items():
            out_dict_k = None
            for dict_k, dict_v in variable_dict.items():
                if dict_v['name'] == k and dict_k.startswith('surf'):
                    out_dict_k = 'remotap_' + dict_k

            out_dims = v.dimensions
            out_data = v[()]

            if (
                out_dims[0] == 'bins_along_track'
                and out_dims[1] == 'bins_across_track'
            ):
                out_dims = ('Npix',) + v.dimensions[2:]
                out_data = out_data.reshape((npix,) + out_data.shape[2:])

            _ = writevariablefromname(
                grp_surface,
                out_dict_k,
                out_dims,
                out_data
            )

        _alb_wls = []
        alb_data = []
        for k, v in in_nc.variables.items():
            if k.startswith('albedo'):
                _alb_wls.append(v.getncattr('central wavelength'))
                alb_data.append(v[()].reshape(-1))

        alb_wls = np.array(_alb_wls)
        alb_data = np.stack(alb_data)
        alb_order = np.argsort(alb_wls)
        alb_wls = alb_wls[alb_order]
        alb_data = alb_data[alb_order]
        grp_surface.createDimension('Nbands_combine', len(alb_wls))
        _ = writevariablefromname(
            grp_surface,
            'remotap_surf_bands_combine',
            'Nbands_combine',
            alb_wls
        )
        _ = writevariablefromname(
            grp_surface,
            'remotap_surf_xbrdf_wave_combine',
            ('Nbands_combine', 'Npix'),
            alb_data
        )

        # Ocean (dummy for now)
        grp_ocean = out_nc.createGroup('ocean')
        _ = writevariablefromname(
            grp_ocean,
            'remotap_ocn_Xchl',
            'Npix',
            -999*np.ones(npix)
        )
        for par in ('uwnd', 'vwnd'):
            _ = writevariablefromname(
                grp_ocean,
                f'remotap_ocn_{par}',
                'Npix',
                np.zeros(npix)
            )

        # Cirrus (dummy for now)
        grp_cirrus = out_nc.createGroup('cirrus')
        for par in (
            'aspect_ratio',
            'cgd_mean',
            'cod_median',
            'cth_mean',
            'effective_radius',
            'roughness'
        ):
            _ = writevariablefromname(
                grp_cirrus,
                f'remotap_cir_{par}',
                'Npix',
                9.96921e36*np.ones(npix)
            )
        _ = writevariablefromname(
            grp_cirrus,
            'remotap_cir_cirrus_fraction',
            'Npix',
            np.zeros(npix)
        )

        # Cloud (dummy for now)
        grp_cloud = out_nc.createGroup('cloud')
        # Aerosol species 3 is water
        rri_water = in_nc.groups['aerosol'].variables['RRI species 3'][:]
        iri_water = in_nc.groups['aerosol'].variables['IRI species 3'][:]
        wl_water = in_nc.groups['aerosol'].variables['wavelength species 3'][:]
        nwave_water = len(wl_water)

        grp_cloud.createDimension('nwave_water_refr', nwave_water)
        _ = writevariablefromname(
            grp_cloud,
            'remotap_cld_wavelength_water_refr',
            'nwave_water_refr',
            wl_water
        )
        _ = writevariablefromname(
            grp_cloud,
            'remotap_cld_RRI_water_refr',
            'nwave_water_refr',
            rri_water
        )
        _ = writevariablefromname(
            grp_cloud,
            'remotap_cld_IRI_water_refr',
            'nwave_water_refr',
            iri_water
        )

        for par in (
            'cloud_reff',
            'cloud_optical_thickness',
            'cloud_top_height_asl',
            'cloud_phase'
        ):
            _ = writevariablefromname(
                grp_cloud,
                f'remotap_cld_{par}',
                'Npix',
                9.96921e36*np.ones(npix)
            )

        _ = writevariablefromname(
            grp_cloud,
            'remotap_cld_cloud_fraction',
            'Npix',
            np.zeros(npix)
        )


def create_remotap_aux(config: dict) -> None:
    with nc.Dataset(config['io_files']['input_sgm_geo']) as in_nc, \
         nc.Dataset(config['io_files']['output_aux'], 'w') as out_nc:
        npix = (
            in_nc.dimensions['bins_along_track'].size
            * in_nc.dimensions['bins_across_track'].size
        )
        out_nc.createDimension('Npix', npix)

        nlay = in_nc.dimensions['number_layers'].size
        nlev = in_nc.dimensions['number_levels'].size

        out_nc.createDimension('nlay', nlay)
        out_nc.createDimension('nlev', nlev)

        temperature = in_nc.variables['temperature'][:].reshape(npix, nlay)
        zlay = in_nc.variables['zlay'][:].reshape(npix, nlay)
        temp_grad = (
            (temperature[:, 1:]-temperature[:, 0:-1])
            / (zlay[:, 1:]-zlay[:, 0:-1])
        )
        ind5k = np.argmin(np.abs(zlay-5000), axis=1)
        ind20k = np.argmin(np.abs(zlay-20000), axis=1)
        for ipix in range(npix):
            temp_grad[ipix, ind5k[ipix]:] = 1e36
            temp_grad[ipix, :ind20k[ipix]] = 1e36

        trop_ind = np.argmin(np.abs(temp_grad), axis=1)
        trop_h = np.minimum(
            20000,
            np.maximum(5000, zlay[np.arange(npix), trop_ind])
        )

        _ = writevariablefromname(
            out_nc,
            'remotap_aux_pixel_id',
            'Npix',
            np.arange(npix) + 1
        )
        for par in (
            'cirrus_fraction',
            'cloud_top_height_asl',
            'cloud_fraction'
        ):
            _ = writevariablefromname(
                out_nc,
                f'remotap_aux_{par}',
                'Npix',
                np.zeros(npix)
            )
        _ = writevariablefromname(
            out_nc,
            'remotap_aux_land_fraction',
            'Npix',
            np.ones(npix)
        )
        for src in ('dem', 'meteo'):
            _ = writevariablefromname(
                out_nc,
                f'remotap_aux_surface_elevation_{src}',
                'Npix',
                in_nc.variables['surface elevation'][:].reshape(npix)
            )
        _ = writevariablefromname(
            out_nc,
            'remotap_aux_surface_elevation_std',
            'Npix',
            np.zeros(npix)
        )
        _ = writevariablefromname(
            out_nc,
            'remotap_aux_zlev',
            ('Npix', 'nlev'),
            in_nc.variables['zlev'][:].reshape(npix, nlev)
        )
        _ = writevariablefromname(
            out_nc,
            'remotap_aux_plev',
            ('Npix', 'nlev'),
            in_nc.variables['pressure_levels'][:].reshape(npix, nlev)
        )

        _ = writevariablefromname(
            out_nc,
            'remotap_aux_temperature',
            ('Npix', 'nlay'),
            temperature
        )
        _ = writevariablefromname(
            out_nc,
            'remotap_aux_tropopause_height',
            'Npix',
            trop_h
        )

        _ = writevariablefromname(
            out_nc,
            'remotap_aux_surface_pressure_meteo',
            'Npix',
            in_nc.variables['surface pressure'][:].reshape(npix)
        )

        for ab in 'ab':
            _ = writevariablefromname(
                out_nc,
                f'remotap_aux_hybrid_{ab}',
                'nlev',
                np.zeros(nlev)
            )

        dcol_air = in_nc["dcol_air"][:]  # molec m^-2
        for sp in ('ch4', 'co2', 'no2', 'o3'):
            dcol = in_nc[f'dcol_{sp}'][:]  # molec m^-2
            vmr = dcol/dcol_air
            mass = getattr(constants, 'M'+sp.upper())  # kg mol^-1
            mmr = vmr * mass / constants.MDRYAIR

            _ = writevariablefromname(
                out_nc,
                f'remotap_aux_MMR_{sp}',
                ('Npix', 'nlay'),
                mmr.reshape(npix, nlay)
            )

        dcol = in_nc['dcol_h2o'][:]
        _ = writevariablefromname(
            out_nc,
            'remotap_aux_specific_humidity',
            ('Npix', 'nlay'),
            (dcol/dcol_air).reshape(npix, nlay)
        )


def fill_ini_template(datadir: str,
                      geomfile: str,
                      inputfile: str,
                      auxfile: str,
                      resultdir: str) -> str:
    ini_template = """
* RemoTAP configuration (Berlin regional scenario, CO2I)
********************************************************************************
1 1
* Definition of the spectral bands to simulate.
1
1 2 1590.00 1675.0 0 3 3 1.32e-7 3 2.025e5    #1:[photon/m2/s/microns/sr]; \
    2:[mol/m2/s/nm/sr; 3:photons/cm2/s/nm/sr
  1 1 1.0 0.0 2.0 1.0   !mol hitran id(e: water vapor id=1), (e: from here to \
    end just copy the values)fit flag, apri_err, par_min, par_max, prior
  {datadir}/XSDB/joostadb/LUT/hit08_01_Voigt_5900.00-6350.00_spec_res_0.005.nc
  2 1 1.0 0.0 2.0 1.0   !mol hitran id(e: co2 id=2), fit flag, apri_err, \
    par_min, par_max, prior
  {datadir}/XSDB/joostadb/LUT/hit08_02_LM__full_5900.00-6350.00_spec_res_0.005.nc
  6 1 1.0 0.0 2.0 1.0   !mol hitran id (e: methan id=6), fit flag, apri_err, \
    par_min, par_max, prior
  {datadir}/XSDB/joostadb/LUT/hit08_06_Voigt_5900.00-6350.00_spec_res_0.005.nc
  0.00265 10 10 0.3 0.085945399    !reso, wvbd,ntau,fwhm,reso_meas      ! \
    reso (0.00265, 0.01)
  3 3 1 3
********************************************************************************
* File containing observation geometry (geo-location, solar and viewing \
    angles).
{geomfile}
********************************************************************************
* Surface BPDF configuration.
M      # F = Facet model, M = Maignan model
********************************************************************************
* Surface BRDF configuration.
KER    # KER = Kernel-based, RPV = Rahman-Pinty-Verstraete model
0      # Add Fresnel term to R(1,1) => 0 = NO, 1 = YES
********************************************************************************
* File containing the input parameters.
{inputfile}
* File containing the auxiliary parameters.
{auxfile}
********************************************************************************
* Path to the location where output files will be written.
{resultdir}
********************************************************************************
* Path to Mie/T-Matrix LUTs.
{datadir}/MIE_TMATRIX/MIE_TMATRIX_created20230525/
********************************************************************************
* Path to miscellaneous small data files.
{datadir}/SMALL_DATA/SMALL_DATA_created20230525/
********************************************************************************
* Path to underlight LUTs.
{datadir}/OCEAN_BODY/NN_TABLE_created20230525/
********************************************************************************
* Path to solar irradiance LUT.
{datadir}/SOLAR/irrcmasun2000_co2m_created20230525.nc
********************************************************************************
* Path to cirrus LUTs.
{datadir}/CIRRUS/phase_matrix_cirrus_full_created20230525.nc
********************************************************************************
"""
    return ini_template.format(
        datadir=datadir,
        geomfile=geomfile,
        inputfile=inputfile,
        auxfile=auxfile,
        resultdir=resultdir
    )


def create_remotap_ini(config: dict) -> None:
    if not config['io_files']['dir_remotap_data'].endswith('/'):
        config['io_files']['dir_remotap_data'] += '/'
    if not config['io_files']['dir_remotap_result'].endswith('/'):
        config['io_files']['dir_remotap_result'] += '/'

    ini_string = fill_ini_template(
        datadir=config['io_files']['dir_remotap_data'],
        geomfile=config['io_files']['output_geometry'],
        inputfile=config['io_files']['output_input'],
        auxfile=config['io_files']['output_aux'],
        resultdir=config['io_files']['dir_remotap_result']
    )

    if not os.path.isdir(config['io_files']['dir_remotap_result']):
        os.makedirs(config['io_files']['dir_remotap_result'])
        print(
            'Created output directory',
            config['io_files']['dir_remotap_result']
        )

    with open(config['io_files']['output_ini'], "w") as f:
        print(ini_string, file=f)


class TransformCoords:
    """A class to transform coordinates"""

    def __init__(self, origin):
        """Initialize class based on the Origin

        Parameters
        ----------
        origin : [lat, lon] Float64
            Origin in lat-lon coordinates
        """
        self.rd = np.pi / 180.0
        phi0 = origin[0] * self.rd
        self.ld0 = origin[1] * self.rd
        self.s_p0 = np.sin(phi0)
        self.c_p0 = np.cos(phi0)
        self.fact = 6371.0
        self.factmts = self.fact * 1000

    def latlon2xykm(self, lat, lon):
        """latlon2xykm Convert from lat lon to xy in km

        Parameters
        ----------
        lat : Matrix/Vector/Float64
            Latitude
        lon : Matrix/Vector/Float64
            Longitude

        Returns
        -------
        x  : Matrix/Vector/Float64
            x in km
        y  : Matrix/Vector/Float64
            y in km
        """
        ld = lon * self.rd
        phi = lat * self.rd
        s_p = np.sin(phi)
        c_p = np.cos(phi)
        ll = ld - self.ld0
        c_l = np.cos(ll)
        s_l = np.sin(ll)
        c_pl = c_p * c_l
        w = np.sqrt(2.0 / (np.maximum(
            1.0 + self.s_p0 * s_p + self.c_p0 * c_pl, 1.0e-10)))
        x = c_p * s_l * w
        y = (self.c_p0 * s_p - self.s_p0 * c_pl) * w
        return x * self.fact, y * self.fact

    def latlon2xymts(self, lat, lon):
        """latlon2xymts Convert from lat lon to xy in km

        Parameters
        ----------
        lat : Matrix/Vector/Float64
            Latitude
        lon : Matrix/Vector/Float64
            Longitude

        Returns
        -------
        x  : Matrix/Vector/Float64
            x in mts
        y  : Matrix/Vector/Float64
            y in mts
        """
        ld = lon * self.rd
        phi = lat * self.rd
        s_p = np.sin(phi)
        c_p = np.cos(phi)
        ll = ld - self.ld0
        c_l = np.cos(ll)
        s_l = np.sin(ll)
        c_pl = c_p * c_l
        w = np.sqrt(2.0 / (np.maximum(
            1.0 + self.s_p0 * s_p + self.c_p0 * c_pl, 1.0e-10)))
        x = c_p * s_l * w
        y = (self.c_p0 * s_p - self.s_p0 * c_pl) * w
        return x * self.factmts, y * self.factmts

    def xykm2latlon(self, x1, y1):
        """xykm2latlon Convert from x, y to lat-lon

        Parameters
        ----------
        x1 : Matrix/Vector/Float64
            x coordinate in km
        y1 : Matrix/Vector/Float64
            y coordinate in km

        Returns
        -------
        lat  : Matrix/Vector/Float64
            latitude
        lon  : Matrix/Vector/Float64
            Longitude
        """
        x, y = x1 / self.fact, y1 / self.fact
        p = np.maximum(np.sqrt(x**2 + y**2), 1.0e-10)
        c = 2.0 * np.arcsin(p / 2.0)
        s_c = np.sin(c)
        c_c = np.cos(c)
        phi = np.arcsin(c_c * self.s_p0 + y * s_c * self.c_p0 / p)
        ld = self.ld0 + np.arctan2(x * s_c,
                                   (p * self.c_p0 * c_c - y * self.s_p0 * s_c))
        lat = phi / self.rd
        lon = ld / self.rd
        if isinstance(lat, np.ndarray):
            lat[lat > 90.0] -= 180.0
            lat[lat < -90.0] += 180.0
            lon[lon > 180.0] -= 360.0
            lon[lon < -180.0] += 360.0
        else:
            if abs(lat) > 90.0:
                if lat > 0:
                    lat = lat - 180.0
                else:
                    lat = lat + 180.0
            if abs(lon) > 180.0:
                if lon > 0:
                    lon = lon - 360.0
                else:
                    lon = lon + 360.0
        return lat, lon

    def xymts2latlon(self, x1, y1):
        """xymts2latlon Convert from x, y to lat-lon

        Parameters
        ----------
        x1 : Matrix/Vector/Float64
            x coordinate in m
        y1 : Matrix/Vector/Float64
            y coordinate in m

        Returns
        -------
        lat  : Matrix/Vector/Float64
            latitude
        lon  : Matrix/Vector/Float64
            Longitude
        """
        x, y = x1 / self.factmts, y1 / self.factmts
        p = np.maximum(np.sqrt(x**2 + y**2), 1.0e-10)
        c = 2.0 * np.arcsin(p / 2.0)
        s_c = np.sin(c)
        c_c = np.cos(c)
        phi = np.arcsin(c_c * self.s_p0 + y * s_c * self.c_p0 / p)
        ld = self.ld0 + np.arctan2(x * s_c,
                                   (p * self.c_p0 * c_c - y * self.s_p0 * s_c))
        lat = phi / self.rd
        lon = ld / self.rd
        if isinstance(lat, np.ndarray):
            lat[lat > 90.0] -= 180.0
            lat[lat < -90.0] += 180.0
            lon[lon > 180.0] -= 360.0
            lon[lon < -180.0] += 360.0
        else:
            if abs(lat) > 90.0:
                if lat > 0:
                    lat = lat - 180.0
                else:
                    lat = lat + 180.0
            if abs(lon) > 180.0:
                if lon > 0:
                    lon = lon - 360.0
                else:
                    lon = lon + 360.0
        return lat, lon


def convert_geoscene_to_remotap(config: dict) -> None:
    """Convert Tango SGM geoscene to RemoTAP input file format."""
    # get data
    atm_org = get_geosgm_data(config['io_files']['input_sgm_geo'])
    # get gm data
    gm_org = get_gm_data(config['io_files']['input_gm'])
    # create a transform method
    trans = TransformCoords(atm_org.co2_src_zlatlon[1:])
    # convert lat-lon of gm to x-y and get bounds
    gm_org.xpos, gm_org.ypos = trans.latlon2xymts(gm_org.lat, gm_org.lon)
    # The orginal gm data for SZA, SAA. VZA, VAA are extrapolated to
    # the atmospheric mesh.
    gm_data = libNumTools.expand_geometry(atm_org, gm_org)

    print('Generating RemoTAP geometry... ', end='', flush=True)
    create_remotap_geometry(config,
                            gm_data.satellite_altitude,
                            gm_data.julian_day,
                            gm_data.saa,
                            gm_data.vaa,
                            gm_data.sza,
                            gm_data.vza,
                            atm_org.lat,
                            atm_org.lon)
    print('done.')
    print('Generating RemoTAP input... ', end='', flush=True)
    create_remotap_input(config)
    print('done.')
    print('Generating RemoTAP aux... ', end='', flush=True)
    create_remotap_aux(config)
    print('done.')
    print('Generating RemoTAP ini... ', end='', flush=True)
    create_remotap_ini(config)
    print('done.')
