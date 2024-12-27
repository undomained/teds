# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

from netCDF4 import Dataset
import os
import pickle
import sys
import numpy as np
import yaml
from tqdm import tqdm
from copy import deepcopy

from ..lib import libNumTools, libRT, libSURF, libATM
from ..lib.libWrite import writevariablefromname
from ..lib.libNumTools import TransformCoords, convert
from scipy.interpolate import RegularGridInterpolator,griddata

class Emptyclass:
    """Empty class. Data container."""

    pass


class Dict2Class:
    """Convert a dictionaly to a class."""

    def __init__(self, arg_dict):
        self.__dict__.update(arg_dict)

def get_gm_data(filename):

    names = [
        'sza', 'saa', 'vza', 'vaa', 'lat', 'lon',
    ]

    input = Dataset(filename, mode='r')

    gm_data = Emptyclass()

    for name in names:
        gm_data.__setattr__(name, input[name][:])

    input.close()

    return gm_data

def radsgm_output(filename_rad, rad_output, gm_data):
    # write radiances
    nalt, nact, nwav = rad_output['radiance'].shape
    # open file
    nc = Dataset(filename_rad, mode='w')
    nc.title = 'Tango Carbon E2ES SGM radiometric scene'
    nc.product_type = 'SGM'
    nc.createDimension('wavelength', nwav)     # spectral axis
    nc.createDimension('across_track_sample', nact)     # across track axis
    nc.createDimension('along_track_sample', nalt)     # along track axis

    # wavelength
    wavelength = np.zeros((nact, nwav))
    for i in range(nact):
        wavelength[i, :] = rad_output['wavelength']
    _ = writevariablefromname(nc,
                              "wavelength",
                              ('wavelength'),
                              rad_output['wavelength'])
    # solar irradiance
    _ = writevariablefromname(nc,
                              "solarirradiance",
                              ('wavelength'),
                              rad_output['solar irradiance'])
    # radiance

    _ = writevariablefromname(
        nc,
        'radiance_sgm',
        ('along_track_sample', 'across_track_sample', 'wavelength'),
        rad_output['radiance'])

    # add coordinates to SGM atmosphere
    _ = writevariablefromname(
        nc,
        'latitude',
        ('along_track_sample', 'across_track_sample'),
        gm_data.lat)

    _ = writevariablefromname(
        nc,
        'longitude',
        ('along_track_sample', 'across_track_sample'),
        gm_data.lon)

    return

def sgm_output_atm_ref(filename, atm, albedo, gm_data, gases):

    #define dimensions

    dim_alt,dim_act, dim_lay = atm.zlay.shape
    dim_lev = atm.zlev.shape[2]

    output_atm = Dataset(filename, mode='w')
    output_atm.title = 'Tango Carbon E2ES SGM atmospheric scene'
#    output_atm.createDimension('along_gm_org.__getattribute__(gm_para)track', dim_alt)      # along track axis
    output_atm.createDimension('along_track_sample', dim_alt)      # along track axis
    output_atm.createDimension('across_track_sample', dim_act)     # across track axis
    output_atm.createDimension('number_layers', dim_lay)         # layer axis
    output_atm.createDimension('number_levels', dim_lev)         # level axis
    _dims = ('along_track_sample', 'across_track_sample', 'number_layers')
    # central layer height
    _ = writevariablefromname(output_atm, 'central_layer_height', _dims, atm.zlay)
    # columndensity_co2
    _ = writevariablefromname(output_atm, 'subcol_density_co2', _dims, atm.CO2)
    # columndensity_ch4
    _ = writevariablefromname(output_atm, 'subcol_density_ch4', _dims, atm.CH4)
    # columndensity_h2o
    _ = writevariablefromname(output_atm, 'subcol_density_h2o', _dims, atm.H2O)
    # level height
    _dims = ('along_track_sample', 'across_track_sample', 'number_levels')
    _ = writevariablefromname(output_atm, 'levelheight', _dims, atm.zlev)

    xco2 = np.sum(atm.CO2,axis=2)/atm.air*1.e6  #[ppm]
    xch4 = np.sum(atm.CH4,axis=2)/atm.air*1.e9  #[ppb]
    xh2o = np.sum(atm.H2O,axis=2)/atm.air*1.e6  #[ppm]

    _dims = ('along_track_sample', 'across_track_sample')
    # albedo
    _ = writevariablefromname(output_atm, 'albedo', _dims, albedo)
    # column_co2
    _ = writevariablefromname(output_atm, 'XCO2', _dims, xco2)
    # column_ch4
    _ = writevariablefromname(output_atm, 'XCH4', _dims, xch4)
    # column_h2o
    _ = writevariablefromname(output_atm, 'XH2O', _dims, xh2o)
    # column_air
    _ = writevariablefromname(output_atm, 'column_air', _dims, atm.air)

    # add coordinates to SGM atmosphere
    _ = writevariablefromname(output_atm, 'latitude', _dims, gm_data.lat)
    _ = writevariablefromname(output_atm, 'longitude', _dims, gm_data.lon)
    output_atm.close()

def get_geosgm_data(filename):

    input = Dataset(filename, mode='r')

    names = ['col_air', 'dcol_ch4', 'dcol_co2', 'dcol_h2o', 'lat' , 'lon',
             'XCH4', 'XCO2', 'XH2O', 'zlay', 'zlev','xpos', 'ypos']

    atm_data = Emptyclass()

    for name in names:
        atm_data.__setattr__(name, input[name][:])

    atm_data.__setattr__('albedo', input['albedo B11'][:])

    sources = [x for x in input.variables.keys() if 'source location' in x]
    emissions = [x for x in input.variables.keys() if 'emission' in x]

    for source in sources:
        attrib = source.removesuffix(' source location').lower()+'_src_zlatlon'
        atm_data.__setattr__(attrib, input[source][:])

    for emission in emissions:
        attrib = emission.removesuffix(' emission').lower()+'_src_kgps'
        atm_data.__setattr__(attrib, input[emission][:])

    input.close()

    return atm_data

def extract_atm(atm,ialt,iact):

    atm_ext = Emptyclass()

    atm_ext.__setattr__('lat',atm.lat[ialt,iact])
    atm_ext.__setattr__('lon',atm.lon[ialt,iact])
    atm_ext.__setattr__('CO2',atm.CO2[ialt,iact,:])
    atm_ext.__setattr__('CH4',atm.CH4[ialt,iact,:])
    atm_ext.__setattr__('H2O',atm.H2O[ialt,iact,:])
    atm_ext.__setattr__('zlay',atm.zlay[ialt,iact,:])
    atm_ext.__setattr__('zlev',atm.zlev[ialt,iact,:])
    atm_ext.__setattr__('air',atm.air[ialt,iact])

    return atm_ext

def Carbon_radiation_scene_generation(config: dict) -> None:
    """Generate TOA spectra for the Tango-Carbon instrument.

    Args:
      config
        Configuration dictionary

    # The function is setup in two parts:
    # First the objects alebdo and atm are defined. There are three different ways to do this
    # 1. We take the geometry from the input file and convolve and intepolate the geosscene data 
    # to this mesh (runset[0] = 0)
    # 2. We take the atmosphere one-by-one from the geoscene file and interpolate and expand the 
    # geometry this mesh (runset[0] = 1)
    # 3. In the case that gm information and geo-scene is given on the same grid, just use them 
    # both as is (runset[0] = 2)
    # Second, we perform the radiative transfer simualtionon the grid provided  by alebdo, atm. 
    # This time, we have two options:
    # 1. line-by-line spectra are provided in the output (runset[1] = 0)
    # 2. To save memory, we convolve the spectra with the ISRF and sample the data with high 
    # oversampling (runset[1] = 1)
    # The different setting for the SGM-RAD simualtions are defined per configuration. 
    """
    if(config['convolve_atm_input'] & (config['profile']=='orbit')):
        runset = [0,0]
    if((not config['convolve_atm_input']) & (config['profile']=='orbit')):
        runset = [1,1]
    if((not config['convolve_atm_input']) & (config['profile']=='individual_spectra')):
        runset = [2,0]
    # Other configurations are possible are not considered as baseline so far. 

    #get data
    atm_org = get_geosgm_data(config['io_files']['input_sgm_geo'])
    #get gm data
    gm_org  = get_gm_data(config['io_files']['input_gm'])
    # create a transform method

    trans = TransformCoords(atm_org.co2_src_zlatlon[1:])
    # convert lat-lon of gm to x-y and get bounds
    gm_org.xpos, gm_org.ypos = trans.latlon2xymts(gm_org.lat, gm_org.lon)

    #prepare the input data to calculate the radiation scene

    if(runset[0]==0):

        # This option takes the gm data from gm_org and convolves the atmospheric
        # and albedo data with the SEDF. It provide one atmosphere for each observation grid point
        # It has the advantage to have a well-defined reference scene for each observation.
        # Its disaadvantage is that it cannot be used to study sub-pixel effects

        gm_data = gm_org
        #convolution with instrument spatial response
        atm_conv = libNumTools.convolvedata(atm_org, config)
        #interpolate convolved data to gm grid
        albedo, atm = libNumTools.interpolate_data_regular(atm_conv, gm_data, config["selected_gases"])
        #store sgm data that are used for RT simulations
        sgm_output_atm_ref(config['io_files']['output_geo_ref'], atm, albedo, gm_data, config["selected_gases"])

    if(runset[0]==1):

        # This option uses the atmospheric grid given by atm.org and extrapolates the geometry given
        # in gm_org to the atmospheric mesh

        #first we take over and reformate the model atmosphere form atm_org for the selected gases.
        albedo, atm = convert(atm_org, config["selected_gases"])

        #The orginal gm data for SZA, SAA. VZA, VAA are extrapolated to the atmospheric mesh
        gm_data = libNumTools.expand_geometry(atm_org, gm_org) 

    if(runset[0]==2):

        # This option uses the atmospheric grid given by atm.org and the geometry given
        # in gm_org assumning that they are given on the same mesh

        #first we take over the model atmosphere form atm_org for the selected gases.

        atm = Emptyclass()
        albedo = deepcopy(atm_org.albedo)
        for gas in config["selected_gases"]:
            atm.__setattr__(gas.upper(), atm_org.__getattribute__("dcol_"+gas))
        atm.__setattr__("lat",  atm_org.__getattribute__("lat"))
        atm.__setattr__("lon",  atm_org.__getattribute__("lon"))
        atm.__setattr__("zlay", atm_org.__getattribute__("zlay"))
        atm.__setattr__("air",  atm_org.__getattribute__("col_air"))
        atm.__setattr__("zlev", atm_org.__getattribute__("zlev"))
     
        #The orginal gm data for SZA, SAA. VZA, VAA are extrapolated to the atmospheric mesh
        gm_data = gm_org
        #store sgm data that are used for RT simulations
        sgm_output_atm_ref(config['io_files']['output_geo_ref'], atm, albedo, gm_data, config["selected_gases"])

        #fig, (ax0, ax1) = plt.subplots(2, 1)
        #im1 = ax0.pcolormesh(gm_org.ypos,gm_org.xpos,gm_org.vza)
        #fig.colorbar(im1, ax=ax0)
        #im2 = ax0.pcolormesh(gm_data.ypos,gm_data.xpos,gm_data,vza)
        #fig.colorbar(im2, ax=ax1)

    # Independently from config['convolve_atm_input'], we have atmospheric data and geometry given on the same grid.
    # for the dataset atm, albedo, gm_data, we perform the RT simulation

    # line-by-line spectral grid
    wave_start = config['spec_lbl_settings']['wave_start']
    wave_end   = config['spec_lbl_settings']['wave_end']
    dwave_lbl  = config['spec_lbl_settings']['dwave']
    wave_lbl   = np.arange(wave_start, wave_end, dwave_lbl)  # nm

    # define reference model atmosphere for which we calculate X sections
    # First, vertical layering
    nlay = config['atmosphere']['nlay']
    dzlay = config['atmosphere']['dzlay']
    psurf = config['atmosphere']['psurf']
    nlev = nlay + 1  # number of levels

    zlay = (np.arange(nlay-1, -1, -1)+0.5)*dzlay  # altitude of layer midpoint
    zlev = np.arange(nlev-1, -1, -1)*dzlay  # altitude of layer interfaces = levels

    # Second, online calculation of cross section or retrieving data from dummy file depending on config
    if ((not os.path.exists(config['io_files']['dump_xsec'])) or config['xsec_forced']):

        iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
        molec = libRT.molecular_data(wave_lbl)

        molec.get_data_HITRAN(config['io_files']['input_hapi'], iso_ids)

        atm_ref = libATM.atmosphere_data(zlay, zlev, psurf)
        atm_ref.get_data_AFGL(config['io_files']['input_afgl'])

        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm_ref)
        # Dump optics.prop dictionary into temporary pkl file
        pickle.dump(optics.prop, open(config['io_files']['dump_xsec'], 'wb'))

    else:

        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)

        # Read optics.prop dictionary from pickle file
        optics.prop = pickle.load(open(config['io_files']['dump_xsec'], 'rb'))

    # Next we perform the RT simuations for the generated input.
    # Two differnent cases depending on config['convolve_atm_input']
    # (1) For spatially convolved scenes, spectra are provided on the line-by-line spectral grid
    # (2) For scenes that are not spatially convolved, we provide spectra convolved with the ISRF
    #
    rad_output = {}

    # solar irradiance spectrum
    sun = libRT.read_sun_spectrum_TSIS1HSRS(config['io_files']['input_sun_reference'])
    sun_lbl= np.interp(wave_lbl, sun['wl'], sun['phsm2nm'])

    # dimensions
    nalt, nact = gm_data.sza.shape
    nwav_lbl = wave_lbl.size

    if(runset[1]==0):
        #line-by-line spectra stored in the output file
        rad_output['wavelength'] = wave_lbl  # nm
        rad_output['solar irradiance'] = sun_lbl

        # Calculate surface data
        surface = libSURF.surface_prop(wave_lbl)
        rad = np.empty([nalt, nact, nwav_lbl])

        print('Radiative tranfer simulation...')
        for ialt in tqdm(range(nalt)):
            for iact in range(nact):
                #extract the atmosphere that belongs to (ialt, iact)
                atm_ext = extract_atm(atm,ialt,iact)
                #calculate optical depths
                optics.set_opt_depth_species(atm_ext, ['molec_01', 'molec_32', 'molec_07'])
                # Earth radiance spectra
                alb = [albedo[ialt, iact]]
                mu_sza = np.cos(np.deg2rad(gm_data.sza[ialt, iact]))
                mu_vza = np.cos(np.deg2rad(gm_data.vza[ialt, iact]))
                surface.get_albedo_poly(alb)
                rad[ialt,iact,:] = libRT.transmission(
                    sun_lbl, optics, surface, mu_sza, mu_vza)
    if(runset[1]==1):
        #line-by-line simulations with subsquently ISRF convolution (line-by-line spectra 
        # are not kept in momory) 
        wave_start = config['spec_conv_settings']['wavestart']
        wave_end   = config['spec_conv_settings']['waveend']
        dwave      = config['spec_conv_settings']['dwave']
        wave_conv  = np.arange(wave_start, wave_end, dwave)  # nm
        nwav_conv  = wave_conv.size

        rad_output['wavelength'] = wave_conv

        #function for ISRF convolution
        isrf_convolution = libNumTools.get_isrf(wave_conv, wave_lbl, config['isrf_settings'])
        #convolved solar spectrum
        rad_output['solar irradiance'] = isrf_convolution(sun_lbl)

        # Calculate surface data
        surface = libSURF.surface_prop(wave_lbl)

        rad = np.empty([nalt, nact, nwav_conv])

        print('Radiative tranfer simulation...')
        for ialt in tqdm(range(nalt)):
            for iact in range(nact):
                atm_ext = extract_atm(atm,ialt,iact)
                optics.set_opt_depth_species(atm_ext, ['molec_01', 'molec_32', 'molec_07'])
                # Earth radiance spectra
                alb = [albedo[ialt, iact]]
                mu_sza = np.cos(np.deg2rad(gm_data.sza[ialt, iact]))
                mu_vza = np.cos(np.deg2rad(gm_data.vza[ialt, iact]))
                surface.get_albedo_poly(alb)
                radiance = libRT.transmission(
                        sun_lbl, optics, surface, mu_sza, mu_vza)

                rad[ialt,iact,:]= isrf_convolution(radiance)

    rad_output['radiance'] = rad

    # sgm output to radiometric file
    radsgm_output(config['io_files']['output_rad'], rad_output, gm_data)

    print('=>Carbon radsgm calculation finished successfully')
    return

if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    Carbon_radiation_scene_generation(config)
