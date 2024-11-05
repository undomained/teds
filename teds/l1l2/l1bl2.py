# ==============================================================================
#     level-1b to level-2 processor
#     This source code is licensed under the 3-clause BSD license found in
#     the LICENSE file in the root directory of this project.
# ==============================================================================

import datetime
import random
import string
import os
import sys
import numpy as np
import pickle as pkl
from tqdm import tqdm
from copy import deepcopy
import netCDF4 as nc
import torch
import yaml

from teds.lib import libNumTools
from teds.lib import libRT
from teds.lib import libATM
from teds.lib import libINV
from teds.lib.libWrite import writevariablefromname
from teds.lib import libRTorCH4_rt as rt

class Emptyclass:
    """Empty class. Data container."""
    
    pass

def get_l1b(l1b_filename):
    # getting l1b data from file
    nc_l1b = nc.Dataset(l1b_filename, mode='r')
        
    l1b_data = {}
    l1b_data['wavelength'] = deepcopy(nc_l1b['observation_data/wavelength'][:])
    l1b_data['radiance'] = deepcopy(nc_l1b['observation_data/radiance'][:])
    l1b_data['noise'] = deepcopy(nc_l1b['observation_data/radiance_stdev'][:])
    l1b_data['sza'] = deepcopy(nc_l1b['geolocation_data/solar_zenith'][:])
    l1b_data['saa'] = deepcopy(nc_l1b['geolocation_data/solar_azimuth'][:])
    l1b_data['vza'] = deepcopy(nc_l1b['geolocation_data/sensor_zenith'][:])
    l1b_data['vaa'] = deepcopy(nc_l1b['geolocation_data/sensor_azimuth'][:])
    l1b_data['latitude'] = deepcopy(nc_l1b['geolocation_data/latitude'][:])
    l1b_data['longitude'] = deepcopy(nc_l1b['geolocation_data/longitude'][:])
    # Extract mask
    nc_var = nc_l1b['observation_data/radiance']
    l1b_data['mask'] = ~nc_var[:].mask
    if not nc_var[:].mask.shape:
        l1b_data['mask'] = np.full(nc_var[:].shape, True)
    
    return (l1b_data)

def get_sgm_atm(filen_sgm_atm):
    data = nc.Dataset(filen_sgm_atm, mode='r')
    atm_sgm = {}
    surf_sgm = {}
    surf_sgm['albedo'] = deepcopy(data['albedo'][:])
    atm_sgm['zlay'] = deepcopy(data['zlay'][:])
    atm_sgm['zlev'] = deepcopy(data['zlev'][:])
    atm_sgm['dcol_co2'] = deepcopy(data['dcol_co2'][:])
    atm_sgm['dcol_ch4'] = deepcopy(data['dcol_ch4'][:])
    atm_sgm['dcol_h2o'] = deepcopy(data['dcol_h2o'][:])
    atm_sgm['col_co2'] = deepcopy(data['XCO2'][:])
    atm_sgm['col_ch4'] = deepcopy(data['XCH4'][:])
    atm_sgm['col_h2o'] = deepcopy(data['XH2O'][:])
    if 'col_air' in data.variables:
        atm_sgm['col_air'] = deepcopy(data['col_air'][:])
    data.close()
    return (surf_sgm, atm_sgm)


def write_gasdata(gas, output_l2, _dims, _dims3d, nalt, nact, nlay, l2product, retrieval_init):
    scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
    # names
    xgas = "X" + gas
    xgas_precision = xgas + " precision"
    xgas_ker = xgas + ' col avg kernel'
    l2_X = np.zeros((nalt, nact))
    l2_X_prec  = np.zeros((nalt, nact))
    l2_X_avgk  = np.zeros((nalt, nact, nlay))
    l2_X_prof  = np.zeros((nalt, nact, nlay))

    for ialt in range(nalt):
        for iact in range(nact):
            tmp = (l2product[ialt, iact][xgas]*scale[gas])
            l2_X[ialt, iact] = tmp
            l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
            l2_X_avgk[ialt, iact, :] = l2product[ialt, iact][xgas_ker][:]
            l2_X_prof[ialt, iact, :] = retrieval_init['trace gases'][gas]['ref_profile']
            
    # write data
    prec_varname   = 'precision'+xgas
    avgker_varname = 'avgkernel'+xgas
    aprof_varname  = 'apriori_profile'+gas

    _ = writevariablefromname(output_l2, xgas,           _dims,   l2_X)
    _ = writevariablefromname(output_l2, prec_varname,   _dims,   l2_X_prec)
    _ = writevariablefromname(output_l2, avgker_varname, _dims3d, l2_X_avgk)
    _ = writevariablefromname(output_l2, aprof_varname,  _dims3d, l2_X_prof)

def write_proxygasdata(gas, output_l2, _dims, nalt, nact, nlay, l2product):
    scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
    # names
    xgas = "X" + gas + " proxy"
    xgas_precision = xgas + " precision"
    l2_X = np.zeros((nalt, nact))
    l2_X_prec = np.zeros((nalt, nact))
    l2_X_accu = np.zeros((nalt, nact))
    l2_X_qav  = np.zeros((nalt, nact), dtype = np.int16)
    
    for ialt in range(nalt):
        for iact in range(nact):
            tmp = (l2product[ialt, iact][xgas]*scale[gas])
            l2_X[ialt, iact] = tmp
            l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
            l2_X_accu[ialt, iact] = 99999.
            l2_X_qav[ialt, iact]  = 100

    # define names in constants_outputvariables
    vargas = "proxyX"+gas
    varprecgas = "precision"+vargas
    varaccugas = "accuracy" +vargas
    varqavgas  = "qa_value" +vargas
    # write data
    _ = writevariablefromname(output_l2, vargas, _dims, l2_X)
    _ = writevariablefromname(output_l2, varprecgas, _dims, l2_X_prec)
    _ = writevariablefromname(output_l2, varaccugas, _dims, l2_X_accu)
    _ = writevariablefromname(output_l2, varqavgas,  _dims, l2_X_qav)

def level2_output(filename, l2product, retrieval_init, l1bproduct, settings):
    # output of level 2 data
    nalt = len(l2product)
    nalt = l1bproduct['latitude'][:, 0].size
    nact = len(l2product[0])
    nlay = len(l2product[0][0]['XCO2 col avg kernel'])

    # filename
    print(filename)
    output_l2 = nc.Dataset(filename, mode='w')
    output_l2.title = 'Tango Carbon E2ES L2 product'
    output_l2.createDimension('number_layers', nlay)     # spectral axis
    output_l2.createDimension('bins_across_track', nact)     # across track axis
    output_l2.createDimension('bins_along_track', nalt)     # along track axis
    output_l2.createDimension('bins_albedo', 1)     # spectral bins albedo
    grp_ns = output_l2.createGroup('non_scattering_retrieval')
    grp_prior = output_l2.createGroup('prior')
    grp_diag = output_l2.createGroup('diagnostics')

    # layer height

    _ = writevariablefromname(output_l2, 'central_layer_height', ('number_layers',),  retrieval_init['zlay'])

    # dimensions
    _dimalt = ('bins_along_track')
    _dims   = ('bins_along_track', 'bins_across_track')
    _dims3d = ('bins_along_track', 'bins_across_track', 'number_layers',)
    # variables
    l2_conv = np.zeros((nalt, nact), dtype = np.int32)
    l2_numb_iter = np.zeros((nalt, nact), dtype = np.int32)
    l2_proc_flag = np.zeros((nalt, nact), dtype = np.int16) + 100
    l2_chi2 = np.zeros((nalt, nact))
    l2_alb = np.zeros((nalt, nact))
    l2_spec_shift  = np.zeros((nalt, nact))
    l2_spec_squeeze = np.zeros((nalt, nact))
    l2_aqutime = np.zeros(nalt) + 99999.   #aquisition time (dummy)
    
    for ialt in range(nalt):
        for iact in range(nact):
            l2_conv[ialt, iact]         = int(l2product[ialt, iact]['convergence'])
            l2_numb_iter[ialt, iact]    = l2product[ialt, iact]['number_iter']
            l2_chi2[ialt, iact]         = l2product[ialt, iact]['chi2']
            l2_alb[ialt, iact]          = l2product[ialt, iact]['alb0']
            l2_spec_shift[ialt, iact]   = l2product[ialt, iact]['spec_shift']
            l2_spec_squeeze[ialt, iact] = l2product[ialt, iact]['spec_squeeze']

    # aquisition time
    _ = writevariablefromname(output_l2, 'aqui_time',       _dimalt, l2_aqutime)
    # processing flag
    _ = writevariablefromname(grp_diag,  'process_flag',    _dims, l2_proc_flag)
    # convergence
    _ = writevariablefromname(grp_diag,  'convergence',     _dims, l2_conv)
    # latitude
    _ = writevariablefromname(output_l2, 'latitude',        _dims, l1bproduct['latitude'][:, :nact])
    # longitude
    _ = writevariablefromname(output_l2, 'longitude',       _dims, l1bproduct['longitude'][:, :nact])
    # maxiterations
    _ = writevariablefromname(grp_diag,  'maxiterations',   _dims, l2_numb_iter)
    # chi2
    _ = writevariablefromname(grp_diag,  'spectralchi2',    _dims, l2_chi2)
    # albedo
    _ = writevariablefromname(output_l2, 'albedo',          _dims, l2_alb)
    # surface pressure
    _ = writevariablefromname(grp_prior, 'surface_pressure',_dims, retrieval_init['surface pressure'])
    # surface elevation
    _ = writevariablefromname(grp_prior, 'surface_elevation',_dims, retrieval_init['surface elevation'])
    # spectral shift
    _ = writevariablefromname(output_l2, 'spectralshift',   _dims, l2_spec_shift)
    # chi2
    _ = writevariablefromname(output_l2, 'spectralsqueeze', _dims, l2_spec_squeeze)
    
    scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
    # names

    gases = ['CO2', 'CH4', 'H2O']
    for gas in gases:
        xgas = "X" + gas
        xgas_precision = xgas + " precision"
        xgas_ker = xgas + ' col avg kernel'

        l2_X = np.zeros((nalt, nact))
        l2_X_prec  = np.zeros((nalt, nact))
        l2_X_avgk  = np.zeros((nalt, nact, nlay))
        l2_X_prof  = np.zeros((nalt, nact, nlay))

        for ialt in range(nalt):
            for iact in range(nact):
                tmp = (l2product[ialt, iact][xgas]*scale[gas])
                l2_X[ialt, iact] = tmp
                l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
                l2_X_avgk[ialt, iact, :] = l2product[ialt, iact][xgas_ker][:]
                l2_X_prof[ialt, iact, :] = retrieval_init['trace gases'][gas]['ref_profile']
            
            # write data
            prec_varname   = 'precision'+xgas
            avgker_varname = 'avgkernel'+xgas
            aprof_varname  = 'apriori_profile'+gas

        _ = writevariablefromname(grp_ns, xgas,           _dims,   l2_X)
        _ = writevariablefromname(grp_ns, prec_varname,   _dims,   l2_X_prec)
        _ = writevariablefromname(output_l2, avgker_varname, _dims3d, l2_X_avgk)
        _ = writevariablefromname(grp_prior, aprof_varname,  _dims3d, l2_X_prof)
    
    #next the main product, proxy
    gases = ['CO2', 'CH4']
    for gas in gases:

        xgas = "X" + gas + " proxy"
        xgas_precision = xgas + " precision"
        l2_X = np.zeros((nalt, nact))
        l2_X_prec = np.zeros((nalt, nact))
        l2_X_accu = np.zeros((nalt, nact))
        l2_X_qav  = np.zeros((nalt, nact), dtype = np.int16)
    
        for ialt in range(nalt):
            for iact in range(nact):
                l2_X[ialt, iact] = (l2product[ialt, iact][xgas]*scale[gas])
                l2_X_prec[ialt, iact] = l2product[ialt, iact][xgas_precision]*scale[gas]
                l2_X_accu[ialt, iact] = 99999.
                l2_X_qav[ialt, iact]  = 100

        # define names in constants_outputvariables
        vargas = "proxyX"+gas
        varprecgas = "precision"+vargas
        varaccugas = "accuracy" +vargas
        varqavgas  = "qa_value" +vargas
        # write data
        _ = writevariablefromname(output_l2, vargas, _dims, l2_X)
        _ = writevariablefromname(output_l2, varprecgas, _dims, l2_X_prec)
        _ = writevariablefromname(output_l2, varaccugas, _dims, l2_X_accu)
        _ = writevariablefromname(output_l2, varqavgas,  _dims, l2_X_qav)

    output_l2.close()


def level2_diags_output(filename, l2product, measurement):
    # this function writes some ectra diagnostics to a file that are not included in the level 2 data product.
    nalt = len(l2product)
    nact = len(l2product[0])
    nwave = len(measurement[0, 0]['wavelength'])
    # file name
    print('diag_output')
    output_l2diag = nc.Dataset(filename, mode='w')
    output_l2diag.title = 'Tango Carbon E2ES L2 diagnostics'
    output_l2diag.createDimension('bins_spectral', nwave)     # spectral axis
    output_l2diag.createDimension('bins_along_track', nalt)     # across track axis
    output_l2diag.createDimension('bins_across_track', nact)     # across track axis

    # create 2d arrays
    l2diag_wave = np.zeros((nalt, nact, nwave))
    l2diag_measurement = np.zeros((nalt, nact, nwave))
    l2diag_sun_irradiance = np.zeros((nalt, nact, nwave))
    l2diag_gainCO2 = np.zeros((nalt, nact, nwave))
    l2diag_gainH2O = np.zeros((nalt, nact, nwave))
    l2diag_gainCH4 = np.zeros((nalt, nact, nwave))
    for ialt in range(nalt):
        for iact in range(nact):
            l2diag_wave[ialt, iact, :] = measurement[ialt, iact]['wavelength'][:]
            l2diag_measurement[ialt, iact, :] = measurement[ialt, iact]['ymeas'][:]
            l2diag_sun_irradiance[ialt, iact, :] = measurement[ialt, iact]['sun'][:]
            l2diag_gainCO2[ialt, iact, :] = l2product[ialt, iact]['gain_XCO2'][:]
            l2diag_gainH2O[ialt, iact, :] = l2product[ialt, iact]['gain_XH2O'][:]
            l2diag_gainCH4[ialt, iact, :] = l2product[ialt, iact]['gain_XCH4'][:]

    _dims = ('bins_along_track', 'bins_across_track', 'bins_spectral')
    # wavelength of measurement
    _ = writevariablefromname(output_l2diag, "l2wavelength", _dims, l2diag_wave)
    # spectral radiance measurements
    _ = writevariablefromname(output_l2diag, "l2measurement", _dims, l2diag_measurement)
    # solar irradiance spectrum
    _ = writevariablefromname(output_l2diag, "l2solarirradiance", _dims, l2diag_sun_irradiance)
    # gain CO2
    _ = writevariablefromname(output_l2diag, "l2gainCO2", _dims, l2diag_gainCO2)
    # gain CH4
    _ = writevariablefromname(output_l2diag, "l2gainCH4", _dims, l2diag_gainCH4)
    # gain H2O
    _ = writevariablefromname(output_l2diag, "l2gainH2O", _dims, l2diag_gainH2O)
    output_l2diag.close()

def level2_output_RTorCH4(
    config,
    l2,
    species_names,
    nalt,
    nact,
    nlay,
    extra_dimensions = [],
    extra_variables = []
):
    """Write level2 output from RTorCH4 retrieval.

    In addition to the dictionary output by
    level1b_to_level2_processor_RTorCH4, "unknown" additional variables may be
    inserted in the output using the extra_* parameters:

    extra_dimensions is a list of (dimension_name, dimension_size)-tuples for
    extra NetCDF dimensions.

    extra_variables is a list of
    (variable_name, parent_group, dimension_list, data, [type])-tuples of the
    extra variables to output. The variable will be placed in the group named
    by parent_group; this must be a nonempty string, and not one of the groups
    created by default (diagnostics/non_scattering_retrieval/prior). All
    unique parent_groups encountered in extra_variables will be created.
    dimension_list is a list of dimension names. type is the NetCDF variable
    type; if it is omitted, 'f8' is assumed.
    """
    # We try to use config['io_files']['output_l2'] as the output file.
    # However, if it already exists or otherwise cannot be opened for writing
    # (e.g. because it's being locked by another program), we will try a few
    # other paths to prevent a long retrieval being lost at the very end.
    outpath = config['io_files']['output_l2']
    outdir = os.path.dirname(outpath)
    outfile = os.path.basename(outpath)
    randoutfile = ''.join([
        random.choice(string.ascii_letters+string.digits)
        for _ in range(8)
    ]) + '_' + outfile

    try_paths = [
        outpath,
        os.path.join(outdir, randoutfile),
        os.path.join('.', randoutfile),
        os.path.join('/tmp', randoutfile)
    ]

    out = None
    for p in try_paths:
        if os.path.exists(p):
            continue

        try:
            out = nc.Dataset(p, mode='w')
            print('Writing level 2 output to', p)
            break
        except:
            continue

    out.title = 'Tango Carbon E2ES L2 product'
    out.createDimension('number_layers', nlay)
    out.createDimension('bins_across_track', nact)
    out.createDimension('bins_along_track', nalt)
    out.createDimension('bins_albedo', 1)
    grp_ns = out.createGroup('non_scattering_retrieval')
    grp_prior = out.createGroup('prior')
    grp_diag = out.createGroup('diagnostics')
    dimalt = ('bins_along_track',)
    dims2d = ('bins_along_track', 'bins_across_track')
    dims3d = ('bins_along_track', 'bins_across_track', 'number_layers')
    dimlay = ('number_layers',)

    # Keys for variables that go in subgroups
    keys_ns = [
        f'X{spn}' for spn in species_names
    ] + [
        f'precisionX{spn}' for spn in species_names
    ]

    keys_prior = [
        f'apriori_profile{spn}' for spn in species_names
    ] + ['surface_elevation', 'surface_pressure']

    keys_diag = [
        'spectralchi2',
        'convergence',
        'maxiterations',
        'process_flag'
    ]

    # Keys for variables with dimensions other than (alt, act)
    keys_alt = [
        'aqui_time'
    ]

    keys_3d = [
        f'avgkernelX{spn}' for spn in species_names
    ] + [
        f'apriori_profile{spn}' for spn in species_names
    ]

    keys_lay = [
        'central_layer_height'
    ]

    # Write all variables
    for k, v in l2.items():
        dims = dims2d
        if k in keys_alt:
            dims = dimalt
        elif k in keys_3d:
            dims = dims3d
        elif k in keys_lay:
            dims = dimlay

        try:
            if k in keys_ns:
                writevariablefromname(grp_ns, k, dims, v)
            elif k in keys_prior:
                writevariablefromname(grp_prior, k, dims, v)
            elif k in keys_diag:
                writevariablefromname(grp_diag, k, dims, v)
            else:
                writevariablefromname(out, k, dims, v)
        except Exception as e:
            print(f'Error writing variable {k}', file=sys.stderr)
            raise e

    for dimname, dimsize in extra_dimensions:
        out.createDimension(dimname, dimsize)

    extra_groups = {}
    for vartup in extra_variables:
        varname, pargroup, dimslist, vardata, vartype = (
            None, None, None, None, 'f8'
        )
        if len(vartup) == 4:
            varname, pargroup, dimlist, vardata = vartup
        elif len(vartup) == 5:
            varname, pargroup, dimlist, vardata, vartype = vartup

        if pargroup not in extra_groups:
            extra_groups[pargroup] = out.createGroup(pargroup)

        v = extra_groups[pargroup].createVariable(varname, vartype, dimlist)
        v[:] = vardata

    out.close()

def level1b_to_level2_processor_RTorCH4(config):
    """Compute level 2 data product using RTorCH4."""

    print("RTorCH4 level 1B to level 2 processor")

    l1b = get_l1b(config['io_files']['input_l1b'])

    # Hardcode the number of albedo coefficients for now
    N_alb = 2

    # Number of iterations
    N_iter = config['retrieval_init']['max_iter']

    # get sgm geo data
    surf_sgm, atm_sgm = get_sgm_atm(config['io_files']['input_sgm'])

    # mask on radiance nan values
    mask = l1b['mask']

    # Internal lbl spectral grid
    wave_start = config['spec_settings']['wavestart']
    wave_end = config['spec_settings']['waveend']
    wave_extend = config['spec_settings']['wave_extend']
    dwave_lbl = config['spec_settings']['dwave']
    wave_lbl = np.arange(
        wave_start-wave_extend,
        wave_end+wave_extend,
        dwave_lbl
    )  # nm
    N_lbl_wls = len(wave_lbl)

    # Geometry dimensions: along_track, across_track
    sza_full = l1b['sza']
    vza_full = l1b['vza']
    N_alt, N_act = sza_full.shape

    # Vertical layering
    N_layers = config['atmosphere']['nlay']
    dzlay = config['atmosphere']['dzlay']
    psurf = config['atmosphere']['psurf']
    nlev = N_layers + 1  # number of levels

    # altitude of layer midpoint
    zlay = (np.arange(N_layers-1, -1, -1)+0.5)*dzlay
    # altitude of layer interfaces = levels
    zlev = np.arange(nlev-1, -1, -1)*dzlay

    # model atmosphere
    atm = libATM.atmosphere_data(zlay, zlev, psurf)
    atm.get_data_AFGL(config['io_files']['input_afgl'])

    sun = libRT.read_sun_spectrum_TSIS1HSRS(
        config['io_files']['input_sun_reference']
    )
    sun_lbl = libRT.interpolate_sun(sun, wave_lbl)
    del sun

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)] # see hapi manual sec 6.6
    molec = libRT.molecular_data(wave_lbl)
    molec.get_data_HITRAN(config['io_files']['input_hapi'], iso_ids)

    # Calculate optical properties
    # Init class with optics.prop dictionary
    optics = libRT.optic_abs_prop(wave_lbl, zlay)
    # If pickle file exists read from file
    if (
        (not os.path.exists(config['io_files']['dump_xsec']))
        or config['xsec_forced']
    ):
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm)
        # Dump optics.prop dictionary into temporary pkl file
        pkl.dump(optics.prop, open(config['io_files']['dump_xsec'], 'wb'))
    else:
        # Read optics.prop dictionary from pickle file
        optics.prop = pkl.load(open(config['io_files']['dump_xsec'], 'rb'))

    # Molecular species (by HITRAN id) to retrieve.
    # The order of these species is fixed from this point onwards.
    species = ['molec_01', 'molec_32', 'molec_07']
    species_names = ['H2O', 'CH4', 'CO2'] # order must match with species!
    N_species = len(species)

    # Create a matrix from the dictionary.
    # Conversion factor: cm^-2 -> m^-2
    xsec = np.stack(
        [optics.prop[sp]['xsec'].T for sp in species]
    )*1e-4

    del molec, optics

    # First guess for trace gas concentrations. N.B. same order as species!
    init_conc = np.array([1e-3, 1.7e-6, 3e-4])
    print('Trace gas concentration priors:')
    for spi, spn in enumerate(species_names):
        key = f'prior_X{spn}'
        if key in config['retrieval_init']:
            init_conc[spi] = config['retrieval_init'][key]
        print(f'{spn}: {init_conc[spi]:e}')

    # Column profile matrix; each pixel can have different profiles if using
    # the SGM as the data source. If AFGL is used, promote the per-granule
    # profile to a per-pixel one.
    col_profiles = np.zeros((N_alt, N_act, N_species, N_layers))

    if config['retrieval_init']['sw_prof_init'] not in ('afgl', 'sgm'):
        raise RuntimeError(
            'sw_prof_init of l1bl2 configuration not set correctly!'
        )

    airsum = None
    for spi, spn in enumerate(species_names):
        if config['retrieval_init']['sw_prof_init'] == 'afgl':
            col_profiles[:,:,spi,:] = np.einsum(
                "xyz,z->xyz",
                np.ones((N_alt, N_act, N_layers)),
                getattr(atm, spn)
            )
            airsum = np.sum(atm.air)*np.ones((N_alt, N_act))
        elif config['retrieval_init']['sw_prof_init'] == 'sgm':
            col_profiles[:,:,spi,:] = (
                atm_sgm[f'dcol_{spn.lower()}'][:N_alt,:N_act,:]
            )
            airsum = atm_sgm['col_air'][:N_alt,:N_act]

    del atm

    init_col_scales_full = np.einsum(
        "xys,xy,s->xys",
        1/np.sum(col_profiles, axis=-1),
        airsum,
        init_conc
    )

    # Wavelength dimensions: across_track, wavelength
    obs_wls = l1b['wavelength']
    N_obs_wls = obs_wls.shape[1]
    obs_fwhm = (
        config['isrf_settings']['fwhm'] * np.ones(N_obs_wls)
    )

    isrf_beta = 2
    if config['isrf_settings']['type'] == 'generalized_normal':
        isrf_beta = 1/config['isrf_settings']['bcoeff']

    # Diagonal of all inverse covariance matrices
    # Shape: along_track, across_track, wavelength
    invcov_diag = 1/np.clip(l1b['noise']**2, 1e-9, None)

    # Ignore pixels where one or more wavelength components are masked
    ignore_pixels_full = np.any(~l1b['mask'], axis=2)

    # If there is a convergence_criteria block in the config file, use them;
    # otherwise, use backwards-compatible chi2_lim.
    #
    # A batch is considered converged if:
    # (
    #   prev_chi2-chi2 < convergence_deltaplus for all pixels AND
    #   prev_chi2-chi2 > convergence_deltaminus for all pixels AND
    #   1-convergence_epsilon < chi2/prev_chi2 < 1+convergence_epsilon
    #       for all pixels AND
    #   1-convergence_mean_epsilon < ⟨chi2/prev_chi2⟩
    #       < 1+congergence_mean_epsilon
    # ) OR chi2 < convergence_chi2 for all pixels
    convergence_deltaplus = 0.05
    convergence_deltaminus = -0.05
    convergence_epsilon = 1e8
    convergence_mean_epsilon = 1e8
    convergence_chi2 = 1e-2
    if 'chi2_lim' in config['retrieval_init']:
        convergence_deltaplus = config['retrieval_init']['chi2_lim']
        convergence_deltaminus = -config['retrieval_init']['chi2_lim']

    if 'convergence_criteria' in config:
        cc = config['convergence_criteria']
        if 'deltaplus' in cc:
            convergence_deltaplus = cc['deltaplus']
        if 'deltaminus' in cc:
            convergence_deltaminus = cc['deltaminus']
        if 'epsilon' in cc:
            convergence_epsilon = cc['epsilon']
        if 'mean_epsilon' in cc:
            convergence_mean_epsilon = cc['mean_epsilon']
        if 'chi2' in cc:
            convergence_chi2 = cc['chi2']

    # Output variables
    out_col_scales = np.nan * np.ones((N_alt, N_act, N_species))
    out_alb = np.nan * np.ones((N_alt, N_act, N_alb))
    out_chi2 = np.nan * np.ones((N_alt, N_act, N_iter+1,))
    out_Ainv = None
    out_B = None
    out_col_kern_base = np.nan * np.ones((N_alt, N_act, N_species, N_layers))
    iterations = np.zeros((N_alt, N_act), dtype=int)
    converged = np.zeros((N_alt, N_act), dtype=bool)

    # Batch size can be set in config. If not given, a while ACT line will be
    # used.
    batch_size = N_act
    if 'batch_size' in config['retrieval_init']:
        batch_size = min(N_act, config['retrieval_init']['batch_size'])
    act_inds = np.arange(N_act)
    act_chunks = np.split(act_inds, act_inds[::batch_size])[1:]
    N_chunks = len(act_chunks)

    print(f'Using {N_chunks} ACT chunk'+('s' if N_chunks != 1 else ''))

    # dealloc_chunk: if True, we do not keep a radtran object for each ACT
    # chunk, instead creating a new object every time it is needed. This saves
    # memory at the cost of greater processing overhead.
    deallocate_chunk = (
        'deallocate_chunk' in config['retrieval_init']
        and config['retrieval_init']['deallocate_chunk']
    )

    # We use Gauss-Newton, so we don't need autograd at all -> no_grad context
    with torch.no_grad():#, torch.profiler.profile(
#            schedule=torch.profiler.schedule(
#                wait=1, warmup=1, active=10, repeat=1
#            ),
#            on_trace_ready=torch.profiler.tensorboard_trace_handler(
#                "/Users/petersk/tmp/tango_retrieval_prof"
#            ),
#            activities=[torch.profiler.ProfilerActivity.CPU],
#            record_shapes=True,
#            with_stack=False,#True,
#            profile_memory=True
#    ) as prof:
        dtype = torch.float32
        if (
            'use_float64' in config['retrieval_init']
            and config['retrieval_init']['use_float64']
        ):
            dtype = torch.float64
        cpu = torch.device('cpu')
        device = torch.device('cuda') if torch.cuda.is_available() else cpu
        T = lambda a: torch.tensor(a, device=device, dtype=dtype)

        obs_wls = T(obs_wls)
        obs_fwhm = T(obs_fwhm)
        wave_lbl = T(wave_lbl)
        xsec = T(xsec)
        col_profiles = T(col_profiles)
        init_col_scales_full = T(init_col_scales_full)
        sun_lbl = T(sun_lbl)
        sza_full = T(sza_full)
        vza_full = T(vza_full)
        radiance_full = T(l1b['radiance'])
        invcov_diag = T(invcov_diag)
        ignore_pixels_full = torch.tensor(
            ignore_pixels_full,
            dtype=bool, device=device
        )

        # The forward model and irradiance at observed wavelengths are
        # computed at the first batch.
        radtran = [None]*N_chunks
        isrfs = [rt.GeneralisedNormalISRF(isrf_beta) for _ in range(N_chunks)]
        sun_obs = [None]*N_chunks

        retstart = datetime.datetime.now()
        print('Begin retrieval at', retstart)

        # ALT scanline main loop
        for batch in tqdm(range(N_alt)):
            # If we have a profiler, step it; if not, be quiet about it.
            try:
                prof.step()
            except NameError:
                pass

            # ACT chunk main loop
            for chunk, chunk_inds in enumerate(act_chunks):
                sza = sza_full[batch,chunk_inds]
                vza = vza_full[batch,chunk_inds]
                radiance = radiance_full[batch,chunk_inds,:]
                valid_mask = ~ignore_pixels_full[batch,chunk_inds]

                N = len(chunk_inds)

                tau_base = torch.einsum(
                    "szl,Nsz->Nszl",
                    xsec,
                    col_profiles[batch,chunk_inds]
                )

                if radtran[chunk] is None:
                    radtran[chunk] = rt.RTForwardModel(
                        obs_wls[chunk_inds,:],
                        obs_fwhm,
                        wave_lbl,
                        tau_base,
                        sun_lbl,
                        sza,
                        vza,
                        tau_offset=None,
                        isrf=isrfs[chunk],
                        dtype=dtype,
                        device=device
                    )
                else:
                    # In the forward model, only the geometry/optical depth
                    # changes per batch, so use setup_measurement to reduce
                    # init time.
                    radtran[chunk].setup_measurement(
                        sza,
                        vza,
                        tau_base=tau_base
                    )

                if sun_obs[chunk] is None:
                    isrfs[chunk].set_parameters(
                        obs_wls[chunk_inds,:],
                        obs_fwhm,
                        wave_lbl,
                        cache_hint=True
                    )
                    sun_obs[chunk] = isrfs[chunk].convolve(
                        torch.outer(
                            torch.ones(N, dtype=dtype, device=device),
                            sun_lbl
                        )
                    )

                idx = torch.argmax(radiance, dim=1)
                pixidx = torch.arange(N)
                alb_first_guess = (
                    torch.pi * radiance[pixidx,idx]
                    / sun_obs[chunk][pixidx,idx]
                    / torch.cos(sza*torch.pi/180.)
                )
                # Use first guess for baseline albedo; set higher-order
                # coefficients to 0.
                alb = torch.stack(
                    (
                        (alb_first_guess,)
                        + (torch.zeros_like(alb_first_guess),) * (N_alb-1)
                    ),
                    dim=-1
                )
                alb_prior = torch.detach_copy(alb)

                col_scales = torch.clone(
                    init_col_scales_full[batch,chunk_inds,:]
                )

                waveshift = None

                # Inverse covariance
                invcov = torch.einsum(
                    "NO,Oo->NOo",
                    invcov_diag[batch,chunk_inds,:],
                    torch.eye(N_obs_wls, dtype=dtype, device=device)
                )

                chi2 = None
                prev_chi2 = None
                Ainv = None
                B = None

                # Gauss-Newton solver
                retrieval = rt.GaussNewtonRetrieval(
                    radtran[chunk],
                    radiance,
                    inverse_covariance = invcov
                )

                # Gauss-Newton fit main loop
                for i in range(N_iter):
                    i1 = i+1

                    new_col_scales, new_alb, new_ws, chi2, Ainv, B = (
                        retrieval.inversion_step(
                            col_scales,
                            alb,
                            waveshift,
                            ignore_pixels=~valid_mask,
                            aux_matrices=True
                        )
                    )
                    # If the step introduces NaNs, infinities, or column
                    # scales below 0, it is mathematically invalid and must be
                    # rejected outright.  We further impose some constraints
                    # on the regime of physically reasonable parameters.
                    valid_mask &= (
                        torch.all(torch.isfinite(new_col_scales), dim=1)
                        & torch.all(torch.isfinite(new_alb), dim=1)
                        & torch.all(new_col_scales >= 0, dim=1)
                        & torch.isfinite(chi2)
                    )

                    chi2 = torch.clamp(
                        chi2.detach(),
                        min=radtran[chunk].epsilon
                    )
                    chi2[~valid_mask] = torch.nan

                    col_scales[valid_mask] = new_col_scales[valid_mask]
                    col_scales[~valid_mask] = torch.nan
                    alb[valid_mask] = new_alb[valid_mask]
                    alb[~valid_mask] = torch.nan

                    out_chi2[batch,chunk_inds,i] = chi2.cpu().detach().numpy()
                    iterations[batch,chunk_inds] += torch.where(
                        valid_mask, 1, 0
                    ).cpu().numpy()

                    # Test for convergence
                    if prev_chi2 is not None:
                        rat_chi = chi2[valid_mask]/prev_chi2[valid_mask]
                        mrat_chi = rat_chi.mean()
                        delta_chi = prev_chi2[valid_mask]-chi2[valid_mask]
                        if (
                            (
                                torch.all(delta_chi < convergence_deltaplus)
                                and torch.all(
                                    delta_chi > convergence_deltaminus
                                )
                                and torch.all(rat_chi > 1-convergence_epsilon)
                                and torch.all(rat_chi < 1+convergence_epsilon)
                                and (mrat_chi > (1-convergence_mean_epsilon))
                                and (mrat_chi < (1+convergence_mean_epsilon))
                            )
                            or torch.all(chi2 < convergence_chi2)
                        ):
                            converged[batch,chunk_inds] = (
                                valid_mask.cpu().detach().numpy()
                            )
                            break
                    prev_chi2 = chi2.detach()

                # We've finished optimisation here,
                # but still want the final χ²
                y, J = radtran[chunk].forward_inst(
                    col_scales,
                    alb,
                    waveshift,
                    grad=True
                )
                dy = y-radiance
                chi2 = torch.einsum(
                    "NO,NOo,No->N", dy, retrieval.inv_cov, dy
                ) / (N_obs_wls-radtran[chunk].N_params)

                if torch.any(torch.isnan(chi2)) and dtype != torch.float64:
                    dy64 = dy.to(torch.float64)
                    chi2 = (
                        torch.einsum(
                            "NO,NOo,No->N", dy64, retrieval.inv_cov64, dy64
                        ) / (N_obs_wls-radtran[chunk].N_params)
                    ).to(dtype)

                for i in range(i1, N_iter+1):
                    out_chi2[batch,chunk_inds,i] = chi2.cpu().detach().numpy()

                out_col_scales[batch,chunk_inds,:] = (
                    col_scales.cpu().detach().numpy()
                )
                out_alb[batch,chunk_inds,:] = alb.cpu().detach().numpy()

                if out_Ainv is None:
                    out_Ainv = np.nan * np.ones(
                        (N_alt, N_act) + tuple(Ainv.shape[1:])
                    )
                out_Ainv[batch,chunk_inds] = Ainv.cpu().detach().numpy()

                if out_B is None:
                    out_B = np.nan * np.ones(
                        (N_alt, N_act) + tuple(B.shape[1:])
                    )
                out_B[batch,chunk_inds] = B.cpu().detach().numpy()

                tau_per_molec, tau_layers = (
                    radtran[chunk].last_output__compute_optical_depth
                )
                rad_lbl, dev_tau_lbl, dev_alb_lbl = (
                    radtran[chunk].rt_solver.last_output__radiance_toa
                )
                out_col_kern_base[batch,chunk_inds,:,:] = torch.einsum(
                    "Nso,Nol,Nszl,Nl->Nsz",
                    B[:,:N_species,:],
                    isrfs[chunk].tensor(),
                    tau_per_molec,
                    dev_tau_lbl
                ).cpu().detach().numpy()

                if deallocate_chunk:
                    radtran[chunk] = None

        max_iter_done = np.amax(iterations[~np.isnan(iterations)])
        out_chi2 = out_chi2[:,:,:max_iter_done+1]
        retend = datetime.datetime.now()
        print('Retrieval completed in', retend-retstart)

        l2 = {}
        l2['albedo'] = out_alb[:,:,0]
        l2['spectralchi2'] = out_chi2[:,:,-1]
        l2['convergence'] = converged.astype(int)
        l2['maxiterations'] = iterations
        l2['aqui_time'] = np.zeros(N_alt) + 99999.
        l2['process_flag'] = (
            np.zeros((N_alt, N_act), dtype=np.int16) + 100
        )
        l2['latitude'] = l1b['latitude']
        l2['longitude'] = l1b['longitude']
        l2['surface_pressure'] = np.zeros((N_alt, N_act)) + 1013.
        l2['surface_elevation'] = np.zeros((N_alt, N_act))
        l2['spectralshift'] = np.zeros((N_alt, N_act))
        l2['spectralsqueeze'] = np.zeros((N_alt, N_act))
        l2['central_layer_height'] = zlay

        col_profiles = col_profiles.cpu().detach().numpy()

        for spi, spn in enumerate(species_names):
            spx = species[spi]

            cscale = out_col_scales[:,:,spi]
            nanmask = np.isnan(cscale)
            colpri = col_profiles[:,:,spi,:]
            colsum = np.sum(colpri, axis=-1)
            l2[f'apriori_profile{spn}'] = colpri

            X = cscale * colsum/airsum
            prec = X * np.sqrt(out_Ainv[:,:,spi,spi])
            prec[nanmask] = np.nan
            l2[f'X{spn}'] = X
            l2[f'X{spn}'][nanmask] = np.nan
            l2[f'precisionX{spn}'] = prec

            delta_prof = np.einsum(
                "lc,lcz->lcz",
                colsum,
                1/colpri
            )
            l2[f'avgkernelX{spn}'] = out_col_kern_base[:,:,spi,:]*delta_prof
            # N.B. column averaging kernels generally differ from Jochen's
            # version due to different priors.

        scaling = {
            'CO2': 1E-6,
            'CH4': 1E-9,
            'H2O': 1E-6
        }
        # Proxy model
        xch4_model = 1.8E-6
        xco2_model = 405.E-6
        l2['proxyXCO2'] = l2['XCO2']/l2['XCH4']*xch4_model/scaling['CO2']
        l2['proxyXCH4'] = l2['XCH4']/l2['XCO2']*xco2_model/scaling['CH4']
        rel_error = np.sqrt(
            (l2['precisionXCO2']/l2['XCO2'])**2
            + (l2['precisionXCH4']/l2['XCH4'])**2
        )
        l2['precisionproxyXCO2'] = rel_error*l2['proxyXCO2']
        l2['precisionproxyXCH4'] = rel_error*l2['proxyXCH4']

        l2['accuracyproxyXCO2'] = 99999. * np.ones_like(l2['proxyXCO2'])
        l2['qa_valueproxyXCO2'] = 100 * np.ones(
            l2['proxyXCO2'].shape, dtype=np.int16
        )

        l2['accuracyproxyXCH4'] = 99999. * np.ones_like(l2['proxyXCH4'])
        l2['qa_valueproxyXCH4'] = 100 * np.ones(
            l2['proxyXCH4'].shape, dtype=np.int16
        )

        for spn in scaling.keys():
            l2[f'X{spn}'] /= scaling[spn]
            l2[f'precisionX{spn}'] /= scaling[spn]

        level2_output_RTorCH4(
            config, l2, species_names, N_alt, N_act, N_layers
        )



def level1b_to_level2_processor(config, sw_diag_output = False):
    # get the l1b files

    print('level 1B to 2 proessor ...')
    l1b = get_l1b(config['io_files']['input_l1b'])
    
    # get sgm geo data
    surf_sgm, atm_sgm = get_sgm_atm(config['io_files']['input_sgm'])

    # mask on radiance nan values
    mask = l1b['mask']

    # Internal lbl spectral grid
    wave_start = config['spec_settings']['wavestart']
    wave_end = config['spec_settings']['waveend']
    wave_extend = config['spec_settings']['wave_extend']
    dwave_lbl = config['spec_settings']['dwave']
    wave_lbl = np.arange(wave_start-wave_extend, wave_end+wave_extend, dwave_lbl)  # nm

    # Vertical layering
    nlay = config['atmosphere']['nlay']
    dzlay = config['atmosphere']['dzlay']
    psurf = config['atmosphere']['psurf']
    nlev = nlay + 1  # number of levels

    zlay = (np.arange(nlay-1, -1, -1)+0.5)*dzlay  # altitude of layer midpoint
    zlev = np.arange(nlev-1, -1, -1)*dzlay  # altitude of layer interfaces = levels

    # model atmosphere

    atm = libATM.atmosphere_data(zlay, zlev, psurf)
    atm.get_data_AFGL(config['io_files']['input_afgl'])

    # XCO2 = np.sum(atm.CO2)/np.sum(atm.air)
    # XCH4 = np.sum(atm.CH4)/np.sum(atm.air)
    # XH2O = np.sum(atm.H2O)/np.sum(atm.air)
    # print(XCO2*1.E6, XCH4*1.E9, XH2O*1E6)

    sun = libRT.read_sun_spectrum_TSIS1HSRS(config['io_files']['input_sun_reference'])
    sun_lbl = libRT.interpolate_sun(sun, wave_lbl)

    # Download molecular absorption parameter
    iso_ids = [('CH4', 32), ('H2O', 1), ('CO2', 7)]  # see hapi manual  sec 6.6
    molec = libRT.molecular_data(wave_lbl)
    molec.get_data_HITRAN(config['io_files']['input_hapi'], iso_ids)

    # Calculate optical properties
    # If pickle file exists read from file
    if ((not os.path.exists(config['io_files']['dump_xsec'])) or config['xsec_forced']):
        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)
        # Molecular absorption optical properties
        optics.cal_molec_xsec(molec, atm)
        # Dump optics.prop dictionary into temporary pkl file
        pkl.dump(optics.prop, open(config['io_files']['dump_xsec'], 'wb'))
    else:
        # Init class with optics.prop dictionary
        optics = libRT.optic_abs_prop(wave_lbl, zlay)
        # Read optics.prop dictionary from pickle file
        optics.prop = pkl.load(open(config['io_files']['dump_xsec'], 'rb'))

    optics.set_opt_depth_species(atm, ['molec_01', 'molec_32', 'molec_07'])

    # Initialization of the least squares fit
    retrieval_init = {}
    retrieval_init['chi2 limit'] = config['retrieval_init']['chi2_lim']
    retrieval_init['maximum iteration'] = config['retrieval_init']['max_iter']
    retrieval_init['zlay'] = zlay
    retrieval_init['trace gases'] = {'CO2': {'init': 300,      'scaling': 1E-6},
                                     'CH4': {'init': 1700,     'scaling': 1E-9},
                                     'H2O': {'init': 1000,     'scaling': 1E-6}}

#    retrieval_init['surface'] = {'alb0': 0.17}
    retrieval_init['wavelength lbl'] = wave_lbl
    retrieval_init['solar irradiance'] = sun_lbl

    nact = l1b['latitude'][0, :].size
    nalt = l1b['latitude'][:, 0].size

    l2product = {}
    xch4_model = 1.8E-6
    xco2_model = 405.E-6

    XCO2 = np.zeros(nact)
    XCO2_true_smoothed = np.zeros(nact)
    XCO2_true = np.zeros(nact)
    XCO2_prec = np.zeros(nact)

    l2product = np.ndarray((nalt, nact), np.object_)
    measurement = np.ndarray((nalt, nact), np.object_)

#    plt.plot(l1b['wavelength'][0,mask[0,0,:]], l1b['radiance'][0,0,mask[0,0,:]])
#    plt.plot(l1b['wavelength'][10,mask[0,9,:]], l1b['radiance'][0,9,mask[0,9,:]])
#    sys.exit()
    # We introduce two types of ialt indices, ialt points to l1b data structure and ilat_l2 to l2 data structure.
    # Later is different to l1b structure because of option for image selection (sw_ALT_select).

    runtime_cum = {}
    runtime_cum['opt'] = 0.
    runtime_cum['rtm'] = 0.
    runtime_cum['conv'] = 0.
    runtime_cum['kern'] = 0.

    for ialt, ialt in enumerate(tqdm(range(nalt))):
        for iact in range(nact):
            # number of 'good' pixels
            numb_spec_points = np.sum(mask[ialt,iact,:])/float(mask[ialt,iact,:].size)
            if(numb_spec_points >=0.9):
                # initialization of pixel retrieval
                xco2_ref = 405.  # ppm
                xco2 = np.sum(atm.CO2) / np.sum(atm.air) * 1.E6
    
                if (config['retrieval_init']['sw_prof_init'] == 'afgl'):
                    retrieval_init['trace gases']['CO2']['ref_profile'] = atm.CO2*xco2_ref/xco2
                    retrieval_init['trace gases']['CH4']['ref_profile'] = atm.CH4
                    retrieval_init['trace gases']['H2O']['ref_profile'] = atm.H2O
                elif (config['retrieval_init']['sw_prof_init'] == 'sgm'):
                    retrieval_init['trace gases']['CO2']['ref_profile'] = atm_sgm['dcol_co2'][ialt, iact, :]*xco2_ref/xco2
                    retrieval_init['trace gases']['CH4']['ref_profile'] = atm_sgm['dcol_ch4'][ialt, iact, :]
                    retrieval_init['trace gases']['H2O']['ref_profile'] = atm_sgm['dcol_h2o'][ialt, iact, :]
                else:
                    sys.exit('sw_prof_init of l1bl2 configuration not set correctly!')
                    
                retrieval_init['surface pressure']  = np.zeros([nalt,nact])+1013.  #these are dummy values for the time being
                retrieval_init['surface elevation'] = np.zeros([nalt,nact])

                wavelength = l1b['wavelength'][iact, mask[ialt,iact, :]].data
                istart = np.argmin(np.abs(wavelength - wave_start))
                iend = np.argmin(np.abs(wavelength - wave_end))
                wave_meas = wavelength[istart:iend+1]  # nm
                
                # define isrf function
    
                isrf_convolution = libNumTools.get_isrf(wave_meas, wave_lbl, config['isrf_settings'])
    
                atm_ret = deepcopy(atm)  # to initialize each retrieval with the same atmosphere
    
                sun = isrf_convolution(sun_lbl)
    
                # Observation geometry
                nwave = wave_meas.size
                measurement[ialt, iact] = {}
                measurement[ialt, iact]['wavelength'] = wave_meas[:]
                measurement[ialt, iact]['mu0'] = np.cos(np.deg2rad(l1b['sza'][ialt, iact]))
                measurement[ialt, iact]['muv'] = np.cos(np.deg2rad(l1b['vza'][ialt, iact]))
                measurement[ialt, iact]['ymeas'] = l1b['radiance'][ialt, iact, mask[ialt,iact, :]][istart:iend+1].data
                measurement[ialt, iact]['Smeas'] = np.eye(
                    nwave)*(l1b['noise'][ialt, iact, mask[ialt,iact, :]][istart:iend+1].data)**2
                measurement[ialt, iact]['sun'] = sun
    
                # derive first guess albedo from the maximum reflectance
                ymeas_max = np.max(measurement[ialt, iact]['ymeas'])

                idx = np.where(measurement[ialt, iact]['ymeas'] == ymeas_max)[0][0]
                alb_first_guess = measurement[ialt, iact]['ymeas'][idx] / \
                    sun[idx]*np.pi/np.cos(np.deg2rad(l1b['sza'][ialt, iact]))
                retrieval_init['surface'] = {'alb0': alb_first_guess, 'alb1': 0.0}
    
                # Non-scattering least squares fit
                l2product[ialt, iact], runtime = libINV.Gauss_Newton_iteration(
                    retrieval_init, atm_ret, optics, measurement[ialt, iact],
                    config['retrieval_init']['max_iter'], config['retrieval_init']['chi2_lim'],
                    isrf_convolution)
    
                for key in runtime.keys():
                    runtime_cum[key] = runtime_cum[key] + runtime[key]
    
                if (not l2product[ialt, iact]['convergence']):
                    print('pixel did not converge (ialt,iact) = ', ialt, iact)
    
                # Define proxy product
                l2product[ialt, iact]['XCO2 proxy'] = l2product[ialt,
                                                                   iact]['XCO2']/l2product[ialt, iact]['XCH4']*xch4_model
                l2product[ialt, iact]['XCH4 proxy'] = l2product[ialt,
                                                                   iact]['XCH4']/l2product[ialt, iact]['XCO2']*xco2_model
                rel_error = np.sqrt((l2product[ialt, iact]['XCO2 precision']/l2product[ialt, iact]['XCO2'])**2 +
                                    (l2product[ialt, iact]['XCH4 precision']/l2product[ialt, iact]['XCH4'])**2)
                l2product[ialt, iact]['XCO2 proxy precision'] = rel_error * l2product[ialt, iact]['XCO2']
                l2product[ialt, iact]['XCH4 proxy precision'] = rel_error * l2product[ialt, iact]['XCH4']

            else:
                l2product[ialt, iact] = level2_nan(retrieval_init, nlay)
                l2product[ialt, iact]['XCO2 proxy'] = float("nan")
                l2product[ialt, iact]['XCH4 proxy'] = float("nan")
                l2product[ialt, iact]['XCO2 proxy precision'] = float("nan")
                l2product[ialt, iact]['XCH4 proxy precision'] = float("nan")

            l2product[ialt, iact]['spec_shift']   = 0.  #dummies for the time beeing
            l2product[ialt, iact]['spec_squeeze'] = 0.

            # XCO2[iact] = l2product[ialt, iact]['XCO2 proxy']*1.E6
            # XCO2_prec[iact] = l2product[ialt, iact]['XCO2 proxy precision']*1.E6
            # XCO2_true_smoothed[iact] = np.dot(l2product[ialt, iact]['XCO2 col avg kernel'],
            #                                   atm_sgm['dcol_co2'][ialt, iact, :])/np.sum(atm.air)*1.E6
            # XCO2_true[iact] = np.sum(atm_sgm['dcol_co2'][ialt, iact, :])/np.sum(atm.air)*1.E6

    # output to netcdf file
    level2_output(config['io_files']['output_l2'], l2product, retrieval_init, l1b, config['retrieval_init'])

#    if(sw_diag_output):
#        level2_diags_output(config['l2_diag'], l2product, measurement)

    print('=> l1bl2 finished successfully')

#    print('cumulative run time: ',runtime_cum)

    return l2product

def level2_nan(retrieval_init, nlay):

    output = {}
    output['chi2'] = float("nan")
    output['convergence'] = False
    output['number_iter'] = 0.

    for l, key in enumerate(retrieval_init['trace gases'].keys()):
        output['X'+key] = float("nan")
        output['X'+key+' precision'] = float("nan")
        output['gain_X'+key] = float("nan")

    m = 0
    while 'alb%d' % (m) in retrieval_init['surface'].keys():
        output['alb%d' % (m)] = float("nan")
        m = m+1
    output['XCO2 col avg kernel'] = np.zeros(nlay)
    output['XCH4 col avg kernel'] = np.zeros(nlay)
    output['XH2O col avg kernel'] = np.zeros(nlay)

    return(output)
    
if __name__ == '__main__':
    config = yaml.safe_load(open(sys.argv[1]))
    level1b_to_level2_processor(config)
