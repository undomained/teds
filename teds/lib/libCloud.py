import numpy as np
import tqdm
import logging
import netCDF4 as nc
import scipy 

from teds.lib import constants
from lib.libWrite import writevariablefromname
from teds.lib import libAMF

logger = logging.getLogger('E2E')

def read_atm(file_atm, slice_alt, slice_act):
    atm = {}
    with nc.Dataset(file_atm) as f:
        for key in f.variables.keys():
            atm[key] = f[key][slice_alt,slice_act]
    return atm


def read_doas(file_doas, slice_alt, slice_act):
    doas = {}
    with nc.Dataset(file_doas) as f:
        doas['o2o2_scd'] = f['doas/o2o2/oxygen_oxygen_dimer_slant_column_density'][slice_alt,slice_act]
        doas['o2o2_R'] = f['doas/o2o2/O2O2_continuum_reflectance_at_reference_wavelength'][slice_alt,slice_act]
        doas['lat'] = f['lat'][slice_alt,slice_act]
        doas['lon'] = f['lon'][slice_alt,slice_act]
        doas['sza'] = f['sza'][slice_alt,slice_act]
        doas['vza'] = f['vza'][slice_alt,slice_act]
        doas['saa'] = f['saa'][slice_alt,slice_act]
        doas['vaa'] = f['vaa'][slice_alt,slice_act]
    return doas

def get_cloud_fraction(cfg, doas, atm):
    # placeholder for cloud classification code
    # for now take cloud fraction from SGM

    cloud_results = {}
    cloud_results['cloud_fraction'] = atm['cloud_fraction'].copy()
    return cloud_results

def get_cloud_parameters(cfg, doas, atm, cld_results):
    # retrieve scattering cloud parameters
    # for now: assume binary cloud mask, no cloud fraction

    ### algo settings
    cld_pres_first_guess = 500 # [hPa] -  first guess of the cloud bottom pressure before iteration
    max_iter = 2 # max number of iterations
    wvl_o2o2_amf = 477.0 # nm
    lapserate = -0.0065 # [K/m]
    cld_thickness = 1000.0 # [m] - fixed from AMF radiative transfer

    # init output
    dictnames = ['cot','cld_top_pres','cld_bot_pres']

    for name in dictnames:
        cld_results[name] = np.ma.masked_all_like(doas['lat'])

    # load NN
    amf_scm_NN = libAMF.read_NN('LUT_AMF_SCM', cfg['LUT_NN_file'])
    cot_NN = libAMF.read_NN('LUT_COT_SCM', cfg['LUT_NN_file'])

    # convert input
    mu = np.cos(np.deg2rad(doas['vza']))
    mu0 = np.cos(np.deg2rad(doas['sza']))
    dphi	= np.abs( 180.0 - np.abs(doas['vaa']-doas['saa'])) # RAA [deg]

    play = atm['pressure_layers'][:,:,::-1] # hPa
    tlay = atm['temperature'][:,:,::-1] # K
    surface_pressure = atm['pressure_levels'][:,:,-1] # hPa


    cld_results['cld_bot_pres'] = np.ones_like(cld_results['cld_bot_pres'])*cld_pres_first_guess

    # iterate because of cot dependency on cloud pressure
    for i in range(max_iter):

        # NN for COT (cloud optical thickness)

        vector_cot_NN = np.array([surface_pressure, atm['albedo'], doas['o2o2_R'], cld_results['cld_bot_pres'], mu0, mu, dphi])

        cld_results['cot'] = libAMF.predict_NN_vector_2D(vector_cot_NN, cot_NN)

        # apply NN bounds and mask
        cld_results['cot'][cld_results['cot']< 0.0] = 0.0
        cld_results['cot'][cld_results['cot']> 256.0] = 256.0
        cld_results['cot'] = np.ma.masked_where(cld_results['cloud_fraction']== 0.0, cld_results['cot'])
        cld_results['cot'] = np.ma.masked_invalid(cld_results['cot'])

        # NN for cloud bottom pressure - iterate over pixels, as each pixel has a cloud pressure iteration

        iterlist = tqdm.tqdm(np.ndindex(doas['lat'].shape), total=doas['lat'].size)
        for idx,idy in iterlist:
            
            # skip if cot did not converge
            if np.ma.is_masked(cld_results['cot'][idx,idy]):
                continue

            #### using bisection method.
            point_cloud = [wvl_o2o2_amf, surface_pressure[idx,idy], atm['albedo'][idx,idy], cld_results['cot'][idx,idy], 0, mu0[idx,idy], mu[idx,idy], dphi[idx,idy] , 0]
            scd_O2O2_cld_list = bisection_o2o2(amf_scm_NN, point_cloud, play[idx,idy], tlay[idx,idy], doas['o2o2_scd'][idx,idy])

            # filter out nans
            interp_subset = np.isnan(scd_O2O2_cld_list)==0

            if interp_subset.sum() == 0:
                cld_results['cld_bot_pres'][idx,idy] = np.ma.masked
                cld_results['cld_top_pres'][idx,idy] = np.ma.masked
                break

            # # interpolate scd_O2O2_cld_list to find cloud pressure

            cld_pres_interpolator = scipy.interpolate.interp1d(scd_O2O2_cld_list[interp_subset],play[idx,idy,interp_subset],bounds_error=False, fill_value=(np.nan,np.nan))
            cld_results['cld_bot_pres'][idx,idy] = cld_pres_interpolator(doas['o2o2_scd'][idx,idy]) #[hPa]

            # only in last iteration
            if np.isnan(cld_results['cld_bot_pres'][idx,idy]) == False and i==max_iter-1:
                # calculate cloud top pressure 
                # option 1
                # use barometric equation to calculate cld_pres + (cloud thickness in m)
                # assume cloud is in troposphere:

                # create PT spline, use log pressure, reverse levels
                pt_spline = scipy.interpolate.InterpolatedUnivariateSpline(np.log(play[idx,idy,::-1]),tlay[idx,idy,::-1], ext='const', k=4)

                T_cldpres = pt_spline(np.log(cld_results['cld_bot_pres'][idx,idy]))

                cld_results['cld_top_pres'][idx,idy] = cld_results['cld_bot_pres'][idx,idy] * \
                    ((T_cldpres+ cld_thickness*lapserate)/T_cldpres)**(-constants.g0*constants.MDRYAIR/(constants.Rgas*lapserate))

    return cld_results


def bisection_o2o2(amf_NN,input_vector,p,T,scd_O2O2_cld_meas):
    # bisection method for strictly increasing relation
    # works well for the o2o2_scd - cld_pres  relation

    # for very low COT (<3) o2o2_scd vs cld_pres curve not always increasing, sometimes c-shape, then there are multiple solutions.
    # --> algorithm becomes unstable (not enough cloud to fit)


    # method searches the cloud pressure list with bisection method
    # goal is to find the cloud pressures that result in modelled scd_o2o2_cld closest to measured scd_O2O2_cld
    # then linear interpolation on the closest pressures to get the final cloud pressure

    # constants
    O2_vmr = 0.21 # volume mixing ratio of O2

    amf_geo  = 1/input_vector[-4] + 1/input_vector[-3]

    amf_cld_list = np.zeros(len(p))
    scd_O2O2_cld_list = np.full(len(p),np.nan)

    # only consider cloud pressures inside cloud pressure LUT/NN
    p_cld = p[p>111]

    check=True

    # start loop
    while True:

        # check if scd is within range
        if check:
            k=0
            interval = 99

        if k<0 or k>=len(p_cld):
            break

        # set cloud pressure
        input_vector[-5] = p_cld[k]

        # get amf cld from NN
        for j in range(len(p)):
            input_vector[-1] = p[j]
            amf_cld_list[j] = libAMF.predict_NN(input_vector, amf_NN)*amf_geo	
            if amf_cld_list[j] < 0.0:
                amf_cld_list[j] = 0.0	

        # integrate with trapezoidal rule
        amf_tot_cld = 0.0
        integrant_cld = np.zeros((2,))

        integrant_cld[0] = amf_cld_list[0] * p[0]*100 / T[0]

        for ilev in range(1,len(p)):
            integrant_cld[1] = amf_cld_list[ilev] * p[ilev]*100 / T[ilev]

            dp = np.abs(p[ilev] - p[ilev-1])*100

            amf_tot_cld += (integrant_cld[0] + integrant_cld[1]) * dp / 2.0

            # shift integrand to avoid double calculation.
            integrant_cld[0] = integrant_cld[1]

        # calculate O2O2 cloud column

        scd_O2O2_cld_list[k] = O2_vmr**2 * constants.Rgas / (constants.MDRYAIR*constants.g0*constants.kboltzmann**2) * amf_tot_cld / (constants.NA**2) # [mol2 m-5]
        
        # check for convergence, break loop
        if (interval==0):
            # print(counter)
            break

        # break early if scd out of range
        if check:
            if scd_O2O2_cld_list[k] < scd_O2O2_cld_meas:
                break
            else:
                check = False
                
                # start in the middle of interval
                k = len(p_cld)//2
                interval = k

                continue

        # if scd model is lower than scd meas go to lower half of interval
        if scd_O2O2_cld_list[k] < scd_O2O2_cld_meas:
            k -= int(np.ceil(interval/2))

        # if scd model is higher than scd meas go to upper half of interval
        elif scd_O2O2_cld_list[k] > scd_O2O2_cld_meas:
            k += int(np.ceil(interval/2))

        # half the interval
        interval =  interval//2


    return scd_O2O2_cld_list


def write_cloud(cfg, cloud, slice_alt, slice_act):

    
    def write_out(var):

        if cloud[var].ndim == 2:
            dim = ('scanline','ground_pixel')
            out = np.ma.masked_all_like(dst['lat'])
            out[slice_alt,slice_act] = cloud[var]
        elif cloud[var].ndim == 3:
            dim = ('scanline','ground_pixel','pressure_layers')
            out = np.ma.masked_all(dst['lat'].shape+(cloud[var].shape[-1],))
            out[slice_alt,slice_act,:] = cloud[var]
        else:
            logging.error('{var} has {var.ndim} dimensions, not recognised.')

        _ = writevariablefromname(dst, vardict[var], dim, out)
        
        return


    vardict = {'cloud_fraction':'cloud_fraction',
               'cot':'cloud_optical_thickness',
               'cld_top_pres':'cloud_top_pressure',
               'cld_bot_pres':'cloud_bottom_pressure'}
    

    with nc.Dataset(cfg['io']['l2'], 'a') as dst:

        for var in vardict:
            write_out(var)

    return