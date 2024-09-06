import numpy as np
import tqdm
import netCDF4 as nc

from teds import log
from teds.lib import constants
from teds.lib.libWrite import writevariablefromname


def read_atm(file_atm, slice_alt, slice_act):
    atm = {}
    with nc.Dataset(file_atm) as f:
        for key in f.variables.keys():
            atm[key] = f[key][slice_alt,slice_act]
    return atm

def read_doas(file_doas, slice_alt, slice_act):
    doas = {}
    with nc.Dataset(file_doas) as f:
        doas['no2_scd'] = f['doas/no2/nitrogendioxide_slant_column_density'][slice_alt,slice_act]
        doas['no2_scd_error'] = f['doas/no2/nitrogendioxide_slant_column_density_precision'][slice_alt,slice_act]
        doas['lat'] = f['lat'][slice_alt,slice_act]
        doas['lon'] = f['lon'][slice_alt,slice_act]
        doas['sza'] = f['sza'][slice_alt,slice_act]
        doas['vza'] = f['vza'][slice_alt,slice_act]
        doas['saa'] = f['saa'][slice_alt,slice_act]
        doas['vaa'] = f['vaa'][slice_alt,slice_act]
    return doas

def read_cloud(file_cloud, slice_alt, slice_act):
    cloud = {}
    with nc.Dataset(file_cloud) as f:
        if 'cloud_fraction' in f.variables:
            cloud['cot'] = f['cloud_optical_thickness'][slice_alt,slice_act]
            cloud['cot_error_alb'] = f['cloud_optical_thickness_error_albedo'][slice_alt,slice_act]
            cloud['cot_error_R'] = f['cloud_optical_thickness_error_reflectance'][slice_alt,slice_act]
            cloud['cloud_bottom_pressure'] = f['cloud_bottom_pressure'][slice_alt,slice_act]
            cloud['cloud_bottom_pressure_error'] = f['cloud_bottom_pressure_error'][slice_alt,slice_act]
            cloud['cloud_fraction'] = f['cloud_fraction'][slice_alt,slice_act]
        else:
            cloud['cloud_fraction'] = np.zeros_like(f['lat'][slice_alt,slice_act])
    return cloud

def get_amf(cfg, doas, atm, cloud):
    # -----------------------------------------------------------------
    # - Calculate NO2 AMF using AMF LUT NN
    # - Vectorised
    # - TODO:
    #   - trop-strat separation
    #   - cloud radiance fraction
    #   - albedo error estimate
    # -----------------------------------------------------------------

    #### placeholder for albedo error
    error_albedo = 0.015 # assumed albedo error
    ###

    # create output fields
    results = {}
    dictnames = ['no2_total_amf','no2_total_amf_error',
                 'no2_total_vcd','no2_total_vcd_error',
                 'no2_total_scd','no2_total_scd_error' ]

    for name in dictnames:
        results[name] = np.ma.masked_all_like(doas['lat'])

    results['no2_averaging_kernel'] = np.ma.masked_all_like(atm['temperature'])
    results['pressure_layers'] = atm['pressure_layers'][:,:,::-1] # hPa

    # load NN
    amf_clear_NN = read_NN('LUT_AMF_clear', cfg['LUT_NN_file'])
    amf_cloud_NN = read_NN('LUT_AMF_SCM', cfg['LUT_NN_file'])

    # convert input
    mu = np.cos(np.deg2rad(doas['vza']))
    mu0 = np.cos(np.deg2rad(doas['sza']))
    dphi	= np.abs( 180.0 - np.abs(doas['vaa']-doas['saa'])) # RAA [deg]

    pressure_levels_midpoint = results['pressure_layers']
    surface_pressure = atm['pressure_levels'][:,:,-1] # hPa

    no2_profile =  atm['dcol_no2'][:,:,::-1] /constants.NA * 1e4 # [molec/cm2] to [mol/m2]

    # get temperature correction
    f_tcorr = get_tcorr(atm['temperature'][:,:,::-1])

    # calculate tropopause layer index
    # tropopause_layer_index = get_tropopause(atm['temperature'][iscan,ipxl,:])

    # NN LUT 
    wvl_no2_amf = 437.5 # nm
    wvl = np.ones_like(doas['lat'])*wvl_no2_amf


    # clear-sky
    vector_clear =  np.array([wvl, surface_pressure, atm['albedo_B01'], mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_clear = get_boxamf_NN(vector_clear, amf_clear_NN, pressure_levels_midpoint, doas, mu, mu0)

    # skip if no cloudy pixels
    if cloud['cloud_fraction'].sum()>0.0:

        # cloud
        vector_cloud =  np.array([wvl, surface_pressure, atm['albedo_B01'], cloud['cot'],
                                  cloud['cloud_bottom_pressure'], mu0, mu, dphi, np.zeros_like(doas['lat']) ])
        boxamf_cloud = get_boxamf_NN(vector_cloud, amf_cloud_NN, pressure_levels_midpoint, doas, mu, mu0)

    else:
        boxamf_cloud = np.zeros_like(boxamf_clear)

    # troposphere only profile:
    # no2_profile_trop = np.copy(no2_profile)
    # no2_profile_trop[(tropopause_layer_index+1):] = 0

    # stratosphere only profile:
    # no2_profile_strat = np.copy(no2_profile)
    # no2_profile_strat[:(tropopause_layer_index+1)] = 0

    results['no2_total_scd'] = doas['no2_scd']
    results['no2_total_scd_error'] = doas['no2_scd_error']

    # calculate total amf
    results['no2_total_amf_clear']= np.sum(boxamf_clear*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)
    results['no2_total_amf_cloud']= np.sum(boxamf_cloud*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)


    results['no2_total_amf'] = cloud['cloud_fraction']*results['no2_total_amf_cloud'] + (1 - cloud['cloud_fraction']) * results['no2_total_amf_clear']

    # calculate scd total by dividing by total amf
    results['no2_total_vcd']= results['no2_total_scd'] / results['no2_total_amf']

    # results_scan['tropopause_layer_index'][iscan,ipxl] = tropopause_layer_index

    # averaging kernel:  ak =  box_amf * temperature_correction / total_amf
    results['no2_averaging_kernel'] = (cloud['cloud_fraction'][...,np.newaxis]*boxamf_cloud + (1-cloud['cloud_fraction'][...,np.newaxis])*boxamf_clear)* f_tcorr  / results['no2_total_amf'][...,np.newaxis]



    ##### error propagation

    # boxamf clear with +/- albedo error

    vector_clear_error_alb_min =  np.array([wvl, surface_pressure, np.clip(atm['albedo_B01']-error_albedo,a_min=0, a_max=None), mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_clear_error_alb_min = get_boxamf_NN(vector_clear_error_alb_min, amf_clear_NN, pressure_levels_midpoint, doas, mu, mu0)

    vector_clear_error_alb_plus =  np.array([wvl, surface_pressure, np.clip(atm['albedo_B01']-error_albedo,a_min=0, a_max=None), mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_clear_error_alb_plus = get_boxamf_NN(vector_clear_error_alb_plus, amf_clear_NN, pressure_levels_midpoint, doas, mu, mu0)

    amf_total_clear_error_alb_min = np.sum(boxamf_clear_error_alb_min*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)
    amf_total_clear_error_alb_plus = np.sum(boxamf_clear_error_alb_plus*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)

    # boxamf cloud with +/- albedo,R and scd errors.

    # albedo error
    vector_cloud_error_alb_min =  np.array([wvl, surface_pressure, np.clip(atm['albedo_B01']-error_albedo,a_min=0, a_max=None),
                                            np.clip(cloud['cot']-cloud['cot_error_alb'],a_min=0, a_max=256),
                                            cloud['cloud_bottom_pressure'],
                                            mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_cloud_error_alb_min = get_boxamf_NN(vector_cloud_error_alb_min, amf_cloud_NN, pressure_levels_midpoint, doas, mu, mu0)

    vector_cloud_error_alb_plus =  np.array([wvl, surface_pressure, np.clip(atm['albedo_B01']-error_albedo,a_min=0, a_max=None),
                                            np.clip(cloud['cot']+cloud['cot_error_alb'],a_min=0, a_max=256),
                                            cloud['cloud_bottom_pressure'], mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_cloud_error_alb_plus = get_boxamf_NN(vector_cloud_error_alb_plus, amf_cloud_NN, pressure_levels_midpoint, doas, mu, mu0)

    # reflectance error
    vector_cloud_error_R_min =  np.array([wvl, surface_pressure, atm['albedo_B01'],
                                          np.clip(cloud['cot']-cloud['cot_error_R'],a_min=0, a_max=256),
                                          cloud['cloud_bottom_pressure'], mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_cloud_error_R_min = get_boxamf_NN(vector_cloud_error_R_min, amf_cloud_NN, pressure_levels_midpoint, doas, mu, mu0)

    vector_cloud_error_R_plus =  np.array([wvl, surface_pressure, atm['albedo_B01'],
                                           np.clip(cloud['cot']+cloud['cot_error_R'],a_min=0, a_max=256),
                                           cloud['cloud_bottom_pressure'],
                                           mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_cloud_error_R_plus = get_boxamf_NN(vector_cloud_error_R_plus, amf_cloud_NN, pressure_levels_midpoint, doas, mu, mu0)

    # SCD error
    vector_cloud_error_scd_min =  np.array([wvl, surface_pressure, atm['albedo_B01'],
                                          cloud['cot'],
                                          cloud['cloud_bottom_pressure']-cloud['cloud_bottom_pressure_error'], mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_cloud_error_scd_min = get_boxamf_NN(vector_cloud_error_scd_min, amf_cloud_NN, pressure_levels_midpoint, doas, mu, mu0)

    vector_cloud_error_scd_plus =  np.array([wvl, surface_pressure, atm['albedo_B01'],
                                           cloud['cot'],
                                           cloud['cloud_bottom_pressure']+cloud['cloud_bottom_pressure_error'],
                                           mu0, mu, dphi, np.zeros_like(doas['lat']) ])
    boxamf_cloud_error_scd_plus = get_boxamf_NN(vector_cloud_error_scd_plus, amf_cloud_NN, pressure_levels_midpoint, doas, mu, mu0)


    # take average of + and - error
    amf_total_cloudy_error_alb_min = np.sum(boxamf_cloud_error_alb_min*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)
    amf_total_cloudy_error_alb_plus = np.sum(boxamf_cloud_error_alb_plus*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)
    
    amf_total_cloudy_error_R_min = np.sum(boxamf_cloud_error_R_min*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)
    amf_total_cloudy_error_R_plus = np.sum(boxamf_cloud_error_R_plus*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)
    
    amf_total_cloudy_error_scd_min = np.sum(boxamf_cloud_error_scd_min*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)
    amf_total_cloudy_error_scd_plus = np.sum(boxamf_cloud_error_scd_plus*no2_profile*f_tcorr, axis=-1) / np.sum(no2_profile, axis=-1)

    amf_total_error_alb_min = results['no2_total_amf'] - (cloud['cloud_fraction']*amf_total_cloudy_error_alb_min+ (1 - cloud['cloud_fraction']) *amf_total_clear_error_alb_min )
    amf_total_error_alb_plus = results['no2_total_amf'] - (cloud['cloud_fraction']*amf_total_cloudy_error_alb_plus+ (1 - cloud['cloud_fraction']) *amf_total_clear_error_alb_plus )
    amf_total_error_alb = (np.abs(amf_total_error_alb_min)+np.abs(amf_total_error_alb_plus))/2

    amf_total_error_R_min = results['no2_total_amf'] - (cloud['cloud_fraction']*amf_total_cloudy_error_R_min+ (1 - cloud['cloud_fraction']) *results['no2_total_amf_clear'] )
    amf_total_error_R_plus = results['no2_total_amf'] - (cloud['cloud_fraction']*amf_total_cloudy_error_R_plus+ (1 - cloud['cloud_fraction']) *results['no2_total_amf_clear'] )
    amf_total_error_R = (np.abs(amf_total_error_R_min)+np.abs(amf_total_error_R_plus))/2

    amf_total_error_scd_min = results['no2_total_amf'] - (cloud['cloud_fraction']*amf_total_cloudy_error_scd_min+ (1 - cloud['cloud_fraction']) *results['no2_total_amf_clear'] )
    amf_total_error_scd_plus = results['no2_total_amf'] - (cloud['cloud_fraction']*amf_total_cloudy_error_scd_plus+ (1 - cloud['cloud_fraction']) *results['no2_total_amf_clear'] )
    amf_total_error_scd = (np.abs(amf_total_error_scd_min)+np.abs(amf_total_error_scd_plus))/2


    # total amf total error
    results['no2_total_amf_error'] = np.sqrt(amf_total_error_alb**2 + amf_total_error_R**2 + amf_total_error_scd**2 )
    
    # total vcd error
    results['no2_total_vcd_error'] = np.sqrt( (results['no2_total_scd_error']/results['no2_total_amf'])**2 + 
                                                     (results['no2_total_amf_error']*results['no2_total_vcd']/results['no2_total_amf'])**2 )

    # copy mask
    results['no2_total_amf_error'].mask +=  results['no2_total_vcd_error'].mask

    return results


def get_amf_iter(cfg, doas, atm):
    # -----------------------------------------------------------------
    # Calculate NO2 AMF using AMF LUT NN
    # iterate over pixels
    # no clear-sky
    # -----------------------------------------------------------------

    results = {}

    # create output fields
    dictnames = ['no2_total_amf','no2_total_vcd','no2_total_scd']

    for name in dictnames:
        results[name] = np.ma.masked_all_like(doas['lat'])

    results['no2_averaging_kernel'] = np.ma.masked_all_like(atm['temperature'])
    results['pressure_layers'] = atm['pressure_layers'][idx,idy,::-1] # hPa

    # load NN
    amf_clear_NN = read_NN('LUT_AMF_clear', cfg['LUT_NN_file'])



    iterlist = tqdm.tqdm(np.ndindex(doas['lat'].shape), total=doas['lat'].size)
    for idx,idy in iterlist:

        if np.ma.is_masked(doas['no2_scd'][idx,idy]):
            log.info(f'Skipping pixel {idx},{idy}: NaN in doas input')
            continue

        mu = np.cos(np.deg2rad(doas['vza'][idx,idy]))
        mu0 = np.cos(np.deg2rad(doas['sza'][idx,idy]))
        dphi	= np.abs( 180.0 - np.abs(doas['vaa'][idx,idy]-doas['saa'][idx,idy])) # RAA [deg]

        pressure_levels_midpoint = results['pressure_layers']
        surface_pressure = atm['pressure_levels'][idx,idy,-1] # hPa

        no2_profile =  atm['dcol_no2'][idx,idy,::-1] /constants.NA * 1e4 # [molec/cm2] to [mol/m2]

        # get temperature correction
        cl = get_tcorr(atm['temperature'][idx,idy,::-1])

        # calculate tropopause layer index
        # tropopause_layer_index = get_tropopause(atm['temperature'][idx,idy,:])

        # NN LUT 
        point_clear =  [437.5, surface_pressure, atm['albedo_B01'][idx,idy], mu0, mu, dphi, 0]

        # loop over pressure levels and get boxamf clear and cloudy
        boxamf_clear = np.zeros_like(pressure_levels_midpoint)

        amf_geo = 1/mu + 1/mu0


        for j in range(len(pressure_levels_midpoint)):
            # clear
            point_clear[-1] = pressure_levels_midpoint[j]
            boxamf_clear[j] = predict_NN(point_clear, amf_clear_NN)*amf_geo

        boxamf_clear[boxamf_clear<0.0] = 0.0

        # troposphere only profile:
        # no2_profile_trop = np.copy(no2_profile)
        # no2_profile_trop[(tropopause_layer_index+1):] = 0

        # stratosphere only profile:
        # no2_profile_strat = np.copy(no2_profile)
        # no2_profile_strat[:(tropopause_layer_index+1)] = 0

        results['no2_total_scd'][idx,idy] = doas['no2_scd'][idx,idy]

        # calculate total amf
        results['no2_total_amf'][idx,idy] = np.sum(boxamf_clear*no2_profile*cl) / np.sum(no2_profile)

        # calculate scd total by dividing by total amf
        results['no2_total_vcd'][idx,idy] = results['no2_total_scd'][idx,idy] / results['no2_total_amf'][idx,idy]

        # results['tropopause_layer_index'][idx,idy] = tropopause_layer_index

        # averaging kernel:  ak =  box_amf * temperature_correction / total_amf

        results['no2_averaging_kernel'][idx,idy,:] = boxamf_clear * cl  / results['no2_total_amf'][idx,idy]

        # log.info('Processed pixel alt {}/{} act {}/{} in {}s'.format(idx,doas['lat'].shape[0],idy,doas['lat'].shape[1],np.round((time.time() - start_time_pixel),2) ))
    return results

def get_tcorr(t):
    # --------------------------
    # Calculate temperature correction factor (see tropomi atbd)
    # ---------------------------

    # temperature correction factor
    t_sigma = 220 #[K] no2 cross-section spectrum temperature

    f_tcorr = 1 - 0.00316* (t - t_sigma) + 3.39E-6*(t - t_sigma)**2

    return f_tcorr

def get_tropopause(temperature, geometric_layer_height):
    # --------------------------
    # # find tropopause. Use WMO tropopause definition (1957) :
    # condition 1: temperature gradient <= 2 K/ km.
    # condition 2: temperaturegradient next layers within 2 km height also <= 2K/km

    # -- > setting condition for individual layers above selected layer results in very spatially non-smooth tropopause layers
    # ---> TM5 might use chemical tropopause definition based on O3, not available
    # ---> use alternative interpretation: mean lapse rate of layers within +2km height below 2K/km
    # ---> still patches on map.. other solution required.
    # ---------------------------
    #
    # TM5 uses smoothing in offl mode, see ./tropomi/TM.py: smooth_ltropo()

    temp_grad = np.abs( np.diff(temperature) / np.diff(geometric_layer_height) )

    tropopause_layer_index = -1
    for i in range(len(temp_grad)):
        if temp_grad[i] <= 2.0:
            delta_height = np.cumsum(np.diff(geometric_layer_height)[i:])
            # print(temp_grad[i:][delta_height<=2.0])

            # if (temp_grad[i:][delta_height<=2.0] > 2.0).any():
            # 	# print('False positive: ',i)
            # 	continue

            # alternative:
            if (temp_grad[i:][delta_height<=2.0].mean()> 2.0):
                # print('False positive: ',i)
                continue
            else:
                # print('Tropopause found; layer: ',i, ' , layer height: ',geometric_layer_height[i],' km')
                tropopause_layer_index = i
                break

    if tropopause_layer_index == -1:
        print('error: tropopause not found')

    return tropopause_layer_index

def read_NN(parameter, NN_file):
    '''
    Read Neural Network file and save as Python dict

    Input:
    - parameter: name of the LUT NN
    - NN_file: path to file

    Output:
    - NN model saved as dict.
    '''

    dict_out = {}

    f_NN = nc.Dataset(NN_file, 'r')

    dict_out['normalization_input'] = f_NN[f'{parameter}/normalization_input'][:]

    for group in f_NN[parameter].groups.keys():
        dict_out[group] = {'bias':f_NN[f'{parameter}/{group}/bias'][:], 'kernel':f_NN[f'{parameter}/{group}/kernel'][:]}

    # check
    if 'layer_4' not in dict_out:
        log.error('All NNs should have 4 layers')
        raise

    return dict_out


def leakyrelu(x,alpha=0.01):
    # leaky relu function
    # https://keras.io/api/layers/activation_layers/leaky_relu/
        return np.where(x >=0 , x, x * alpha)   

def predict_NN(input_vector, NN):
    '''
    Normalize input vector and apply NN model.

    Input:
    - input_vector: input vector for NN. Variables specified in NN file.
    - NN: NN dict
    
    Output:
    - Output variable. Specified in NN file.
    '''

    # normalize input
    norm = NN['normalization_input'][()]
    input_vector_norm = (input_vector - norm[:,0]) / (norm[:,1] - norm[:,0])
    
    # apply NN
    layer1 = leakyrelu(np.dot(input_vector_norm,NN['layer_1']['kernel']) + NN['layer_1']['bias'])
    layer2 = leakyrelu(np.dot(layer1,NN['layer_2']['kernel']) + NN['layer_2']['bias'])
    layer3 = leakyrelu(np.dot(layer2,NN['layer_3']['kernel']) + NN['layer_3']['bias'])
    output= np.dot(layer3,NN['layer_4']['kernel']) + NN['layer_4']['bias'] # output layer uses linear activation function

    return output[0]

def predict_NN_vector_3D(input_vector, NN):

    # normalize input
    norm = NN['normalization_input'][()]

    input_vector_norm = (input_vector - norm[:,0,np.newaxis,np.newaxis,np.newaxis]) / (norm[:,1,np.newaxis,np.newaxis,np.newaxis] - norm[:,0,np.newaxis,np.newaxis,np.newaxis])

    input_vector_norm = np.moveaxis(input_vector_norm, 0, -1)

    # optimize='optimal'
    layer1 = leakyrelu( np.einsum('ijkl,lm->ijkm',input_vector_norm,NN['layer_1']['kernel']) + NN['layer_1']['bias'][np.newaxis,np.newaxis,np.newaxis,...] )
    layer2 = leakyrelu( np.einsum('ijkl,lm->ijkm',layer1,NN['layer_2']['kernel']) + NN['layer_2']['bias'][np.newaxis,np.newaxis,np.newaxis,...] )
    layer3 = leakyrelu( np.einsum('ijkl,lm->ijkm',layer2,NN['layer_3']['kernel']) + NN['layer_3']['bias'][np.newaxis,np.newaxis,np.newaxis,...] )
    output = np.einsum('ijkl,lm->ijkm',layer3,NN['layer_4']['kernel']) + NN['layer_4']['bias'][np.newaxis,np.newaxis,np.newaxis,...]

    return output[:,:,:,0]

def predict_NN_vector_2D(input_vector, NN):

    # normalize input
    norm = NN['normalization_input'][()]

    input_vector_norm = (input_vector - norm[:,0,np.newaxis,np.newaxis]) / (norm[:,1,np.newaxis,np.newaxis] - norm[:,0,np.newaxis,np.newaxis])

    input_vector_norm = np.moveaxis(input_vector_norm, 0, -1)

    # optimize='optimal'
    layer1 = leakyrelu( np.einsum('ijl,lm->ijm',input_vector_norm,NN['layer_1']['kernel']) + NN['layer_1']['bias'][np.newaxis,np.newaxis,...] )
    layer2 = leakyrelu( np.einsum('ijl,lm->ijm',layer1,NN['layer_2']['kernel']) + NN['layer_2']['bias'][np.newaxis,np.newaxis,...] )
    layer3 = leakyrelu( np.einsum('ijl,lm->ijm',layer2,NN['layer_3']['kernel']) + NN['layer_3']['bias'][np.newaxis,np.newaxis,...] )
    output = np.einsum('ijl,lm->ijm',layer3,NN['layer_4']['kernel']) + NN['layer_4']['bias'][np.newaxis,np.newaxis,...]

    return output[:,:,0]


def get_boxamf_NN(vector, NN, pressure_levels, doas, mu, mu0):

    # repear vector along pressure dim
    vector = np.tile(vector[...,np.newaxis],pressure_levels.shape[-1])
    vector[-1,:,:,:] = pressure_levels

    # try calculating amf vectorised, otherwise loop over pixels
    try:
        boxamf= predict_NN_vector_3D(vector, NN)
    except MemoryError as e:
        log.error(e)
        log.error('Not enough memory for vectorised approach, falling back to pixel by pixel approach')
        boxamf = np.ma.masked_all_like(pressure_levels)
        iterlist = tqdm.tqdm(np.ndindex(doas['lat'].shape), total=doas['lat'].size)
        for idx,idy in iterlist:
            for idz in range(pressure_levels.shape[-1]):
                boxamf[idx,idy,idz]= predict_NN(vector[:,idx,idy,idz], NN)

    # multiply with amf geo
    amf_geo = 1/mu + 1/mu0
    boxamf*= amf_geo[...,np.newaxis]
    boxamf[boxamf<0.0] = 0.0

    return boxamf

def write_amf(cfg, amf, slice_alt, slice_act):

    # convert units
    def molm2_to_moleccm2(var):
        var = var *constants.NA / 1e4 #  [mol/m2] to [molec/cm2]
        return var
    
    def write_out(var):

        if amf[var].ndim == 2:
            dim = ('scanline','ground_pixel')
            out = np.ma.masked_all_like(dst['lat'])
            out[slice_alt,slice_act] = amf[var]
        elif amf[var].ndim == 3:
            dim = ('scanline','ground_pixel','pressure_layers')
            out = np.ma.masked_all(dst['lat'].shape+(amf[var].shape[-1],))
            out[slice_alt,slice_act,:] = amf[var]
        else:
            log.error('{var} has {var.ndim} dimensions, not recognised.')

        _ = writevariablefromname(dst, var, dim, out)
        
        return
    
    for var in ['no2_total_vcd','no2_total_scd','no2_total_vcd_error','no2_total_scd_error']:
        amf[var] = molm2_to_moleccm2(amf[var])

    varlist = ['pressure_layers', 'no2_averaging_kernel', 'no2_total_amf', 'no2_total_vcd', 'no2_total_scd',
                'no2_total_vcd_error','no2_total_scd_error','no2_total_amf_error']

    with nc.Dataset(cfg['io']['l2'], 'a') as dst:
        
        if 'pressure_layers' not in dst.dimensions:
            p_dim = dst.createDimension('pressure_layers', amf['pressure_layers'].shape[-1])

        for var in varlist:
            write_out(var)

    return




