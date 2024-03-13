import sys
import numpy as np
import logging
import tqdm
import time

import netCDF4 as nc
from scipy import interpolate

from lib import constants


def read_atm(file_atm):

    atm = {}

    with nc.Dataset(file_atm) as f:
        for key in f.variables.keys():
            atm[key] = f[key][:]

    return atm


def read_doas(file_doas):

    doas = {}

    with nc.Dataset(file_doas) as f:
        doas['no2_scd'] = f['doas/no2/nitrogendioxide_slant_column_density'][:]
        doas['lat'] = f['lat'][:]
        doas['lon'] = f['lon'][:]
        doas['sza'] = f['sza'][:]
        doas['vza'] = f['vza'][:]
        doas['saa'] = f['saa'][:]
        doas['vaa'] = f['vaa'][:]

    return doas

def get_amf(cfg, doas, atm):
    # -----------------------------------------------------------------
    # Calculate NO2 AMF using AMF LUT NN
    # -----------------------------------------------------------------
    # 
    results = {}

    # create output fields
    dictnames = ['amf_total','no2_total_vcd','no2_total_scd']

    for name in dictnames:
        results[name] = np.ma.masked_all_like(doas['lat'])

    results['no2_averaging_kernel'] = np.ma.masked_all_like(atm['temperature'])
    results['pressure_layer'] = atm['pressure_layers'][:,:,::-1]

    results['no2_total_vcd_sgm'] = atm['col_no2'] /constants.NA * 1e4 # [molec/cm2] to [mol/m2]

    # load NN
    amf_clear_NN = read_NN('LUT_AMF_clear', cfg['LUT_NN_file'])



    iterlist = tqdm.tqdm(np.ndindex(doas['lat'].shape), total=doas['lat'].size)
    for idx,idy in iterlist:

        if np.ma.is_masked(doas['no2_scd'][idx,idy]):
            logging.info(f'Skipping pixel {idx},{idy}: NaN in doas input')
            return

        start_time_pixel = time.time()


        mu = np.cos(np.deg2rad(doas['vza'][idx,idy]))
        mu0 = np.cos(np.deg2rad(doas['sza'][idx,idy]))
        dphi	= np.abs( 180.0 - np.abs(doas['vaa'][idx,idy]-doas['saa'][idx,idy])) # RAA [deg]

        pressure_levels_midpoint = atm['pressure_layers'][idx,idy,::-1] # hPa
        surface_pressure = atm['pressure_levels'][idx,idy,-1] # hPa

        no2_profile =  atm['dcol_no2'][idx,idy,::-1] /constants.NA * 1e4 # [molec/cm2] to [mol/m2]

        # get temperature correction
        cl = get_cl(atm['temperature'][idx,idy,::-1])

        # calculate tropopause layer index
        # tropopause_layer_index = get_tropopause(atm['temperature'][idx,idy,:])

        # NN LUT 
        point_clear =  [437.5, surface_pressure, atm['albedo'][idx,idy], mu0, mu, dphi, 0]

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
        results['amf_total'][idx,idy] = np.sum(boxamf_clear*no2_profile*cl) / np.sum(no2_profile)

        # calculate scd total by dividing by total amf
        results['no2_total_vcd'][idx,idy] = results['no2_total_scd'][idx,idy] / results['amf_total'][idx,idy]

        # results['tropopause_layer_index'][idx,idy] = tropopause_layer_index

        # averaging kernel:  ak =  box_amf * temperature_correction / total_amf

        results['no2_averaging_kernel'][idx,idy,:] = boxamf_clear * cl  / results['amf_total'][idx,idy]

        # logging.info('Processed pixel alt {}/{} act {}/{} in {}s'.format(idx,doas['lat'].shape[0],idy,doas['lat'].shape[1],np.round((time.time() - start_time_pixel),2) ))

    return results


def get_cl(t):
    # --------------------------
    # Calculate temperature correction factor (see tropomi atbd)
    # ---------------------------

    # temperature correction factor
    t_sigma = 220 #[K] no2 cross-section spectrum temperature

    cl = 1 - 0.00316* (t - t_sigma) + 3.39E-6*(t - t_sigma)**2

    return cl


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

    return dict_out


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


    # leaky relu function
    # https://keras.io/api/layers/activation_layers/leaky_relu/
    def leakyrelu(x,alpha=0.01):
        return np.where(x >=0 , x, x * alpha)   

    # check
    if 'layer_4' not in NN:
        logging.error('All NNs should have 4 layers')
        sys.exit()
    
    # apply NN
    layer1 = leakyrelu(np.dot(input_vector_norm,NN['layer_1']['kernel']) + NN['layer_1']['bias'])
    layer2 = leakyrelu(np.dot(layer1,NN['layer_2']['kernel']) + NN['layer_2']['bias'])
    layer3 = leakyrelu(np.dot(layer2,NN['layer_3']['kernel']) + NN['layer_3']['bias'])
    output= np.dot(layer3,NN['layer_4']['kernel']) + NN['layer_4']['bias'] # output layer uses linear activation function

    return output[0]


def write_amf(cfg,amf):

    with nc.Dataset(cfg['l2_file'], 'a') as dst:

        dst.amf_config = str(cfg)

        group = 'amf'
        newgroup = dst.createGroup(group)

        p_dim = dst[group].createDimension('pressure_layer', amf['pressure_layer'].shape[-1])
        
        # coord_string = "lat lon"


        def write_field(dictname,fieldname,units=None):

            if amf[dictname].ndim == 2:
                var = dst[group].createVariable(fieldname, float, ('scanline','ground_pixel'), fill_value=9.96921E36)
                # var.coordinates = coord_string
                var[:,:] =   amf[dictname]

            elif amf[dictname].ndim == 3:
                var = dst[group].createVariable(fieldname, float, ('scanline','ground_pixel','pressure_layer'), fill_value=9.96921E36)
                var[:,:,:] =  amf[dictname]

            if units:
                var.units = units
            else:
                var.units = '-'

            return

        write_field('amf_total','amf_total')
        write_field('no2_averaging_kernel','no2_averaging_kernel')
        write_field('pressure_layer','pressure_layer', units='hPa')
        write_field('no2_total_vcd','no2_total_vcd', units='mol m-2')
        write_field('no2_total_vcd_sgm','no2_total_vcd_sgm', units='mol m-2')
        write_field('no2_total_scd','no2_total_scd', units='mol m-2')

    return




