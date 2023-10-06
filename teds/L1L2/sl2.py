import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import yaml
from tqdm import tqdm
from copy import deepcopy
import netCDF4 as nc

from modules.lib import libNumTools
#from modules.lib import libRT
from modules.lib import libATM
#from modules.lib import libINV

# simplified L2 processor

def get_gm_data(path, filename):

    file = path+filename+'.nc'
    input = nc.Dataset(file, mode='r')

    gm_data = {}
    gm_data['sza'] = deepcopy(input['sza'][:, :])
    gm_data['saa'] = deepcopy(input['saa'][:, :])
    gm_data['vza'] = deepcopy(input['vza'][:, :])
    gm_data['vaa'] = deepcopy(input['vaa'][:, :])
    gm_data['lat'] = deepcopy(input['lat'][:, :])
    gm_data['lon'] = deepcopy(input['lon'][:, :])

    input.close()

    return(gm_data)

# get geometry data

project_path = '/home/jochen/TANGO_E2ES/EndtoEndProject/' 
gm_filen = 'Tango_Carbon_gm_S2_microHH_run0001'
gm_path = project_path + 'data/interface_data/gm/'

gm_data = get_gm_data(gm_path, gm_filen)

libNumTools.print_attributes(gm_data)
microHH_path = project_path + 'data/microHH/Jaenschwalde_simulation1/'

config_microHH = {}
config_microHH['time_stamp']= '0036000'
config_microHH['lat_lon_src'] =  [51.83472, 14.46277]                   #latitude, longitude coordinates for source in microHH

config_kernel_parameter={}                                              #parameter to configure the 2D blurring of the scene
config_kernel_parameter['type']= '2D Gaussian'                          #type pf kernel, current version has only one option  
config_kernel_parameter['fwhm_x']=  300                                 #FWHM in ACT direction
config_kernel_parameter['fwhm_y']=  300                                 #FWHM in flight direction
config_kernel_parameter['size_factor']= 2                               #the size of the convolution kernel in units of fwhm

microHH = libATM.get_microHH_atm(gm_data['lat'], gm_data['lon'], microHH_path,
                                         config_microHH,
                                         config_kernel_parameter)

AFGL_path = project_path + 'data/AFGL/'

nlay= 30                    #number of atmospheric layers
dzlay= 1000                 #geometrical thickness of the layers [m]

atm_std = libATM.get_AFGL_atm_homogenous_distribtution(AFGL_path, nlay, dzlay)

atm = libATM.combine_microHH_standard_atm(microHH, atm_std)

pickle.dump(atm, open(path_data+'microHH_plot.pkl', 'wb'))



