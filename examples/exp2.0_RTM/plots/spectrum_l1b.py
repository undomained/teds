#compare different l1b products
import os
import sys
import numpy as np
from copy import deepcopy
import netCDF4 as nc
import matplotlib.pyplot as plt
import yaml

def get_l1b(filename):

    #getting l1b data from file
    
    print(filename)

    input = nc.Dataset(filename, mode='r')

    l1b_data = {}
    l1b_data['sza'] = deepcopy(input['GEOLOCATION_DATA']['sza'][:])
    l1b_data['saa'] = deepcopy(input['GEOLOCATION_DATA']['saa'][:])
    l1b_data['vza'] = deepcopy(input['GEOLOCATION_DATA']['vza'][:])
    l1b_data['vaa'] = deepcopy(input['GEOLOCATION_DATA']['vaa'][:])
    l1b_data['latitude']   = deepcopy(input['GEOLOCATION_DATA']['lat'][:])
    l1b_data['longitude']  = deepcopy(input['GEOLOCATION_DATA']['lon'][:])
    l1b_data['wavelength'] = deepcopy(input['OBSERVATION_DATA']['wavelength'][:])
    l1b_data['radiance']   = deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    l1b_data['noise']      = deepcopy(input['OBSERVATION_DATA']['radiance_noise'][:])

    input.close()

    return(l1b_data)
  
if __name__ == "__main__":
    # end to end global config file
    
    filename1 = '/home/jochen/TANGO_E2ES/EndtoEndProject/data/interface_data/level1b/Tango_Carbon_l1b_exp2.0.nc'
    
    l1b_v1 = get_l1b(filename1)
    
    fig1 = True
    if(fig1):
        ialt = 0
        iact = 0
  
        fig = plt.figure(figsize=(10, 6), dpi=100)
        plt.subplot(1, 1, 1)
        plt.plot(l1b_v1['wavelength'][iact,:], l1b_v1['radiance'][ialt,iact,:], color = 'blue')
        plt.title('level 1B radiance')
        plt.ylabel('I [ph/(nm m$^2$ sr s)]')    
        plt.xlabel('wavelength [nm]')    
        plt.show()       

        plt.savefig('./plots/radiance.png')
