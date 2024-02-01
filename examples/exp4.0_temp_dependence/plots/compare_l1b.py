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
    l1b_data['wavelength'] = deepcopy(input['OBSERVATION_DATA']['wavelengths'][:])
    l1b_data['radiance']   = deepcopy(input['OBSERVATION_DATA']['radiance'][:])
    l1b_data['noise']      = deepcopy(input['OBSERVATION_DATA']['radiance_noise'][:])

    input.close()

    return(l1b_data)
  
if __name__ == "__main__":
    # end to end global config file
    
    filename1 = '/home/jochen/TANGO_E2ES/EndtoEndProject/data/interface_data/level1b/Tango_Carbon_l1b_exp4.0.nc'
    filename2 = '/home/jochen/TANGO_E2ES/EndtoEndProject/data/interface_data/level1b/Tango_Carbon_l1b_exp4.0_siml1b.nc'
    
    l1b_v1 = nc.Dataset(filename1, mode='r')
    l1b_v2 = nc.Dataset(filename2, mode='r')
    
    sza = deepcopy(l1b_v1['GEOLOCATION_DATA']['sza'][:])
    nalt = len(sza[:,0])
    nact = len(sza[0,:])
    chi2 = np.zeros((nalt,nact))
    
    fig2 = False
    fig1 = True
    if(fig1):
        ialt = 0
        iact = 7
  
        wave1 = deepcopy(l1b_v1['OBSERVATION_DATA']['wavelengths'][:])
        rad1  = deepcopy(l1b_v1['OBSERVATION_DATA']['radiance'][:])
        
        wave2 = deepcopy(l1b_v2['OBSERVATION_DATA']['wavelength'][:])
        rad2  = deepcopy(l1b_v2['OBSERVATION_DATA']['radiance'][:])
        l1b_v1.close()
        l1b_v2.close()
        
        fig = plt.figure(figsize=(15, 8), dpi=100)
        plt.subplot(1, 1, 1)
        plt.plot(wave1[iact,:], rad1[ialt,iact,:], color = 'blue', label = 'L1AL1B')
        plt.plot(wave2[iact,:], rad2[ialt,iact,:], color = 'green', label = 'simplified')
        plt.legend()
        plt.title('iact = 10')
        plt.ylabel('I [ph/(nm m2 sr s)]')    
        
        # plt.subplot(2, 1, 2)

        # delta = (l1b_v2['radiance'][ialt,iact,:] -l1b_v1['radiance'][ialt,iact,:])/l1b_v1['radiance'][ialt,iact,:]*100.

        # plt.plot(l1b_v1['wavelength'][iact,:], delta, color = 'black')
        # plt.xlabel('$\lambda$ [nm]')
        # plt.ylabel('$\delta$ I [%]')  
        # plt.ylim([-3.,5.])
#        plt.savefig('../plots/radiance_iact00.png')
        plt.show()       
        sys.exit()
    if(fig2):
        ialt = 0
        iwave = 200
        iact_index = np.arange(100)
        fig = plt.figure(figsize=(15, 8), dpi=100)
        plt.subplot(2, 1, 1)
        plt.plot(iact_index, l1b_v1['radiance'][ialt,:, iwave], color = 'blue', label = 'binning 2')
        plt.plot(iact_index, l1b_v2['radiance'][ialt,:, iwave], color = 'green', label = 'binning 1')
        plt.legend()
        plt.title(l1b_v1['wavelength'][0,iwave])
        plt.ylabel('I [ph/(nm m2 sr s)]')    
        
        plt.subplot(2, 1, 2)

        delta = (l1b_v2['radiance'][ialt,:,iwave] -l1b_v1['radiance'][ialt,:, iwave])/l1b_v1['radiance'][ialt,:, iwave]*100.

        plt.plot(iact_index, delta, color = 'black')
        plt.xlabel('ACT index')
        plt.ylabel('$\delta$ I [%]')  
        plt.ylim([-10.,10.])
#        plt.savefig('../plots/radiance_iact00.png')
        plt.show()       
        sys.exit()

        plt.savefig('./plots/chi_noisefree.png')
