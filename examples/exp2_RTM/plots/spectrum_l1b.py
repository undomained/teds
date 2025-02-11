#compare different l1b products
import os
import sys
import numpy as np
from copy import deepcopy
import netCDF4 as nc
import matplotlib.pyplot as plt
import yaml
import pandas as pd

def get_sgm_rad_data(filename):
    input = nc.Dataset(filename, mode='r')
    sgm_data = {}
    sgm_data['wavelength line-by-line'] = input['wavelength'][:]
    sgm_data['solar irradiance line-by-line'] = input['solar_irradiance'][:]
    sgm_data['radiance line-by-line'] = input['radiance'][:, :, :]
    input.close()
    return(sgm_data)

if __name__ == "__main__":
    # end to end global config file
    
    filename1 = '/home/jochen/TANGO_E2ES/EndtoEndProject/data/interface_data/sgm/Tango_Carbon_sgm_radiance_exp2.0.nc'
        
    l1b = get_sgm_rad_data(filename1)
    
    conv =1.E-4 #1/m2 -> 1/cm2
    # dl1b = {'wavelength [nm]': l1b['wavelength line-by-line'],
    #         'radiance L_ref [ph./(cm2 nm sr s)]': l1b['radiance line-by-line'][0,0,:]*conv,
    #         'radiance L_{HL dark} [ph./(cm2 nm sr s)]': l1b['radiance line-by-line'][0,1,:]*conv,
    #         'radiance L_{Tr bright} [ph./(cm2 nm sr s)]': l1b['radiance line-by-line'][0,2,:]*conv,
    #         'radiance L_{ML nom} [ph./(cm2 nm sr s)]': l1b['radiance line-by-line'][0,3,:]*conv,
    #         'solar irradiance [ph./(cm2 nm s)]': l1b['solar irradiance line-by-line']*conv}
    # df = pd.DataFrame(data=dl1b)
    # filename = 'SWIR_reference.xlsx'
    # df.to_excel(filename, sheet_name='SWIR reference spectra')
    # sys.exit()
        
    filen0 = '/home/jochen/TANGO_E2ES/EndtoEndProject/data/TANGO reference spectra/TANGO_spectra/LBL_REF_alb0.15_sza70_vza0.dat'
    data = np.genfromtxt(filen0, skip_header=2)
    wave = data[:, 0]
    dark = data[:, 1]
    sun  = data[:, 2]

    fig1 = True
    if(fig1):
        ialt = 0
        iact = 0
  
        fig = plt.figure(figsize=(10, 6), dpi=100)
        plt.subplot(1, 1, 1)
#        plt.plot(l1b['wavelength line-by-line'][:], l1b['radiance line-by-line'][ialt,iact,:]*conv, color = 'blue')
#        plt.plot(wave, dark, color = 'green')
        plt.plot(l1b['wavelength line-by-line'][:], l1b['radiance line-by-line'][ialt,iact,:]/l1b['solar irradiance line-by-line'][:], color = 'blue')
        plt.plot(wave, dark/sun, color = 'green')
        plt.title('level 1B radiance')
        plt.ylabel('I [ph/(nm m$^2$ sr s)]')    
        plt.xlabel('wavelength [nm]')    
        plt.show()       

        plt.savefig('./plots/radiance.png')
