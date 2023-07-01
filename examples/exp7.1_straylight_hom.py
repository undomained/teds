
# exp7.1: straylight analysis up to level 0 for a homogenous light source. We conpare L1A for 
# 1. no stray light
# 2. single kernel simulations
# 3. full kernel simulations

# define  path to search for module
import sys
path = "../"
if(str(path) not in sys.path):
    sys.path.append(path)

path = "../end_to_end/lib/"
if(str(path) not in sys.path):
    sys.path.append(path)

from end_to_end.lib import paths
from end_to_end.GM.gm import geometry_module
from end_to_end.SGM.sgm import scene_generation_module
from end_to_end.IM.create_im_configuration_file import im_configuration

import yaml
import sys
import subprocess
import netCDF4 as nc
import numpy as np

class Emptyclass:
    pass

#run id
run_id = 'exp7.1'

# paths and file names
locations = Emptyclass()
locations.__setattr__('gm', {})
locations.gm['output'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'

locations.__setattr__('sgm', {})
locations.sgm['gm_input'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'
locations.sgm['S2_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_S2_' + run_id + '.npy'
locations.sgm['afgl_input'] = paths.project + \
    paths.data_afgl + 'prof.AFGL.US.std'
locations.sgm['microHH_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_microHH_' + run_id + '.pkl'
locations.sgm['microHH_data_path'] = paths.project + \
    paths.data_microHH + 'Jaenschwalde_simulation1/'
locations.sgm['xsec_dump'] = paths.project + \
    paths.data_tmp + 'Tango_Carbon_xsec_' + run_id + '.pkl'
locations.sgm['sun_reference'] = paths.project + paths.data_sol_spec + \
    'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.sgm['rad_output'] = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_radiance_' + run_id + '.nc'
locations.sgm['geo_output'] = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_atmosphere_' + run_id + '.nc'
locations.sgm['hapi_path'] =  paths.project + paths.data_harpi

locations.__setattr__('im', {})
locations.im['ckd_input'] = paths.project + paths.data_interface + \
    paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
locations.im['rad_input'] = paths.project + paths.data_interface + \
    paths.interface_sgm + 'Tango_Carbon_sgm_radiance_' + run_id + '.nc'
locations.im['binning_table'] = paths.project + \
    paths.data_interface + paths.interface_ckd + 'binning_table.nc'
locations.im['ckd_expo_time'] = paths.project + \
    paths.data_interface + paths.interface_ckd + 'optimal_exposure.txt'
locations.im['ckd_coadding'] = paths.project + \
    paths.data_interface + paths.interface_ckd + 'optimal_coadding.txt'
locations.im['IM_path'] = paths.project + paths.IM_module
locations.im['output'] = paths.project + paths.data_interface + \
    paths.interface_l1a + 'Tango_Carbon_l1a_' + run_id + '.nc'
locations.im['gm_input'] = paths.project + paths.data_interface + \
    paths.interface_gm + 'Tango_Carbon_gm_' + run_id + '.nc'

#scene specification for profile single_pixel and swath
scene_spec = {}
scene_spec['numb_atm_scenes']= 2
scene_spec['scene_trans_index'] = [0,50, 100]
scene_spec['sza'] = [70., 70.]
scene_spec['saa']= [0., 0.]
scene_spec['vza']= [0., 0.]
scene_spec['vaa']= [0., 0.]
scene_spec['albedo']= [0.15, 0.70]

# =============================================================================
#
#  Current version of the Tango E2E simulator can handle the following profiles
#   1.  individual_spectra: In this case only one single spectrum is processed.
#       The IM and L1 processor cannot be used in this case but are replaced by
#       the simplified moduel siml1. scene_spec required
#   2.  single_swath: This profile acounts for a single swath and a homogenous
#       reference scene. It is meant to test the IM and the L1 processor
#       scene_spec required
#   3.  S2_microHH: This simulates a full data granuale of 30 x 30 km2 using
#       flexibale geometry information and Sentinel 2 albedo data and
#       microHH CO2 plume simulations. scene_spec not required
# =============================================================================

profile= 'single_swath'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']    = True
settings['sgm']   = True
settings['im']    = True
           
if __name__ == "__main__":
    # end to end global config file

    # ======= geometry module ======================================
    # choose baseline GM config
    if(settings['gm']):

        config= yaml.safe_load(open(paths.project+paths.GM_module+'gm_config_baseline.yaml'))
        gm_config = {**locations.gm, **config, **scene_spec}
        gm_config['profile'] = profile
        
        geometry_module(gm_config)


    # ======= scene generator module ===============================
    # choose baseline SGM config
    sgm_config= yaml.safe_load(open(paths.project+paths.SGM_module + "sgm_config_baseline.yaml"))
    sgm_config = {**locations.sgm, **sgm_config, **scene_spec}
    sgm_config['profile'] = profile
    if(settings['sgm']):

        scene_generation_module(sgm_config)

    #overwrite sgm radiance file wit a constant radiance value = max(radiance)
    path_sgm_data = paths.project+paths.data_interface + paths.interface_sgm                                         
    file_rad = sgm_config['rad_output']
    sgm_set = nc.Dataset(file_rad,'r+')
    sgm_set.variables['radiance'][:]
    sgm_set['radiance'][:]=np.max(sgm_set['radiance'])
    sgm_set.close()     

    # ======= The instrument model =================================        
    if(settings['im']):

        #with constant stray light kernels
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 1
        locations.im['ckd_input'] = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel.nc'
        locations.im['output'] = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_1kernel_' + run_id + '.nc'
        im_configuration(locations.im, im_config)
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #with constant stray light kernel
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 1
        locations.im['ckd_input'] = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
        locations.im['output'] = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_21kernel_' + run_id + '.nc'
        im_configuration(locations.im, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #with varying stray light kenrels
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 0
        
        locations.im['ckd_input'] = paths.project + paths.data_interface + \
            paths.interface_ckd + 'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel.nc'
        locations.im['output'] = paths.project + paths.data_interface + \
            paths.interface_l1a + 'Tango_Carbon_l1a_without_strayl_' + run_id + '.nc'
        im_configuration(locations.im, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

    print('Experiment 7.1 sucessfully performed.')
    sys.exit()
