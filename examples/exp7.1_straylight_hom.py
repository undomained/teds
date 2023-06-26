
# exp7.1: straylight analysis up to level 0 for a homogenous light source. We conpare L1A for 
# 1. no stray light
# 2. single kernel simulations
# 3. full kernel simulations

import context
from end_to_end.GM.gm import geometry_module
from end_to_end.SGM.sgm import scene_generation_module
from end_to_end.IM.create_im_configuration_file import im_configuration
from end_to_end.L0L1B.create_l01b_configuration_file import l0l1b_configuration
from end_to_end.L0L1B.pixel_mask import determine_pixel_mask
from end_to_end.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor
from end_to_end.L1L2.l1bl2 import level1b_to_level2_processor

import yaml
import sys
import shutil
import subprocess
import netCDF4 as nc
import numpy as np

if __name__ == "__main__":
    # end to end global config file
    
    import paths
    import constants
    
    global_config = yaml.safe_load(open("./exp7.1_straylight_hom.yaml"))

    # ======= geometry module ======================================
    # choose baseline GM config
    if(global_config['setup']['flag_gm']):

        shutil.copyfile(paths.project+ paths.GM_module+'gm_config_baseline.yaml', \
                        paths.project+ paths.GM_module+'gm_config.yaml',)

        gm_config = yaml.safe_load(open(paths.project+paths.GM_module+'gm_config.yaml'))
        geometry_module(paths, global_config, gm_config)

    # ======= scene generator module ===============================
    # choose baseline SGM config
    if(global_config['setup']['flag_sgm']):
        shutil.copyfile(paths.project+paths.SGM_module+'sgm_config_baseline.yaml', \
                        paths.project+paths.SGM_module+'sgm_config.yaml',)
        sgm_config = yaml.safe_load(open(paths.project+paths.SGM_module+ "sgm_config.yaml"))
        #scene_generation_module(paths, global_config, sgm_config)

    #overwrite sgm radiance file wit a constant radiance value = max(radiance)
    path_sgm_data = paths.project+paths.data_interface + paths.interface_sgm                                         
    file_rad = sgm_config['sgm_filename']['filename_rad_output']  + '_' + global_config['profile']+'_'+global_config['run_id']+'.nc'
    sgm_set = nc.Dataset(path_sgm_data+file_rad,'r+')
    sgm_set.variables['radiance'][:]
    sgm_set['radiance'][:]=np.max(sgm_set['radiance'])
    sgm_set.close()     

    # ======= The instrument model =================================        
    if(global_config['setup']['flag_im']):

        #with constant stray light kernels
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 1
        im_config['filename']['ckd'] =    'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel'
        im_config['filename']['level1a'] ='Tango_Carbon_l0_1kernel_'                   #level 1a input
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #with constant stray light kernel
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 1
        im_config['filename']['ckd'] =    'OWL640S_low-gain_radiation_ckd_dose0.0_21kernel'
        im_config['filename']['level1a'] ='Tango_Carbon_l0_21kernel_'                   #level 1a input
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

        #with varying stray light kenrels
        im_config = yaml.safe_load(open(paths.project+paths.IM_module+"im_config.yaml"))
        im_config['noise']['switch'] = 0
        im_config['settings']['bin_id']              = 1
        im_config['settings']['sw_stray']            = 0
        im_config['filename']['ckd'] =    'OWL640S_low-gain_radiation_ckd_dose0.0_1kernel'
        im_config['filename']['level1a'] ='Tango_Carbon_l0_without_strayl_'                   #level 1a input
        im_configuration(paths,global_config, im_config)     
        cmd_str = paths.project+paths.IM_module + 'tango_ckd_model/build/ckdmodel im_config.cfg'
        subprocess.run(cmd_str, shell=True, cwd=paths.project+paths.IM_module)

    print('Experiment 7.1 sucessfully performed.')
    sys.exit()
