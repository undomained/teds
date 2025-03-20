"""
     wrapper script for Tango's end-to-end simulator
     full chain implementation using Sentinel 2 + microHH scene
     This source code is licensed under the 3-clause BSD license found in
     the LICENSE file in the root directory of this project.
"""
# define  path to search for module
import sys
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp2.0_RTM/")
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/")
    
# import E2ES modules 
import paths
from teds.GM.gm import geometry_module
from teds.GM.create_gm_yaml_file import create_gm_config_file
from teds.SGM.sgm import scene_generation_module
from teds.SGM.create_sgm_yaml_file import create_sgm_config_file
from teds.IM.create_im_configuration_file import im_configuration
from teds.L1AL1B.create_l1a1b_configuration_file import l1al1b_configuration
from teds.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor
from teds.L1L2.l1bl2 import level1b_to_level2_processor
from teds.L1L2.create_l2_yaml_file import create_l1bl2_config_file

#import other modules
import yaml
import subprocess
import numpy as np

# ====================configuration part ======================================

class Emptyclass:
    pass

# paths and file names
path_project = '/home/jochen/TANGO_E2ES/EndtoEndProject/'

path_interface = path_project  + 'data/interface_data/'
path_tmp       = path_project  + 'data/tmp/'
path_afgl      = path_project  + 'data/AFGL/'
path_sol_spec  = path_project  + 'data/solar_spectra/'
path_harpi     = path_project  + 'data/harpi/'
path_IM        = path_project  + 'end_to_end/teds/IM/'
path_L1AL1B    = path_project  + 'end_to_end/teds/L1AL1B/'

locations = Emptyclass()
locations.__setattr__('gm', {})
locations.gm['output']        = path_interface + 'gm/Tango_Carbon_gm_exp4.1.nc'

locations.__setattr__('sgm', {})
locations.sgm['geometry']     = path_interface + 'gm/Tango_Carbon_gm_exp4.1.nc'
locations.sgm['s2_albedo']      = path_tmp + 'Tango_Carbon_S2_exp4.1.npy'
locations.sgm['afgl']   = path_afgl + 'prof.AFGL.US.std'
locations.sgm['meteo']   = path_tmp + 'Tango_Carbon_meteo_exp4.1.pkl'
locations.sgm['xsec_dump']    = path_tmp + 'Tango_Carbon_xsec_exp4.1.pkl'
locations.sgm['sun_reference']= path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.sgm['radiance']   = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp4.1.nc'
locations.sgm['geometry']   = path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_exp4.1.nc'
locations.sgm['hapi']    = paths.project  + paths.data_harpi

locations.__setattr__('im', {})
locations.im['ckd']     = path_interface + 'ckd/OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
locations.im['sgm']     = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp4.1.nc'
locations.im['binning_table'] = path_interface + 'ckd/binning_table.nc'
locations.im['ckd_expo_time'] = path_interface + 'ckd/optimal_exposure.txt'
locations.im['ckd_coadding']  = path_interface + 'ckd/optimal_coadding.txt'
locations.im['IM_path']       = path_IM
locations.im['l1a']        = path_interface + 'level1a/Tango_Carbon_l1a_exp4.1.nc'
locations.im['geometry']      = path_interface + 'gm/Tango_Carbon_gm_exp4.1.nc'

locations.__setattr__('l1al1b', {})
locations.l1al1b['ckd'] = path_interface + 'ckd/OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
locations.l1al1b['binning_table'] = path_interface +'ckd/binning_table.nc'
locations.l1al1b['l1a'] = path_interface + 'level1a/Tango_Carbon_l1a_exp4.1.nc'
locations.l1al1b['l1b']= path_interface + 'level1b/Tango_Carbon_l1b_exp4.1.nc'
locations.l1al1b['geometry']  = path_interface + 'gm/Tango_Carbon_gm_exp4.1.nc'
locations.l1al1b['L1AL1B_path']   = path_L1AL1B

locations.__setattr__('siml1b', {})
locations.siml1b['sgm_input']  = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp4.1.nc'
locations.siml1b['gm_input']   = path_interface + 'gm/Tango_Carbon_gm_exp4.1.nc'
locations.siml1b['l1b_output'] = path_interface + 'level1b/Tango_Carbon_l1b_exp4.1_siml1b.nc'

locations.__setattr__('l1bl2', {})
locations.l1bl2['l1b_input']  = path_interface + 'level1b/Tango_Carbon_l1b_exp4.1.nc'
locations.l1bl2['pixel_mask'] = ''  #needs to be specified only if mask =True
locations.l1bl2['afgl_input'] = path_afgl      + 'prof.AFGL.US.std'
locations.l1bl2['sun_reference'] = path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.l1bl2['l2_output'] = path_interface  + 'level2/Tango_Carbon_l2_exp4.1.nc'
locations.l1bl2['xsec_dump'] = path_tmp        + 'Tango_Carbon_xsec_l2_exp4.1.pkl'
locations.l1bl2['l2_diags']  = ''
locations.l1bl2['hapi_path'] =  path_harpi
locations.l1bl2['sgm_input'] = path_interface  + 'sgm/Tango_Carbon_sgm_atmosphere_exp4.1.nc'

#scene specification for profile single_pixel and swath
#scene specification for profile single_pixel and swath
scene = {}
scene['scene_spec'] = {}
scene['scene_spec']['numb_atm']= 1
scene['scene_spec']['scene_trans_index'] = [0,100]
scene['scene_spec']['sza'] = [70.]
scene['scene_spec']['saa']= [0.]
scene['scene_spec']['vza']= [0.]
scene['scene_spec']['vaa']= [0.]
scene['scene_spec']['albedo']= [0.15]

# =============================================================================
profile = 'single_swath'  

settings= {}
settings['gm']        = False
settings['sgm']       = False
settings['im']        = True
settings['l1al1b']    = True
settings['l1bl2']     = True
settings['siml1b']    = False
settings['save_yaml'] = False
# ====================main part ===============================================
if __name__ == "__main__":

    # ======= geometry module =================================================

    if(settings['gm']):

        config= yaml.safe_load(open('./settings/gm_config_baseline.yaml'))
        gm_config = {**locations.gm, **config, **scene}
        gm_config['profile'] = profile
        if(settings['save_yaml']):
            gm_yaml = './save_yaml/gm_config_exp2.0.yaml'
            create_gm_config_file(gm_yaml, gm_config)
        geometry_module(gm_config)

    # ======= scene generator module ==========================================

    if(settings['sgm']):

        sgm_config= yaml.safe_load(open('./settings/sgm_config_baseline.yaml'))
        sgm_config = {**locations.sgm, **sgm_config, **scene}
        sgm_config['profile'] = profile
        if(settings['save_yaml']):
            sgm_yaml = './save_yaml/sgm_config_exp2.0.yaml'
            create_sgm_config_file(sgm_yaml, sgm_config)
        scene_generation_module(sgm_config)

    dimt = 21
    reference_temperature = 10.
    temperature_range = (np.arange(dimt)*0.01 -0.1) + reference_temperature
    scale = 10.
    for target_temperature in temperature_range:
     
        print('===============================================================')
        print('temperature: ', f"{target_temperature: 04.2f}")
        print('===============================================================')
        # ======= The instrument model =======================================

        if(settings['im']):

            # with constant stray light kernel
            im_config= yaml.safe_load(open('./settings/im_config.yaml'))
            im_config['noise']['switch']= 0
            im_config['settings']['bin_id']= 5
            im_config['settings']['sw_stray']= 0
            locations.im['ckd_input']  = path_interface + \
                'ckd/ckd_temperature_dependence/' +\
                'OWL640S_low_gain_21kernel_temp_'+f"{target_temperature:05.2f}" + \
                    '_scale_' + f"{scale:04.1f}" + ".nc"
            im_configuration(locations.im, im_config)
            cmd_str= '../../teds/IM/tango_ckd_model/build/ckdmodel im_config.cfg'
            subprocess.run(cmd_str, shell=True, cwd=path_IM)
        # ======= The L0L1 pocessor ==========================================
        if(settings['l1al1b']):

            l1al1b_config= yaml.safe_load(open('./settings/l1al1b_config.yaml'))
            l1al1b_config['settings']['van_cittert_steps']= 0
            locations.l1al1b['ckd_input'] = path_interface + \
                'ckd/ckd_temperature_dependence/' +\
                'OWL640S_low_gain_21kernel_temp_'+f"{reference_temperature:05.2f}" + \
                    '_scale_' + f"{scale:04.1f}" +  ".nc"
            l1al1b_configuration(locations.l1al1b, l1al1b_config)
            cmd_str= '../../teds/L1AL1B/tango_l1b/build/tango_l1b l1al1b_config.cfg'
            subprocess.run(cmd_str, shell=True, cwd=path_L1AL1B)

        if(settings['siml1b']):
        
            siml1b_config= yaml.safe_load(open('./settings/siml1b_config_baseline.yaml'))
            siml1b_config = {**locations.siml1b, **siml1b_config}
            simplified_instrument_model_and_l1b_processor(siml1b_config)

        if(settings['l1bl2']):
    
            # ======= L1 to L2 processor ===================================
            # choose baseline L1BL2 config
            locations.l1bl2['l2_output'] = path_interface + 'level2/Tango_Carbon_l2_exp4.1_'+\
                'temp_'+f"{target_temperature:05.2f}" + '_scale_' + f"{scale:04.1f}" + ".nc"

            l1bl2_config= yaml.safe_load(open('./settings/l1bl2_config_baseline.yaml'))
            l1bl2_config['retrieval_init']['sw_pixel_mask']= False
            l1bl2_config['isrf_settings']['type']= 'Gaussian'  # type of ISRF, currently only Gaussian or generalized_normal
            l1bl2_config['isrf_settings']['fwhm']=  0.45       # fwhm  [nm]
            l1bl2_config['isrf_settings']['bcoeff']=  0.45     # fwhm  [nm]

            l1bl2_config = {**locations.l1bl2, **l1bl2_config}
            level1b_to_level2_processor(l1bl2_config)

