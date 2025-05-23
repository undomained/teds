"""
     wrapper script for Tango's end-to-end simulator
     full chain implementation using Sentinel 2 + microHH scene
     This source code is licensed under the 3-clause BSD license found in
     the LICENSE file in the root directory of this project.
"""
# define  path to search for module
import sys
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp5.0_rad_offset/")
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/")
    
# import E2ES modules 
from teds.GM.gm import geometry_module
from teds.GM.create_gm_yaml_file import create_gm_config_file
from teds.SGM.sgm import scene_generation_module
from teds.SGM.create_sgm_yaml_file import create_sgm_config_file
from teds.IM.create_im_configuration_file import im_configuration
from teds.L1AL1B.create_l1a1b_configuration_file import l1al1b_configuration
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
locations.gm['output']        = path_interface + 'gm/Tango_Carbon_gm_exp5.0.nc'

locations.__setattr__('sgm', {})
locations.sgm['gm_input']     = path_interface + 'gm/Tango_Carbon_gm_exp5.0.nc'
locations.sgm['S2_dump']      = path_tmp + 'Tango_Carbon_S2_exp5.0.npy'
locations.sgm['afgl_input']   = path_afgl + 'prof.AFGL.US.std'
locations.sgm['meteo_dump']   = path_tmp + 'Tango_Carbon_meteo_exp5.0.pkl'
locations.sgm['xsec_dump']    = path_tmp + 'Tango_Carbon_xsec_exp5.0.pkl'
locations.sgm['sun_reference']= path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.sgm['rad_output']   = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp5.0.nc'
locations.sgm['geo_output']   = path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_exp5.0.nc'
locations.sgm['hapi_path']    = path_harpi

locations.__setattr__('im', {})
locations.im['ckd_input']     = path_interface + 'ckd/OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
locations.im['rad_input']     = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp5.0.nc'
locations.im['binning_table'] = path_interface + 'ckd/binning_table.nc'
locations.im['ckd_expo_time'] = path_interface + 'ckd/optimal_exposure.txt'
locations.im['ckd_coadding']  = path_interface + 'ckd/optimal_coadding.txt'
locations.im['IM_path']       = path_IM
locations.im['output']        = path_interface + 'level1a/Tango_Carbon_l1a_exp5.0.nc'
locations.im['gm_input']      = path_interface + 'gm/Tango_Carbon_gm_exp5.0.nc'

locations.__setattr__('l1al1b', {})
locations.l1al1b['ckd_input'] = path_interface + 'ckd/OWL640S_low-gain_radiation_ckd_dose0.0_21kernel.nc'
locations.l1al1b['binning_table'] = path_interface +'ckd/binning_table.nc'
locations.l1al1b['l1a_input'] = path_interface + 'level1a/Tango_Carbon_l1a_exp5.0.nc'
locations.l1al1b['l1b_output']= path_interface + 'level1b/Tango_Carbon_l1b_exp5.0.nc'
locations.l1al1b['gm_input']  = path_interface + 'gm/Tango_Carbon_gm_exp5.0.nc'
locations.l1al1b['L1AL1B_path']   = path_L1AL1B

locations.__setattr__('l1bl2', {})
locations.l1bl2['l1b_input']  = path_interface + 'level1b/Tango_Carbon_l1b_exp5.0.nc'
locations.l1bl2['pixel_mask'] = ''  #needs to be specified only if mask =True
locations.l1bl2['afgl_input'] = path_afgl      + 'prof.AFGL.US.std'
locations.l1bl2['sun_reference'] = path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.l1bl2['l2_output'] = path_interface  + 'level2/Tango_Carbon_l2_exp5.0.nc'
locations.l1bl2['xsec_dump'] = path_tmp        + 'Tango_Carbon_xsec_l2_exp5.0.pkl'
locations.l1bl2['l2_diags']  = ''
locations.l1bl2['hapi_path'] =  path_harpi
locations.l1bl2['sgm_input'] = path_interface  + 'sgm/Tango_Carbon_sgm_atmosphere_exp5.0.nc'

#scene specification for profile single_pixel and swath
scene = {}
scene['scene_spec'] = {}
nact = 100
scene['scene_spec']['sza']    = (np.zeros(nact)).tolist()
scene['scene_spec']['saa']    = (np.zeros(nact)).tolist()
scene['scene_spec']['vza']    = (np.zeros(nact)).tolist()
scene['scene_spec']['vaa']    = (np.zeros(nact)).tolist()

zero_angles = (np.zeros(nact)).tolist()
zero_albedos = (np.zeros(nact)).tolist()
# =============================================================================

profile= 'single_swath'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']        = True
settings['sgm']       = True
settings['im']        = True
settings['l1al1b']    = True
settings['l1bl2']     = True
settings['save_yaml'] = False
# ====================main part ===============================================
if __name__ == "__main__":

    albedos = np.arange(41)*0.01 + 0.01
    angles  = [70.]
    
    for alb in albedos:

        print('======================================================')
        print('albedo: ',alb)
        print('======================================================')
        scene['scene_spec']['albedo'] = [x+alb for x in zero_albedos] 

        for sza in angles:

            scene['scene_spec']['sza']= [x+sza for x in zero_angles]

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

            # ======= The instrument model ========================================

            if(settings['im']):

            # with constant stray light kernel
                im_config= yaml.safe_load(open('./settings/im_config.yaml'))
                im_config['noise']['switch']= 1
                im_config['settings']['sw_stray']= 0
                im_configuration(locations.im, im_config)
                cmd_str= '../../teds/IM/tango_ckd_model/build/ckdmodel im_config.cfg'
                subprocess.run(cmd_str, shell=True, cwd=path_IM)

            # ======= The L0L1 pocessor ===========================================
            if(settings['l1al1b']):

                print('L0L1B processor')
                l1al1b_config= yaml.safe_load(open('./settings/l1al1b_config.yaml'))
                locations.l1al1b['l1b_output'] = path_interface + 'level1b/Tango_Carbon_l1b_exp5.0_sza'+"{:.1f}".format(sza)+'_alb'+"{:.3f}".format(alb)+'.nc'
                l1al1b_config['settings']['van_cittert_steps']= 0
                l1al1b_configuration(locations.l1al1b, l1al1b_config)
                cmd_str= '../../teds/L1AL1B/tango_l1b/build/tango_l1b l1al1b_config.cfg'
                subprocess.run(cmd_str, shell=True, cwd=path_L1AL1B)

            # ======= L1 to L2 processor ======================================

            if(settings['l1bl2']):
 
                # choose baseline L1BL2 config
                locations.l1bl2['l2_output'] = path_interface  + 'level2/Tango_Carbon_l2_exp5.0_sza'+"{:.1f}".format(sza)+'_alb'+"{:.3f}".format(alb)+'.nc'
                locations.l1bl2['l1b_input'] = path_interface + 'level1b/Tango_Carbon_l1b_exp5.0_sza'+"{:.1f}".format(sza)+'_alb'+"{:.3f}".format(alb)+'.nc'
                l1bl2_config= yaml.safe_load(open('./settings/l1bl2_config_baseline.yaml'))
                l1bl2_config['retrieval_init']['sw_pixel_mask']= False
                l1bl2_config = {**locations.l1bl2, **l1bl2_config}
                level1b_to_level2_processor(l1bl2_config)

    sys.exit()