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
from teds.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor

#import other modules
import yaml
import subprocess

# ====================configuration part ======================================

class Emptyclass:
    pass

# paths and file names

path_project = '/home/jochen/TANGO_E2ES/EndtoEndProject/'

path_interface = path_project  + 'data/interface_data/'
path_tmp       = path_project  + 'data/tmp/'
path_afgl      = path_project  + 'data/AFGL/'
path_sol_spec  = path_project  + 'data/solar_spectra/'

locations = Emptyclass()
locations.__setattr__('gm', {})
locations.gm['output']        = path_interface + 'gm/Tango_Carbon_gm_exp2.0.nc'

locations.__setattr__('sgm', {})
locations.sgm['geometry']     = path_interface + 'gm/Tango_Carbon_gm_exp2.0.nc'
locations.sgm['s2_albedo']      = path_tmp + 'Tango_Carbon_S2_exp2.0.npy'
locations.sgm['afgl']   = path_afgl + 'prof.AFGL.US.std'
locations.sgm['xsec_dump']    = path_tmp + 'Tango_Carbon_xsec_exp2.0.pkl'
locations.sgm['sun_reference']= path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.sgm['radiance']   = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp2.0.nc'
locations.sgm['atmosphere']   = path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_exp2.0.nc'
locations.sgm['hapi']    = paths.project  + paths.data_harpi

locations.__setattr__('siml1b', {})
locations.siml1b['sgm_input']  = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp2.0.nc'
locations.siml1b['gm_input']   = path_interface + 'gm/Tango_Carbon_gm_exp2.0.nc'
locations.siml1b['l1b_output'] = path_interface + 'level1b/Tango_Carbon_l1b_exp2.0.nc'

#scene specification for profile single_pixel and swath
scene = {}
scene['scene_spec'] = {}
scene['scene_spec']['sza']    = [70., 70., 0., 50.]        
scene['scene_spec']['saa']    = [0., 0., 0., 0. ] 
scene['scene_spec']['vza']    = [0., 0., 0., 0. ] 
scene['scene_spec']['vaa']    = [0., 0., 0., 0. ] 
scene['scene_spec']['albedo'] = [0.15, 0.10, 0.80, 0.25] 

# =============================================================================

profile= 'individual_spectra'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']        = True
settings['sgm']       = True
settings['siml1b']    = True
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

    # ======= The instrument model ============================================

    if(settings['siml1b']):
        
        siml1b_config= yaml.safe_load(open('./settings/siml1b_config_baseline.yaml'))
        siml1b_config = {**locations.siml1b, **siml1b_config}
        siml1b_config['profile'] = profile
        simplified_instrument_model_and_l1b_processor(siml1b_config)
    
