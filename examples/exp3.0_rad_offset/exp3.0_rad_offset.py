"""
     wrapper script for Tango's end-to-end simulator
     full chain implementation using Sentinel 2 + microHH scene
     This source code is licensed under the 3-clause BSD license found in
     the LICENSE file in the root directory of this project.
"""
# define  path to search for module
import sys
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp3.0_rad_offset/")
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/")
    
# import E2ES modules 
import paths
from teds.GM.gm import geometry_module
from teds.GM.create_gm_yaml_file import create_gm_config_file
from teds.SGM.sgm import scene_generation_module
from teds.SGM.create_sgm_yaml_file import create_sgm_config_file
from teds.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor
from teds.SIML1B.radiance_offset import radiance_offset
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

locations = Emptyclass()
locations.__setattr__('gm', {})
locations.gm['output']        = path_interface + 'gm/Tango_Carbon_gm_exp3.0.nc'

locations.__setattr__('sgm', {})
locations.sgm['gm_input']     = path_interface + 'gm/Tango_Carbon_gm_exp3.0.nc'
locations.sgm['S2_dump']      = path_tmp + 'Tango_Carbon_S2_exp3.0.npy'
locations.sgm['afgl_input']   = path_afgl + 'prof.AFGL.US.std'
locations.sgm['meteo_dump']   = path_tmp + 'Tango_Carbon_meteo_exp3.0.pkl'
locations.sgm['xsec_dump']    = path_tmp + 'Tango_Carbon_xsec_exp3.0.pkl'
locations.sgm['sun_reference']= path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.sgm['rad_output']   = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp3.0.nc'
locations.sgm['geo_output']   = path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_exp3.0.nc'
locations.sgm['hapi_path']    = paths.project  + paths.data_harpi

locations.__setattr__('siml1b', {})
locations.siml1b['sgm_input']  = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp3.0.nc'
locations.siml1b['gm_input']   = path_interface + 'gm/Tango_Carbon_gm_exp3.0.nc'
locations.siml1b['l1b_output'] = path_interface + 'level1b/Tango_Carbon_l1b_exp3.0.nc'

locations.__setattr__('l1bl2', {})
locations.l1bl2['l1b_input']  = path_interface + 'level1b/Tango_Carbon_l1b_exp3.0.nc'
locations.l1bl2['pixel_mask'] = ''  #needs to be specified only if mask =True
locations.l1bl2['afgl_input'] = path_afgl      + 'prof.AFGL.US.std'
locations.l1bl2['sun_reference'] = path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.l1bl2['l2_output'] = path_interface  + 'level2/Tango_Carbon_l2_exp3.0.nc'
locations.l1bl2['xsec_dump'] = path_tmp        + 'Tango_Carbon_xsec_l2_exp3.0.pkl'
locations.l1bl2['l2_diags']  = ''
locations.l1bl2['hapi_path'] =  path_harpi
locations.l1bl2['sgm_input'] = path_interface  + 'sgm/Tango_Carbon_sgm_atmosphere_exp3.0.nc'

#scene specification for profile single_pixel and swath
scene = {}
scene['scene_spec'] = {}
scene['scene_spec']['sza']    = [70., 60, 50, 40, 30, 20, 10, 0] 
scene['scene_spec']['saa']    = [0.,  0., 0., 0., 0., 0., 0., 0] 
scene['scene_spec']['vza']    = [0.,  0., 0., 0., 0., 0., 0., 0] 
scene['scene_spec']['vaa']    = [0.,  0., 0., 0., 0., 0., 0., 0] 
scene['scene_spec']['albedo'] = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15] 


# =============================================================================

profile= 'individual_spectra'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']        = True
settings['sgm']       = True
settings['siml1b']    = False
settings['l1bl2']     = False
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

    if(settings['l1bl2']):

        rad_offsets = np.arange(11)*0.005
        for rad_offset in rad_offsets:
    
            filename_in  = path_interface + 'level1b/Tango_Carbon_l1b_exp3.0.nc'
            filename_out = path_interface + 'level1b/Tango_Carbon_l1b_exp3.0_rad_offset.nc'     
            radiance_offset(filename_in, filename_out, rad_offset)
    
            # ======= L1 to L2 processor ===================================
            # choose baseline L1BL2 config
            locations.l1bl2['l1b_input']  = filename_out
            locations.l1bl2['l2_output']  = filename_out
            locations.l1bl2['l2_output'] = path_interface + 'level2/Tango_Carbon_l2_exp3.0_rad_offset'+ \
                "{:.3f}".format(rad_offset)+'.nc'
            l1bl2_config= yaml.safe_load(open('./settings/l1bl2_config_baseline.yaml'))
            l1bl2_config['retrieval_init']['sw_pixel_mask']= False
            l1bl2_config['isrf_settings']['type']= 'generalized_normal'  # type of ISRF, currently only Gaussian or generalized_normal
            l1bl2_config['isrf_settings']['fwhm']=  0.45  # fwhm  [nm]
            l1bl2_config['isrf_settings']['bcoeff']=  0.45  # fwhm  [nm]

            l1bl2_config = {**locations.l1bl2, **l1bl2_config}
            level1b_to_level2_processor(l1bl2_config)

