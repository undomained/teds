"""
     wrapper script for Tango's end-to-end simulator
     full chain implementation using Sentinel 2 + microHH Ja4nschwalde scene
     This source code is licensed under the 3-clause BSD license found in
     the LICENSE file in the root directory of this project.
"""
# define  path to search for module
import sys
from typing_extensions import TypeVarTuple
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp1.0_full_chain/")
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/")
    
# import E2ES modules 

from teds.GM.gm import geometry_module
from teds.GM.create_gm_yaml_file import create_gm_config_file
from teds.SGM.sgm_geoscene import geoscene_generation
from teds.SGM.sgm_Carbon_radscene import Carbon_radiation_scene_generation
from teds.SGM.create_sgm_yaml_file import create_sgm_config_file
from teds.IM.create_im_configuration_file import im_configuration
from teds.L1AL1B.create_l1a1b_configuration_file import l1al1b_configuration
from teds.SIML1B.siml1b import simplified_instrument_model_and_l1b_processor
from teds.L1L2.l1bl2 import level1b_to_level2_processor
from teds.L1L2.create_l2_yaml_file import create_l1bl2_config_file
from teds.L1L2.sl2 import simplified_level2
from teds.L2L4.l2l4 import level2_to_level4_processor

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
path_harpi     = path_project  + 'data/harpi/'
path_IM        = path_project  + 'end_to_end/teds/IM/'
path_L1AL1B    = path_project  + 'end_to_end/teds/L1AL1B/'

locations = Emptyclass()
#locations.__setattr__('gm', {})
#locations.gm['output']        = path_interface + 'gm/Tango_Carbon_gm_exp1.0.nc'

locations.__setattr__('l1bl2', {})
locations.l1bl2['l1b_input']  = path_interface + 'level1b/Tango_Carbon_l1b_exp1.0.nc'
locations.l1bl2['pixel_mask'] = ''  #needs to be specified only if mask =True
locations.l1bl2['afgl_input'] = path_afgl      + 'prof.AFGL.US.std'
locations.l1bl2['sun_reference'] = path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.l1bl2['l2_output'] = path_interface  + 'level2/Tango_Carbon_l2_exp1.0.nc'
locations.l1bl2['xsec_dump'] = path_tmp        + 'Tango_Carbon_xsec_l2_exp1.0.pkl'
locations.l1bl2['l2_diags']  = ''
locations.l1bl2['hapi_path'] =  path_harpi
locations.l1bl2['sgm_input'] = path_interface  + 'sgm/Tango_Carbon_sgm_atmosphere_exp1.0.nc'

locations.__setattr__('sl2', {})
locations.sl2['sgm_input']   = path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_exp1.0.nc'
locations.sl2['l2_output']   = path_interface + 'level2/Tango_Carbon_l2_exp1.0_impl.nc'

locations.__setattr__('l2l4', {})
locations.l2l4['sgm_input']  = path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_exp1.0.nc'
locations.l2l4['l2_input']   = path_interface + 'level2/Tango_Carbon_l2_exp1.0_simpl.nc'

profile= 'orbit'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']           = False
settings['geosgm']       = False
settings['radsgm_Carbon']= False
settings['im']           = False
settings['l1al1b']       = False
settings['siml1b']       = False
settings['l1bl2']        = True
settings['save_yaml']    = False
settings['sl2']          = False
settings['l2l4']         = False

# ====================main part ================================================
if __name__ == "__main__":

    # ======= geometry module ======================================

    if(settings['gm']):

        gm_config= yaml.safe_load(open('./settings/gm_config_baseline.yaml'))
        gm_config['profile'] = profile
        if(settings['save_yaml']):
            gm_yaml = './save_yaml/gm_config_exp1.0.yaml'
            create_gm_config_file(gm_yaml, gm_config)
        geometry_module(gm_config)

    # ======= scene generator module ===============================

    if(settings['geosgm']):
        sgm_config= yaml.safe_load(open('./settings/geosgm_config_baseline.yaml'))
        sgm_config['profile'] = profile
        geoscene_generation(sgm_config)

    if(settings['radsgm_Carbon']):

        sgm_config= yaml.safe_load(open('./settings/radsgm_config_baseline.yaml'))
        sgm_config['profile'] = profile
        Carbon_radiation_scene_generation(sgm_config)

    # ======= The instrument model =================================

    if(settings['im']):

        # with constant stray light kernel
        cmd_str= '../../teds/IM/tango_im/build/tango_im.x /home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp1.0_full_chain/settings/im_config.yaml'
        subprocess.run(cmd_str, shell=True, cwd=path_IM)

    # ======= The L0L1 pocessor ====================================
    if(settings['l1al1b']):

        cmd_str= '../../teds/L1AL1B/tango_l1b/build/tango_l1b.x /home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp1.0_full_chain/settings/l1al1b_config.yaml'
        subprocess.run(cmd_str, shell=True, cwd=path_L1AL1B)

    # ======= The simplified IM and L1B model ======================
    if(settings['siml1b']):

        siml1b_config= yaml.safe_load(open('./settings/siml1b_config_baseline.yaml'))        
        siml1b_config['profile'] = profile
        simplified_instrument_model_and_l1b_processor(siml1b_config)

    # ======= L1 to L2 processor ===================================
    if(settings['l1bl2']):
        # choose baseline L1BL2 config

        l1bl2_config= yaml.safe_load(open('./settings/l1bl2_config_baseline.yaml'))
        level1b_to_level2_processor(l1bl2_config)

    # ======= simplifed L1 to L2 processor ===================================
    if(settings['sl2']):

        sl2_config = {}
        sl2_config['instrument'] = 'TANGO_Carbon'
        sl2_config['precision relative'] = 0.005
        sl2_config['precision constant'] = 0.00
        sl2_config['seed'] = 10
        sl2_config = {**locations.sl2, **sl2_config}

        simplified_level2(sl2_config)

    # ======= L2 to L4 processor ===================================
    if(settings['l2l4']):
        # choose baseline L1BL2 config

        l2l4_config= yaml.safe_load(open('./settings/l2l4_config_baseline.yaml'))
        l2l4_config = {**locations.l2l4, **l2l4_config}
        level2_to_level4_processor(l2l4_config)
