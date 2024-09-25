"""
     wrapper script for Tango's end-to-end simulator
     full chain implementation using Sentinel 2 + microHH Ja4nschwalde scene
     This source code is licensed under the 3-clause BSD license found in
     the LICENSE file in the root directory of this project.
"""
# define  path to search for module
import os.path
import sys
#from typing_extensions import TypeVarTuple
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp1.0_full_chain/")
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/")
    
# import E2ES modules 

from teds.gm import geometry_module
from teds.sgm import geoscene_generation
from teds.sgm import Carbon_radiation_scene_generation
from teds.siml1b import simplified_instrument_model_and_l1b_processor
from teds.l1l2.l1bl2 import level1b_to_level2_processor, \
    level1b_to_level2_processor_RTorCH4
from teds.l1l2.sl2 import simplified_level2
from teds.l2l4.l2l4 import level2_to_level4_processor

#import other modules
import yaml
import subprocess

# ====================configuration part ======================================

# paths and file names

path_IM        = '/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/teds/IM/'
path_L1AL1B    = '/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/teds/L1AL1B/'

profile= 'orbit'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']           = False
settings['geosgm']       = False
settings['radsgm_Carbon']= False
settings['im']           = False
settings['l1al1b']       = False
settings['siml1b']       = False
settings['l1bl2']        = True
settings['sl2']          = False
settings['l2l4']         = False

# ====================main part ================================================
if __name__ == "__main__":

    # ======= geometry module ======================================

    if(settings['gm']):

        gm_config= yaml.safe_load(open('./settings/gm_config_baseline.yaml'))
        gm_config['profile'] = profile
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
        level1b_to_level2_processor_RTorCH4(l1bl2_config)

    # ======= simplifed L1 to L2 processor ===================================
    if(settings['sl2']):

        sl2_config = {}
        sl2_config['instrument'] = 'TANGO_Carbon'
        sl2_config['precision relative'] = 0.005
        sl2_config['precision constant'] = 0.00
        sl2_config['seed'] = 10

#       sl2_config = {**locations.sl2, **sl2_config}
#        simplified_level2(sl2_config)

    # ======= L2 to L4 processor ===================================
#    if(settings['l2l4']):
#        # choose baseline L1BL2 config
#
#        l2l4_config= yaml.safe_load(open('./settings/l2l4_config_baseline.yaml'))
#        l2l4_config = {**locations.l2l4, **l2l4_config}
#        level2_to_level4_processor(l2l4_config)
