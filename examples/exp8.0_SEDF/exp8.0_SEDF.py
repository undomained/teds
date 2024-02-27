"""
     wrapper script for Tango's end-to-end simulator
     full chain implementation using Sentinel 2 + microHH scene
     This source code is licensed under the 3-clause BSD license found in
     the LICENSE file in the root directory of this project.
"""
# define  path to search for module
import sys
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/examples/exp8.0_detector_mapping/")
sys.path.insert(1, "/home/jochen/TANGO_E2ES/EndtoEndProject/end_to_end/")
    
# import E2ES modules 
from teds.GM.gm import geometry_module
from teds.GM.create_gm_yaml_file import create_gm_config_file
from teds.SGM.sgm import scene_generation_module_new
from teds.SGM.create_sgm_yaml_file import create_sgm_config_file
from teds.IM.create_im_configuration_file import im_configuration
from teds.L1AL1B.create_l1a1b_configuration_file import l1al1b_configuration
from teds.L1L2.l1bl2 import level1b_to_level2_processor
from teds.L1L2.create_l2_yaml_file import create_l1bl2_config_file
#import other modules
import yaml
import subprocess
import numpy as np
import os
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
locations.gm['output']        = path_interface + 'gm/Tango_Carbon_gm_exp8.0.nc'

locations.__setattr__('sgm', {})
locations.sgm['gm_input']     = path_interface + 'gm/Tango_Carbon_gm_exp8.0.nc'
locations.sgm['S2_dump']      = path_tmp + 'Tango_Carbon_S2_exp8.0.npy'
locations.sgm['afgl_input']   = path_afgl + 'prof.AFGL.US.std'
locations.sgm['meteo_dump']   = path_tmp + 'Tango_Carbon_meteo_exp8.0.pkl'
locations.sgm['xsec_dump']    = path_tmp + 'Tango_Carbon_xsec_exp8.0.pkl'
locations.sgm['sun_reference']= path_sol_spec + 'hybrid_reference_spectrum_c2021-03-04_with_unc.nc'
locations.sgm['rad_output']   = path_interface + 'sgm/Tango_Carbon_sgm_radiance_exp8.0.nc'
locations.sgm['geo_output']   = path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_exp8.0.nc'
locations.sgm['hapi_path']    = path_harpi

#scene specification for profile single_pixel and swath
scene = {}
scene['scene_spec'] = {}
nact = 100
scene['scene_spec']['sza']    = (np.zeros(nact)+50.).tolist()
scene['scene_spec']['saa']    = (np.zeros(nact)).tolist()
scene['scene_spec']['vza']    = (np.zeros(nact)).tolist()
scene['scene_spec']['vaa']    = (np.zeros(nact)).tolist()
alb = np.zeros(nact) + 0.10
alb[::2] = 0.60
scene['scene_spec']['albedo'] = alb.tolist()

# =============================================================================

profile= 'orbit'   #needed to initialize gm and sgm consistently

settings= {}
settings['gm']        = False
settings['sgm']       = True

microHH_sim = ['src1_20180523_0900',  'src1_20180523_1015',  'src1_20180523_1130',  
               'src1_20180523_1245',  'src1_20180523_0915',  'src1_20180523_1030',  
               'src1_20180523_1145',  'src1_20180523_1300',  'src1_20180523_0930',
               'src1_20180523_1045',  'src1_20180523_1200',  'src1_20180523_0945',
               'src1_20180523_1100',  'src1_20180523_1215',  'src1_20180523_1000', 
               'src1_20180523_1115',  'src1_20180523_1230']
# ====================main part ===============================================
if __name__ == "__main__":

    # ======= geometry module =================================================

    if(settings['gm']):

        config= yaml.safe_load(open('./settings/gm_config_baseline.yaml'))
        gm_config = {**locations.gm, **config, **scene}
        gm_config['profile'] = profile
        geometry_module(gm_config)

    # ======= scene generator module ==========================================

    if(settings['sgm']):

        sgm_config= yaml.safe_load(open('./settings/sgm_config_baseline.yaml'))
        sgm_config = {**locations.sgm, **sgm_config, **scene}
        sgm_config['profile'] = profile
        
        for ic, sim in enumerate(microHH_sim):
            sgm_config['meteo']['filesuffix']=sim  # sourcename_yyyymmdd_hhmm
            sgm_config['geo_output_raw']= path_interface + 'sgm/Tango_Carbon_sgm_atmosphere_raw_exp8.0_'+ sim +'.nc'
            scene_generation_module_new(sgm_config, sw_raw_geo_data_only= False)

            sys.exit()
