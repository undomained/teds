import numpy as np
import sys
import os
import yaml

lib_path = os.path.abspath(os.path.join('../lib'))
if lib_path not in sys.path:
    sys.path.append(lib_path)
    
import libSGM 
import constants

path_gm = '../GM/gm_data/'
file_gm = 'Tango_Carbon_gm_20230213'

gm_data = libSGM.get_gm_data(path_gm,file_gm)

libSGM.get_sentinel2_albedo(gm_data)