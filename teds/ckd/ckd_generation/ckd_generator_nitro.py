"""
TANGO NITRO CKD generator of the TANGO nitro instrument.

Individuals CKD and their calculation can be found in /nitro_ckds/ folder
"""

import os
import sys
import yaml
from glob import glob
from teds.ckd.ckd_generation.nc_container import *

def load_config(config_path):
    with open(config_path) as file:
        return yaml.safe_load(file)

    
def flatten_dict(d, parent_key='', sep='/'):
    """
    Flatten a nested dictionary with a specified separator.
    """
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def generate_ckds(ncc, cfg):
    varpaths = glob(f"{cfg['paths']['dir_nitro']}/**/*.py", recursive=True)
    for varpath in varpaths:
        if 'skip' in varpath:
            continue
        mod = import_mod(varpath)
        mod.generate(ncc)

def main(cfg): 
    os.makedirs(cfg['paths']['dir_nitro'], exist_ok=True)
    ncc = NC_container(cfg)
    generate_ckds(ncc, cfg)
    ncc.close()
    
if __name__ == "__main__":
    # Set the current working directory to the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    #os.chdir(script_dir)
    if len(sys.argv) < 2:
        config_path = '{}/cfg/nitro/full_config.yaml'.format(script_dir)
        print(f"No argument given, using {config_path}")
    else:
        config_path = sys.argv[1]

    full_config = load_config(config_path)
    cfg = full_config['ckd']
    main(cfg)