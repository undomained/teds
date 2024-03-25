#!/bin/env python

"""
TANGO NITRO
CKD generator of the TANGO nitro instrument. Needs a config.yaml file
to run, can be used as an input argument in shell. The config file
includes all dimensions, global attributes and some single value 
parameters. 
CKD is generated using the folder structure in the dir_nitro folder as 
groups. Variables are created using .py files in the folder using the 
same name as the file. 
In the .py files is:
 - generate() method: returns an instance of the generator_class
   which is used to standardize the output and add the variable to the 
   ckd nc file. Needs config yaml as argument
   
 - (optional) import_data() method: when data is imported from external
   file which also contains a lot of other data. This method can then 
   export a file with only the useful data for this ckd. 
   Needs config yaml as argument
 - (optional) calculate() method if an elaborate or time consuming
   calculation is needed. It is useful if this exports a file with the 
   data which the gerenate() method can read to prevent recalculating
   all ckds every time. (Can be toggled in config yaml)
"""

from netCDF4 import Dataset
import yaml
import os
from datetime import datetime
from glob import glob
from teds.CKD.ckd_generation.generator_class import *

# set cwd to directory of this file
os.chdir(os.path.dirname(os.path.realpath(__file__))) 

if len(sys.argv) < 2:
    print('usage: ckd_generator.py [ckd.yaml]')
    if os.path.isfile('ckd_nitro.yaml'):
        conf = yaml.safe_load(open('./ckd_nitro.yaml'))
        print('no argument given, using ckd_nitro.yaml in this folder')
    else:
        exit(0)
else: # Load yaml in argument
    conf = yaml.safe_load(open(sys.argv[1]))

# Create folders if not already exists
for d in ['dir_nitro', 'dir_external', 'dir_input']:
    dpath = conf['paths'][d]
    if not os.path.exists(dpath):
        os.mkdir(dpath)

# If recalculation is toggled on, remove all input_files
# to reimport and recalculate data. 
input_files = glob(conf['paths']['dir_input']+'/*')
if conf['recalculate_ckds']:
    [os.remove(f) for f in input_files]

# Initialize netcdf file
nc_ckd = Dataset(conf['paths']['ckd_nitro'], 'w', format="NETCDF4")

# Global detector dimensions
for name in list(conf['dimensions'].keys()):
    dimsz = conf['dimensions'][name]
    nc_ckd.createDimension(name, dimsz)

# Global attributes
for name in list(conf['attributes'].keys()):
    nc_ckd.setncattr(name, conf['attributes'][name])
# Add date of creation
attrval = datetime.today().strftime('%Y-%m-%d')
nc_ckd.setncattr('date_created', attrval)

# Create groups, folder structure is used
dirckd = conf['paths']['dir_nitro']
conf_flat = flatten_dict(conf)
keylist = list(conf_flat.keys())  # Get keys in YAML as paths
for root, subnames, files in os.walk(dirckd):
    if root == dirckd or 'pycache' in root:
        continue
    group = root.replace(dirckd + '/', '')
    nc_ckd.createGroup(group)

    # Check if there are dimensions of this group in cfg file
    indices_in_cfg = [i for i, s in enumerate(keylist) if group in s]
    groupdims = [keylist[i] for i in indices_in_cfg]
    # Create group dimensions
    for groupdim in groupdims:
        dimname = groupdim.split('/')[-1]
        nc_ckd[group].createDimension(dimname, conf_flat[groupdim])

# Generate ckds, loop over folder structure
# every .py file will be a variable with that name
varpaths = glob(conf['paths']['dir_nitro'] + '/**/*[!pycache]*.py', recursive=True)
for varpath in varpaths:
    mod = import_mod(varpath)
    if conf['recalculate_ckds']: # first recalculate if chosen to
        try:
            mod.calculate(conf)
        except:
            pass
    gen = mod.generate(conf)  # Create generator
    nc_ckd = gen.add_to_nc(nc_ckd, conf) # Add to nc ckd file


# Set skips (temporary?)
skips = conf['skips']        
for skipname in skips.keys():
    skipvar = nc_ckd.createVariable(skipname, 'u1', [])
    skipvar = skips[skipname]

# TODO
# - implement temperature dependence/interpolation if we have data for detector
print('[done] >> {} created'.format(Path(conf['paths']['ckd_nitro'])))

nc_ckd.close()

# Copy file with date for version control
#dt = datetime.today().strftime('%Y%m%d')
#ckd_current_path = conf['paths']['ckd_nitro']
#ckd_version_path = ckd_current_path.replace('.nc', '_'+dt+'.nc')
#os.popen('cp {} {}'.format(ckd_current_path, ckd_version_path)) 