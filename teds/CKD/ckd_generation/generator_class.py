"""
TANGO NITRO
Every ckd uses an instance of this class. Dimensions are 
added beforehand with the main ckd_generator_nitro.py and 
can be used as input by using the same names.

add_to_nc()
Adds variable and its attributes to the netcdf ckd file.

"""

from pathlib import Path
import inspect
import numpy as np
import importlib.util
import sys
import os


class ckd_generator():
    def __init__(self):
        self.data = []          # data of the variable
        self.dtype = 'u1'       # data type
        self.dim_names = []     # names of the dimensions
        self.attr_names = []    # attributes names, can be anything
        self.attr_vals = []     # values of the attributes, same 
                                # order as names
        
        # Get name of current variable by filename
        frame = inspect.currentframe()
        outframe = inspect.getouterframes(frame)
        self.path = outframe[1].filename
        self.name = Path(outframe[1].filename).stem


    # Function to add to ckd NETCDF dataset
    def add_to_nc(self, dataset, cfg):
        # Get group name

        path = os.path.relpath(self.path)
        dirpath = cfg['paths']['dir_nitro'].replace('./','')
        group = path.replace(dirpath,'').replace(self.name+'.py', '')
        group = group.replace('/','')
        print("[ckd_generator] >> generating {}".format(self.name))
        # Set Variables
        if len(group):  # Group Variables
            newvar = dataset[group].createVariable(self.name, self.dtype, self.dim_names)
        else:  # Global Variables
            newvar = dataset.createVariable(self.name, self.dtype, self.dim_names) 
        # Set attributes
        for a, attrname in enumerate(self.attr_names):  
            newvar.setncattr(attrname, self.attr_vals[a])

        self.data = np.array(self.data)
        if len(self.data.shape) <= 2:
            newvar[:] = self.data  # Add data to variable
        elif len(self.data.shape) == 3:
            newvar[:,:] = self.data
        else:
            newvar = self.data
        return dataset
        


def import_mod(varpath_py):
    """
    Import a module from path
    """
    module_name = Path(varpath_py).stem
    spec = importlib.util.spec_from_file_location(module_name, varpath_py)
    foo = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = foo
    spec.loader.exec_module(foo)
    return foo


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