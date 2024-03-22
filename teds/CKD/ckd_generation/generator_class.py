"""
TANGO NITRO
Every ckd uses an instance of this class. Dimensions are 
added beforehand with the main ckd_generator_nitro.py and 
can be used as input by using the same names.

add_to_nc()
Adds variable and its attributes to the netcdf ckd file.

make_plot()
Makes a simple 1D or 2D plot and saves it to figures folder

"""

import matplotlib.pyplot as plt
from defs import save_figure
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
        if len(self.data.shape) == 2:
            newvar[:] = self.data  # Add data to variable
        elif len(self.data.shape) == 3:
            newvar[:,:] = self.data
        else:
            newvar = self.data
        return dataset
    
        
    # Plot function
    def make_plot(self, cfg):
        if not len(self.data):
            print("[generator_class] >> No data added to generator {}".format(self.name))
            return 0
        plt.rcParams.update({"text.usetex": False})
        print("[generator_class] >> plotting {}.png".format(self.name))
        fig, ax = plt.subplots(1,1, figsize = (8,8))
        if len(self.dim_names) == 1:
            ax.plot(self.data)
            ax.set_xlabel(self.dim_names[0])
            ax.set_ylabel(self.name)
            ax.set_title(self.name)
            save_figure(fig, self.name, cfg, formats = ['png'], v = False)
        elif len(self.dim_names) == 2:
            im = ax.pcolormesh(self.data, cmap = 'cubehelix', rasterized = True)
            ax.set_xlabel(self.dim_names[1])
            ax.set_ylabel(self.dim_names[0])
            ax.set_title(self.name)
            cbar = fig.colorbar(im)
            if 'units' in self.attr_names:
                ix = np.flatnonzero(np.array(self.attr_names) == 'units')[0]
                units = self.attr_vals[ix]
                unitlbl = '' if units == '1' else units
                cbar.set_label(unitlbl, rotation=270, va = 'bottom')
            save_figure(fig, self.name, cfg, formats = ['png'], v = False)
        else:
            print(self.name, 'cannot plot for this shape')
        


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