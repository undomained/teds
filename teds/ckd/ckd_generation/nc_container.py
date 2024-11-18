import numpy as np
import inspect
import os
from datetime import datetime
from netCDF4 import Dataset
import importlib.util

class NC_container:
    def __init__(self, cfg):
        self.nc = Dataset(cfg['paths']['ckd_nitro'], 'w', format="NETCDF4")
        self.cfg = cfg

        # Create groups
        dirckd = cfg['paths']['dir_nitro']
        for root, _, _ in os.walk(dirckd):
            if root == dirckd or '__pycache__' in root:
                continue
            group_name = root.replace(f"{dirckd}/", '')
            self.nc.createGroup(group_name)
        
        # Create dimensions
        self.create_dims_from_dict(self.nc, cfg['dimensions'])
        for name, value in cfg['attributes'].items():
            self.nc.setncattr(name, value)
            self.nc.setncattr('date_created', datetime.today().strftime('%Y-%m-%d'))
            
        # Create dictionary with created varnames and grouppaths
        self.vars = []
    
    def create_dims_from_dict(self, ncgroup, dic):
        for name, val in dic.items():
            if isinstance(val, dict):
                if name not in ncgroup.groups:
                    print(f'No folder/group named {name}')
                    continue
                self.create_dims_from_dict(ncgroup[name], val)
            else:
                ncgroup.createDimension(name, val)
    
    def create_dims_auto(self, dims):
        groupname = self.get_group_name()
        ncgroup = self.nc[groupname]
        self.create_dims_from_dict(ncgroup, dims)

    def create_dim(self, dimname, size):
        groupname = self.get_group_name()
        ncgroup = self.nc[groupname]
        ncgroup.createDimension(dimname, size)


    def create_var(self, name, dims, values, attrs = {}, dtype = 'u1'):
        group_name = self.get_group_name()
        varpath = f"{group_name}/{name}"
        if varpath in self.vars: # no double calculations
            return
        
        print(f'Generating {name}...')    
        var = self.nc.createVariable(varpath, dtype, dims)
        self.vars.append(varpath)
        for attr_name, attr_val in attrs.items():
            var.setncattr(attr_name, attr_val)
        var[:] = values
        self.name = var
        
        
    def create_var_auto(self, dims, values, attrs = {}, dtype = 'u1'):
        name = self.get_varname_from_file()
        self.create_var(name, dims, values, attrs, dtype)
    
    def get_var(self, varpath):
        if varpath not in self.vars:
            dirpath = self.cfg['paths']['dir_nitro'].replace('./', '')
            varfilepath = f'{dirpath}/{varpath}.py'
            mod = import_mod(varfilepath)
            mod.generate(self)
        var = self.nc[varpath]
        return np.array(var[:])

    def get_varname_from_file(self):
        path = self.get_path_ckd_file()
        j = path.rfind('/') + 1
        name = path[j:].replace('.py','')
        return name

    def get_path_ckd_file(self):
        frame = inspect.currentframe()
        outframe = inspect.getouterframes(frame)
        dirpath = self.cfg['paths']['dir_nitro'].replace('./', '')
        for i, outf in enumerate(outframe):
            fname = os.path.relpath(outf.filename)
            if dirpath in fname:
                i = fname.find(dirpath) + len(dirpath)
                return fname[i:]

    def get_group_name(self):
        path = self.get_path_ckd_file()
        j = path.rfind('/')
        group = path[:j].replace('/','')
        return group
    
    def get_shape(self, dims):
        dimshape = []
        for dim in dims:
            if dim not in self.cfg['dimensions'].keys():
                grp = self.get_group_name()
                dimshape.append(self.cfg['dimensions'][grp][dim])
            else:
                dimshape.append(self.cfg['dimensions'][dim])
        return dimshape
        #return tuple(self.cfg['dimensions'][dim] for dim in dims)
         
    def close(self):
        print(f"[done] >> {self.cfg['paths']['ckd_nitro']} created")
        self.nc.close()
        
    

def import_mod(varpath_py):
    module_name = os.path.basename(varpath_py).split('.')[0]
    spec = importlib.util.spec_from_file_location(module_name, varpath_py)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module