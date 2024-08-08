## Usage
Can be run as: 

`python ckd_generator_nitro.py config.yaml`. 

If no argument is given, the ckd_nitro.yaml in the CKD folder is used.

## Config
The config file (ckd_nitro.yaml) contains:

- Dimensions of ckd netcdf file, same structure as folders (see ckd_nitro.yaml as example)

- Global Attributes

- Some Single Value Parameters, to  be used in one or more ckds


## TANGO NITRO CKD GENERATOR
CKDs are generated using the folder structure in the nitro_ckds folder as 
groups. Variables are created using .py files in the folder using the 
same name as the file. A new CKD can easily be added by creating a .py
file in the group/folder of choice. Dimensions can be added to config file or using methods from nc_container.py

Outputs ckd_nitro.nc to CKD folder. 
