## Usage
Can be run as: 

`python ckd_generator_nitro.py config.yaml`. 

If no argument is given, the ckd_nitro.yaml in the CKD folder is used.

## Config
The config file (ckd_nitro.yaml) contains:

- Dimensions of ckd netcdf file, same structure as folders (see ckd_nitro.yaml as example)

- Global Attributes

- Some Single Value Parameters, to  be used in one or more ckds

- Toggles for making figures, or recalculating all ckds if new ckds have been added. 


## TANGO NITRO CKD GENERATOR
CKDs are generated using the folder structure in the nitro_ckds folder as 
groups. Variables are created using .py files in the folder using the 
same name as the file. A new CKD can easily be added by creating a .py
file in the group/folder of choice. Dimensions can be added to config file.
In the .py files:

 - **generate()** method: returns an instance of the generator_class
   which is used to standardize the output and add the variable to the 
   ckd nc file. Here you can add the dimensions, attributes and values.
   Needs the config yaml as argument

 - **(optional) import_data()** method: when data is imported from external
   file which also contains a lot of other data. This method can be used to 
   export a file with only the useful data for this ckd. Ideally this method
   would be obsolete if the external data would be either standardized or final.
   Needs config yaml as argument

 - **(optional) calculate()** method if an elaborate or time consuming
   calculation is needed. It is useful if this exports a file with the 
   data which the gerenate() method can read to prevent recalculating
   all ckds every time. (Can be toggled in config yaml)

Outputs ckd_nitro.nc to CKD folder. 
