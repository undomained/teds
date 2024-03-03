# Tango end-to-end simulator
Code to run the E2E simulation

License
-------

This project is distributed under the terms of the 3-Clause BSD License. See LICENSE in the root directory of the project or go to https://opensource.org/licenses/BSD-3-Clause.

## The Code
Checkout the code:
`git xxxxxxx`

## Working with virtual environment
Create a virtual environment:
`python -m venv venv_E2E`  (the last venv is name of directory)

### Starting virtual environment:
`source venv_E2E/bin/activate`

### Deactivate virtual environment:
`deactivate`

### PIP
While in virtual environment:
Install/update pip:
`python3 -m pip install --upgrade pip`
Check it:
`python3 -m pip --version`

### Install pakages in virtual env:
    - after activation just use pip install to install whatever package you want
      OR use requirement.txt:
        - Create a requirements.txt file:
          after activation of virtual environment and installing all packages needed:
          `pip freeze > requirements.txt`
        - Use the requirements.txt:
          after activation of virtual environment:
          `pip install -r requirements.txt`

### Set python path:
From the teds directory:
`export PYTHONPATH=$PWD/..:$PYTHONPATH`

## How to run a nominal case for Tango Nitro

### Configuration files
Configuration files can be found in cfg/nitro.  
The configuration file for the nominal run is in this directory.  
The different cases are in subdirectories.  
All settings for the different steps are combined in 1 yaml file.  
This ensures that the settings between te steps will be consistent.  
Note: the .cfg files needed to run the C++ code in the IM and L1AL1B part are generated from the yaml input when step is `im`, `l1al1b` or `all`.  
Where the .cfg file is written is determined in the yaml file settings, but python scripts expects this at the moment to be in the same location as the yaml file.

### Building the executables.
See `build_instructions.md` in the teds directory

### Running the E2E processor
In the teds directory the script `run_E2E.py` is situated.  
This can be run for a single step or all steps after oneother.
Get help on how to run:  
`python run_E2E.py -h`
#### Single step
For instance running only geometry step.  
`python run_E2E.py ../cfg/nitro/full_config.yaml gm`

#### Running the full simulation:  
`python run_E2E.py ../cfg/nitro/full_config.yaml all`

### Output and input
Output files of the different steps is written in dorectory `data/no2` (or someother directory. you can change the location by updating the yaml file).  
The output file(s) of one step can be input files to the next step.  
There are also input files that are not produced by the different steps of the E2E processor. These are
* `ckd.nc`
* `binning_tables.nc`
They can be found in `data/no2/ckd` directory (or some other directory. you can change the location by updating the yaml file).

### Creation of input files that come from outside the E2E processor (like ckd and binning table)
#### Creating binning table
`python CKD/create_binning_table.py`

#### Creating CKD file
TBW


### DISAMAR radiative transfer
The Nitro SGM makes use of the KNMI DISAMAR software suite. This code is open source and available on https://gitlab.com/KNMI-OSS/disamar. See the DISAMAR readme for compiling instructions. The DISAMAR executable path then has to be given in the Nitro SGM yaml file.
