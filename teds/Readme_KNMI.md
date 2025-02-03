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
The different cases (scenarios) are in (scenario) subdirectories in directory cfg/nitro/scenarios.  
All settings for the different steps are combined in 1 yaml file (`full_config.yaml`).  
This ensures that the settings between te steps will be consistent.  
For the IM and L1AL1B steps a special config yaml file is created called `im_config_temp.yaml` and `l1al1b_config_temp.yaml` respectively.
Note: The nitro code uses for IM and L1B `driver_nitro`, which makes use of proctable and not of cal level!

### Building the executables.
The IM and L1B C++ code can be build for the first time running the following commands from the root directory:  
`mkdir build && cd build`  
`cmake -C ../initial_cache.cmake ..`  
`make -j `  
After that usually running only `make -j` in the `build` folder suffices.


### Running the E2E processor
In the teds directory the script `run_E2E.py` is situated.  
This can be run for a single step, multiple steps or all steps after oneother.
Get help on how to run:  
`python run_E2E.py -h`
#### Single step
For instance running only geometry step.  
`python run_E2E.py ../cfg/nitro/full_config.yaml gm`

#### Multiple steps
For instance running only geometry step.  
`python run_E2E.py ../cfg/nitro/full_config.yaml im l1al1b`

#### Running the full simulation:  
`python run_E2E.py ../cfg/nitro/full_config.yaml all`

#### Running the Python IM version
Update the config file and in the IM chapter set `do_python` to true.
`python run_E2E.py ../cfg/nitro/full_config.yaml im`
Output file is found in data/no2 directory and named `<l1a_name>_python.nc`, where `<l1a_name>` is defined in the configuration file.
The output of the inbetween steps can also be found as `im_l1x_<algo>.nc`.

### Output and input
Output files of the different steps is written in dorectory `data/no2` (or someother directory. you can change the location by updating the yaml file).  
The output file(s) of one step can be input files to the next step.  
There are also input files that are not produced by the different steps of the E2E processor. These are
* `ckd_nitro.nc`
* `binning_tables_no2.nc`
They can be found in `data/no2/ckd` directory (or some other directory. you can change the location by updating the yaml file).

### Simulation scenes
In the folder `cfg/nitro/scenes` the config files for the following MicroHH simulation scenes are given: Belchatow, Jaenschwalde, Lipetsk and Matimba.

### Running scenarios
When investigating effects of for instance changing temperature it is required to run a non nominal scenario.
First it is required to create a scenario file.
- In de main directory go to directory scenarios/nitro
- create a scenario yaml file: `vim <scenario_name>.yaml`
- Add the following information
  - title. Title for this scenario
  - description. Description for this scenario
  - `scenario_dir`. Sub directory in base output directory (data/no2) where the scenario output files can be found in sub directory scenarios.
  - steps. The steps that need to be rerun. For instance, if the scenario is running with different CKD files for 
    IM and L1B than there is no need to rerun GM and SGM. The nominal GM and SGM output files can be used.
    But from IM onwards all steps need to be rerun.
  - `scenario_config`. This holds all config settings that need to be different wrt nominal configuration. For instance ckd and `ckd_im`
    Structure needs to be the same as in nominal configuration file.
- in the teds directory create and run the new scenario: `python create_scenarios.py <nominal_config_file> <scenario_file>`, 
  where `<nominal_config_file>` is probably `cfg/nitro/full_config.py` and `<scenario_file>` is the scenario yaml file that is in directory `scenarios/nitro`
- output data is found in directory `data/no2/scenarios/<scenario_dir>`

### Creation of input files that come from outside the E2E processor (like ckd and binning table)
#### Creating binning table
Stand alone:
`python ckd/create_binning_table.py ckd/binning.yaml`

As part of E2E:
`python run_E2E.py ../cfg/nitro/full_config.yaml bin`

#### Creating CKD file
**nitro**
Stand alone:
`python ckd/ckd_generation/ckd_generator_nitro.py`, creates ckd_nitro.nc in CKD directory. 
A custom `config.yaml` can be added as an argument, otherwise the ckd_nitro.yaml in the ckd_generation directory is used.
More detailed information about the generator can be found in `ckd/ckd_generation/readme.md`.

As part of E2E:
`python run_E2E.py ../cfg/nitro/full_config.yaml ckd`

### DISAMAR radiative transfer
The Nitro SGM makes use of the KNMI DISAMAR software suite. This code is open source and available on https://gitlab.com/KNMI-OSS/disamar. See the DISAMAR readme for compiling instructions or run the provided building script in `sgm/build_disamar.sh`. The DISAMAR executable path has to be set in the Nitro configuration yaml file.


To run DISAMAR built with ifort on KNMI workstation the following shell script needs to be sourced in order to set the environment variables required for ifort:
`. /opt/intel/oneapi/compiler/latest/env/vars.sh`
