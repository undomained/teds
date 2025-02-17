import os, sys
import argparse
import yaml

import run_E2E

from teds import log
from teds.im.Python.input.input_yaml import Input_yaml
import teds.lib.lib_utils as Utils
import lib.data_netcdf.data_netcdf as dn

def cmdline(arguments):
    """             
        Get the command line arguments
        - arguments: command line arguments
                
        return: 
        - nominal_file: the nominal configuration file
        - scenario_file: the scenario file
    """         
                
    usage = """Run the scenario creator.
               The the nominal configuration is copied and updated with scenario information
            """

    nominal_help = """The nominal configuration file """
    parser = argparse.ArgumentParser(description= usage)
    parser.add_argument( "nominal_file", metavar="FILE", help=nominal_help)
    parser.add_argument( "scenario_file", metavar='FILE', help="The scenario file for which to create the configuration file for (and to run the E2E simulator.")
                    
    args = parser.parse_args(arguments)
    nominal_file = args.nominal_file
    scenario_file = args.scenario_file
    
    return nominal_file, scenario_file

def update_config(config, settings):
    """
        Find the settings in the config and update the config
    """

    for key, value in settings.items():
        if key not in config:
            # New info???????? Or wrong info?
            error_msg = f"Key:{key} not found in config. Can not be updated. Is it new info? Can not continue at the moment. Exiting"
            log.error(error_msg)
            sys.exit(error_msg)
        if isinstance(value, dict):
            update_config(config[key], value)
        else:
            config[key] = value

    return

def create_scenario_config(nominal_config, scenario):
    """
        Create scenario configuration based on nominal configuration
        Copy the nominal configuration.
        Add the scenario information
        overwrite nominal settings with scenario settings
    """

    # Copy nominal configuration
    scenario_config = nominal_config.copy()
    # Add scenario information
    scenario_config['scenario'] = {}
    scenario_config['scenario']['title'] = scenario['title']
    scenario_config['scenario']['description'] = scenario['description']
    scenario_config['scenario']['subdir'] = scenario['scenario_dir']
    scenario_config['scenario']['steps'] = scenario['steps']

    # Overwrite nominal settings with scenario settings
    scenario_specific_settings = scenario['scenario_config']
    update_config(scenario_config, scenario_specific_settings)

    return scenario_config

def build(nominal_config_file, scenario_file):
    """
        Read in the nominal configuration file
        Read in the scenario file
        Update the nominal configuration with the config settinges from the scenario file.
        The scenario configuration is written to file in subdirectory <cfg_path>/scenarios/scenario_string
        Where scenario_string is the name of the scenario file
        The scenario_config file also includes the steps that need to be rerun.
        The run_E2E build function is run for each requested step.
    """

    nominal_config_input = Input_yaml(nominal_config_file)
    nominal_config = nominal_config_input.read()

    scenario_input = Input_yaml(scenario_file)
    scenario = scenario_input.read()

    scenario_config = create_scenario_config(nominal_config, scenario)

    cfg_path = os.path.split(nominal_config_file)[0]

    scenario_string = os.path.splitext(os.path.basename(scenario_file))[0]
    scenario_config_path = os.path.join(cfg_path,'scenarios',scenario_string)
    scenario_config['header']['path'] = scenario_config_path
    if not os.path.exists(scenario_config_path):
        os.makedirs(scenario_config_path)
    scenario_config_file = os.path.join(scenario_config_path,f"full_config_{scenario_string}.yaml")
    print(f"Creating scenario config file: {scenario_config_file}")
    with open(scenario_config_file,'w') as outfile:
        yaml.dump(scenario_config, outfile)

    # run scenario:
    main_attribute_dict = Utils.get_main_attributes(scenario_config)
    steps = scenario_config['scenario']['steps']
    for step in steps:
        run_E2E.build(scenario_config, step, scenario_config_path, main_attribute_dict)

    return

if __name__ == "__main__":

    # TODO
    # txt file with scenarion names?
    # just the name of a scenario file?

    # Get input arguments
    nominal_file, scenario_file =  cmdline(sys.argv[1:])

    build(nominal_file, scenario_file)

