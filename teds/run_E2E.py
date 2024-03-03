import os, sys
import argparse
import logging
import yaml
import importlib
import subprocess
import lib.data_netcdf.data_netcdf as dn

def cmdline(arguments):
    """             
        Get the command line arguments
        - arguments: command line arguments
                
        return: 
        - cfgFile: the configuration file
    """         
                
    usage = """Run the TANGO E2E processor.
               The configuration file contains the settings for each step in the E2E processor."""

    cfgHelp = """The configuration file needed to run the E2E processor. Possible choices: Geomery
              (gm), Scene Generation (sgm), Instrument model (im), L1A to L1B (l1al1b), L1B to L2
              (l1l2) or all steps (all)."""
    parser = argparse.ArgumentParser(description= usage)
    parser.add_argument( "cfgFile", metavar="FILE", help=cfgHelp)
    parser.add_argument( "step", metavar='STEP', choices=['gm','sgm','im','l1al1b','l1l2','all'],
                       help="The steps that the E2E processor has to run.")
                    
    args = parser.parse_args(arguments)
    cfgFile = args.cfgFile
    step = args.step
    
    return cfgFile, step

def get_logger():
    """
       Gets or creates a logger
    """

    log_level = logging.INFO
    log_format = '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s'
    log_formatter = logging.Formatter(log_format)
    date_format = '%d/%m/%Y %H:%M:%S'

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    c_handler = logging.StreamHandler(sys.stdout)
    c_handler.setFormatter(log_formatter)
    logger.addHandler(c_handler)

    return logger

def getConfig(logger, cfgFile):
    """
        Get the config information from the configuration file
       - logger: Reference to the program logger
       - cfgFile: configuration file
       return:
       - configuration: configuration info
    """
    stream =  open(cfgFile, 'r')
    config = yaml.safe_load(stream)

   #TODO do we need the below capability?
   # Fill in variables in the configuration
    for key in config:
        if isinstance(config[key], str):
            config[key] = config[key].format(**config)
        if isinstance(config[key], list):
            for i, item in enumerate(config[key]):
                if isinstance(item, str):
                    config[key][i] = config[key][i].format(**config)

    config_string = "Configuration used in this step of the analysis:\n"
    for key in config:
        config_string += f"    {key} = {config[key]}\n"

    logger.info(config_string)

    return config

def get_main_attributes(config):
    """
        Define some info to add as attributes to the main of the ouput netCDF file
    """

    attribute_dict = {}

    attribute_dict['git_hash'] = str(subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip())
    attribute_dict['git_hash_short'] = str( subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip())

    attribute_dict['cfg_file'] = config['header']['file_name']
    attribute_dict['cfg_version'] = config['header']['version']
    config.pop('header')
    attribute_dict['E2E_configuration'] = str(config)
    #TODO: Add other information that might be handy to have in the attributes of the netCDF output file.

    return attribute_dict

def add_attributes_to_output(logger, output_file, attribute_dict):
    """
        Add attributes to the output file
    """

    out_data = dn.DataNetCDF(logger, output_file, mode='r')
    for name, value in attribute_dict.items():
        out_data.add(name, value=value, kind='attribute')
    out_data.write()
    return

def build(logger, config, step, cfg_path, attribute_dict):
    """
        Run E2E processor.
        - logger: Reference to the program logger
        - config: configuration file containing the settings for the different steps in the E2E processor
        - step: indicating which step in the E2E processor to run. 
        - cfg_path: configuration path
        - attribute_dict: Dictionary with attributes to be added to main of output netCDF file
    """

    if step == 'gm' or step == 'all':

        E2EModule = importlib.import_module("GM.gm")
        E2EModule.geometry_module(config)
        # add attributes to the output file
        add_attributes_to_output(logger, config['gm_file'], attribute_dict)

    if step == 'sgm' or step == 'all':
        #TODO need to be filled in
        E2EModule = importlib.import_module("SGM.sgm")
        E2EModule.scene_generator_module(config)
        add_attributes_to_output(logger, config['sgm_file'], attribute_dict) #?

    if step == 'im' or step == 'all':
        # Create cfg file to be used for IM executable
        E2EModule = importlib.import_module("IM.create_im_configuration_file_nitro")
        E2EModule.im_configuration(config)
        # Need to call C++
#        output = subprocess.run(["IM/tango_ckd_model/build/ckdmodel", "../cfg/nitro/im_config.cfg"], stdout = subprocess.PIPE, universal_newlines = True).stdout
#        subprocess.run(["IM/tango_ckd_model/build/ckdmodel", "../cfg/nitro/im_config.cfg"])
        subprocess.run(["IM/tango_ckd_model/build/ckdmodel", f"{cfg_path}/im_config.cfg"])
        add_attributes_to_output(logger, config['l1a_file'], attribute_dict)

    if step == 'l1al1b' or step == 'all':
        # Create cfg file to be used for L1AL1B executable
        E2EModule = importlib.import_module("L1AL1B.create_l1a1b_configuration_file_nitro")
        E2EModule.l1al1b_configuration(config)
        # Need to call C++
#        subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b", "../cfg/nitro/l1al1b_config.cfg"])
        subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b", f"{cfg_path}/l1al1b_config.cfg"])
        add_attributes_to_output(logger, config['l1b_file'], attribute_dict)

    if step == 'l1l2' or step == 'all':
        E2EModule = importlib.import_module("L1L2.l1l2")
        E2EModule.level1b_to_level2_processor(config)
        add_attributes_to_output(logger, config['l2_file'], attribute_dict) #?



if __name__ == "__main__":
    
    cfgFile, step =  cmdline(sys.argv[1:])
    cfg_path, filename = os.path.split(cfgFile)

    build_logger = get_logger()

    config = getConfig(build_logger, cfgFile)

    main_attribute_dict = get_main_attributes(config)

    build(build_logger, config, step, cfg_path, main_attribute_dict)

