import os
import sys
import argparse
import importlib
import subprocess
import teds.lib.lib_utils as Utils


def cmdline(arguments):
    """
        Get the command line arguments
        - arguments: command line arguments

        return:
        - cfgFile: the configuration file
    """

    usage = """Run the TANGO E2E processor.
               The configuration file contains the settings for each step in
    the E2E processor."""

    cfgHelp = """The configuration file needed to run the E2E processor.
    Possible choices: Geomery
              (gm), Scene Generation (sgm), Instrument model (im), L1A to L1B
    (l1al1b), L1B to L2
              (l1l2) or all steps (all)."""
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("cfgFile", metavar="FILE", help=cfgHelp)
    parser.add_argument("step",
                        metavar='STEP',
                        choices=['gm', 'sgm', 'im', 'l1al1b', 'l1l2', 'all'],
                        help="The steps that the E2E processor has to run.")

    args = parser.parse_args(arguments)
    cfgFile = args.cfgFile
    step = args.step

    return cfgFile, step


def build(logger, config, step, cfg_path, attribute_dict):
    """
        Run E2E processor.
        - logger: Reference to the program logger
        - config: configuration file containing the settings for the different
                  steps in the E2E processor
        - step: indicating which step in the E2E processor to run.
        - cfg_path: configuration path
        - attribute_dict: Dictionary with attributes to be added to main of
                          output netCDF files
    """

    if step == 'gm' or step == 'all':

        E2EModule = importlib.import_module("GM.gm")
        E2EModule.geometry_module(logger, config)
        # add attributes to the output file
        Utils.add_attributes_to_output(
            logger, config['gm_file'], attribute_dict)

    if step == 'sgm' or step == 'all':
        E2EModule = importlib.import_module("SGM.sgm_no2")
        E2EModule.scene_generation_module_nitro(logger, config)
        Utils.add_attributes_to_output(
            logger, config['sgm_rad_file'], attribute_dict)
        Utils.add_attributes_to_output(
            logger, config['sgm_atm_file'], attribute_dict)

    if step == 'im' or step == 'all':
        # Create cfg file to be used for IM executable
        E2EModule = importlib.import_module(
            "IM.create_im_configuration_file_nitro")
        E2EModule.im_configuration(config)
        # Need to call C++
        subprocess.run(["IM/tango_ckd_model/build/ckdmodel",
                        f"{cfg_path}/im_config.cfg"])
        Utils.add_attributes_to_output(
            logger, config['l1a_file'], attribute_dict)

    if step == 'l1al1b' or step == 'all':
        # Create cfg file to be used for L1AL1B executable
        E2EModule = importlib.import_module(
            "L1AL1B.create_l1a1b_configuration_file_nitro")
        E2EModule.l1al1b_configuration(config)
        # Need to call C++
        subprocess.run(["L1AL1B/tango_l1b/build/tango_l1b",
                        f"{cfg_path}/l1al1b_config.cfg"])
        Utils.add_attributes_to_output(
            logger, config['l1b_file'], attribute_dict)

    if step == 'l1l2' or step == 'all':
        E2EModule = importlib.import_module("L1L2.l1bl2_no2")
        E2EModule.l1bl2_no2(logger, config)
        Utils.add_attributes_to_output(
            logger, config['l2_file'], attribute_dict)

    if step == 'pam' or step == 'all':
        E2EModule = importlib.import_module("PAM.pam")
        E2EModule.pam_nitro(logger, config)


if __name__ == "__main__":

    # Get logger for run script
    build_logger = Utils.get_logger()

    # Get input arguments
    cfgFile, step = cmdline(sys.argv[1:])
    cfg_path, filename = os.path.split(cfgFile)

    # Get configuration info
    config = Utils.getConfig(build_logger, cfgFile)
    config['header']['path'] = cfg_path

    # Get information (like git hash and config file name and version
    # (if available) that will be added to the output files as
    # attributes
    main_attribute_dict = Utils.get_main_attributes(config)

    build(build_logger, config, step, cfg_path, main_attribute_dict)
