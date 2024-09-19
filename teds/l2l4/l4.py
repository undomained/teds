import xarray as xr
import numpy as np
from l4_tools import *

def level2_to_level4_processor(config):
    """
    Process Level 2 to Level 4 emission data based on the provided configuration.
    
    Parameters:
    config (dict): A dictionary containing file paths and processing options.

    Returns:
    dict: The processed Level 4 product grouped by the estimation method.
    """
    
    sgm_filename = config["io_files"]["input_sgm"]
    sgm_ref_filename = config["io_files"]["input_sgm_ref"]
    l2_filename = config["io_files"]["input_l2"]
    l4_filename = config["io_files"]["output_l4"]
    methods = config["emission_estimation"]["method"]
    species = config["emission_estimation"]["species"]

    l4_product = {}

    # Iterate over all methods in the configuration
    for method in methods:
        print(f"Processing with method: {method}")
        method_group = {}  # Store the results for the current method

        # Emission estimation for each species based on the method
        for sp in species:
            apr_name = sp.split(" ")[0]  # Derive apriori variable name

            if method == "microhh_fit":
                # Run the existing microhh_fit function
                apr_emission, ret_emission, ret_emission_error = microhh_fit(
                    sgm_filename, sgm_ref_filename, l2_filename, sp, apr_name
                )
                method_group[sp] = (apr_emission, ret_emission, ret_emission_error)

            # Additional methods can be added here with elif blocks
            # elif method == "other_method":
            #     # Placeholder for future emission estimation methods
            #     ...

            else:
                print(f"Warning: Method {method} is not implemented yet.")
                continue

        # Add the results of this method to the overall l4_product
        if method_group:
            l4_product[method] = method_group

    # Create the output NetCDF file with grouped data
    create_grouped_dataset(l4_filename, l4_product)

    print('=> Level 2 to Level 4 processing finished successfully')

    return l4_product
