# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
import numpy as np
import numpy.typing as npt
import xarray as xr

from teds.l1l2.l1bl2 import write_l2


# Helper function to generate noise based on type
def generate_noise(intensity: float, base_value: float, noise_type: str) -> (
        float):
    """Generates noise based on the specified type.

    Parameters
    ----------
    intensity
        The standard deviation of the noise.
    base_value
        The base value to apply relative noise to.
    noise_type
        Type of noise, either 'absolute' or 'relative'.

    Returns
    -------
        Noise value to be added to the data.

    """
    if noise_type == 'absolute':
        # Absolute noise with given intensity
        return np.random.normal(0, intensity)
    elif noise_type == 'relative':
        # Relative noise as a percentage of the base value
        return np.random.normal(0, intensity) / 100 * base_value
    else:
        raise ValueError(
            "Invalid noise type. Choose 'absolute' or 'relative'.")


def simplified_level1b_to_level2_processor(config: dict) -> npt.NDArray:

    # 0. Load input data file
    ref_filename = config["io_files"]["input_geo_ref"]
    ds = xr.open_dataset(ref_filename)  # Load the NetCDF dataset

    # Extract noise settings from the configuration
    noise_needed = config['sim_with_noise']
    # Type of noise: 'absolute' or 'relative'
    noise_type = config["noise_model"]["type"]
    # Noise intensities for different gases
    intensity_values = {
        'XCO2': config["noise_model"]['intensity_xco2'],
        'XCH4': config["noise_model"]['intensity_xch4'],
        'XH2O': config["noise_model"]['intensity_xh2o'],
        'XCO2_proxy': config["noise_model"]['intensity_xco2_proxy'],
        'XCH4_proxy': config["noise_model"]['intensity_xch4_proxy']
    }
    seed = config["noise_model"]['seed']  # Random seed for reproducibility

    # Define the output filename
    filename = config["io_files"]["output_siml2"]

    # Prepare the Level 2 product: a 2D array of dictionaries
    # containing retrieval data.
    # Number of along-track bins, across-track bins, and layers
    nalt, nact, nlay = ds["zlay"].shape
    # Initialize the output array
    l2product = np.empty((nalt, nact), dtype=object)

    # Set the random seed for noise generation
    np.random.seed(seed)

    # Iterate over all bins to fill the Level 2 product data
    for ialt in range(nalt):
        for iact in range(nact):
            # Initialize noise for each gas and proxy
            noise = {}
            if noise_needed:
                # If noise is needed, calculate noise for each gas
                # type based on its settings.
                for key in intensity_values:
                    # Determine the base value from dataset
                    base_value = ds[key.split('_')[0]].values[ialt, iact]
                    noise[key] = generate_noise(intensity_values[key],
                                                base_value, noise_type)
            else:
                # No noise case: set noise for each gas type to zero
                noise = {key: 0 for key in intensity_values}

            # Fill the Level 2 product dictionary for each bin
            l2product[ialt, iact] = {
                'convergence': 1,  # Convergence flag (dummy value)
                'number_iter': 1,  # Number of iterations (dummy value)
                'chi2': 1,  # Chi-square value (dummy value)
                # Albedo value from dataset
                'alb0': ds["albedo"].values[ialt, iact],
                'spec_shift': 0,  # Spectral shift (dummy value)
                'spec_squeeze': 0,  # Spectral squeeze (dummy value)
                # XCO2 value with noise
                'XCO2': ds["XCO2"].values[ialt, iact] + noise['XCO2'],
                'XCO2 precision': noise['XCO2'],  # XCO2 noise precision
                # XCO2 column average kernel (dummy array)
                'XCO2 col avg kernel': np.ones(nlay),
                # XCH4 value with noise
                'XCH4': ds["XCH4"].values[ialt, iact] + noise['XCH4'],
                'XCH4 precision': noise['XCH4'],  # XCH4 noise precision
                # XCH4 column average kernel (dummy array)
                'XCH4 col avg kernel': np.ones(nlay),
                # XH2O value with noise
                'XH2O': ds["XH2O"].values[ialt, iact] + noise['XH2O'],
                'XH2O precision': noise['XH2O'],  # XH2O noise precision
                # XH2O column average kernel (dummy array)
                'XH2O col avg kernel': np.ones(nlay),
                # XCO2 proxy value with noise
                'XCO2 proxy': (ds["XCO2"].values[ialt, iact]
                               + noise['XCO2_proxy']),
                # XCO2 proxy noise precision
                'XCO2 proxy precision': noise['XCO2_proxy'],
                # XCH4 proxy value with noise
                'XCH4 proxy': (ds["XCH4"].values[ialt, iact]
                               + noise['XCH4_proxy']),
                # XCH4 proxy noise precision
                'XCH4 proxy precision': noise['XCH4_proxy'],
            }

    # Define retrieval initialization parameters
    retrieval_init = {
        'zlay': ds["zlay"].values[0, 0, :],  # Layer heights (dummy data)
        # Surface pressure calculation
        'surface pressure': ds["col_air"] / 6.02214076e23 * 0.029 * 9.81 / 100,
        # Surface elevation data
        'surface elevation': ds["zlev"].values[:, :, -1],
        'trace gases': {  # Reference profiles for trace gases
            'CO2': {'ref_profile': np.ones(nlay)},  # Dummy profile for CO2
            'CH4': {'ref_profile': np.ones(nlay)},  # Dummy profile for CH4
            'H2O': {'ref_profile': np.ones(nlay)}   # Dummy profile for H2O
        }
    }

    # Define Level 1B product data (latitude and longitude)
    l1bproduct = {
        'latitude': ds["lat"].values,  # Latitude data from the dataset
        'longitude': ds["lon"].values,  # Longitude data from the dataset
    }

    # Additional settings (currently empty but can be used for future
    # extensions)
    # Currently not used, so an empty dictionary is provided
    settings: dict = {}

    # Create the output file using the prepared data
    write_l2(filename,
             l2product,  # type: ignore
             retrieval_init,  # type: ignore
             l1bproduct,
             settings)  # type: ignore

    print('=> l1bl2 finished successfully')

    return l2product
