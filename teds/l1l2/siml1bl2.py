# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
import numpy as np
import numpy.typing as npt
import xarray as xr

from teds.gm.types import Geometry
from teds.l1l2.l1bl2 import write_l2
from teds.l1l2.types import L2
from teds.sgm.atmosphere import Atmosphere


# Helper function to generate noise based on type
def generate_noise(intensity: float,
                   base_values: npt.NDArray[np.float64],
                   noise_type: str) -> npt.NDArray[np.float64]:
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
    noise = np.empty(base_values.shape)
    for i in range(noise.size):
        noise.ravel()[i] = np.random.normal(0, intensity)
    if noise_type == 'absolute':
        # Absolute noise with given intensity
        return noise
    elif noise_type == 'relative':
        # Relative noise as a percentage of the base value
        return noise / 100 * base_values
    else:
        raise ValueError(
            "Invalid noise type. Choose 'absolute' or 'relative'.")


def simplified_level1b_to_level2_processor(config: dict) -> L2:

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
    nalt = ds.sizes['bins_along_track']
    nact = ds.sizes['bins_across_track']
    nlay = ds.sizes['layers']
    # Initialize the output array
    trace_gases = ['CO2', 'CH4', 'H2O']
    l2 = L2(nalt, nact, 0, nlay, 1, trace_gases)

    # Set the random seed for noise generation
    np.random.seed(seed)

    # Chi-square value (dummy value)
    l2.chi2[:] = 1
    # Convergence flag (dummy value)
    l2.converged[:] = 1
    # Number of iterations (dummy value)
    l2.number_iter[:] = 1
    # Albedo value from dataset
    l2.albedo_coeffs[0] = ds["albedo_b11"].values
    for gas in trace_gases:
        scale = {'CO2': 1.E6, 'CH4': 1.E9, 'H2O': 1.E6}
        l2.mixing_ratios[gas] = ds['x' + gas.lower()].values / scale[gas]
        l2.precisions[gas][:] = 0
        l2.col_avg_kernels[gas][:] = 1
    l2.proxys[gas] = l2.mixing_ratios[gas]
    l2.proxy_precisions[gas] = l2.precisions[gas]

    # If noise is needed, calculate noise for each gas type based on
    # its settings.
    if noise_needed:
        for gas in trace_gases:
            # Determine the base value from dataset
            noise = generate_noise(
                intensity_values['X' + gas],
                l2.mixing_ratios[gas],
                noise_type)
            l2.mixing_ratios[gas] += noise
            l2.precisions[gas][:] = noise
            if gas != 'H2O':
                l2.proxys[gas] += noise
                l2.proxy_precisions[gas][:] = noise

    atm = Atmosphere.from_empty()
    atm.zlay = ds["zlay"].values

    # Define retrieval initialization parameters
    retrieval_init = {
        # Surface pressure calculation
        'surface_pressure': ds["col_air"] / 6.02214076e23 * 0.029 * 9.81 / 100,
        # Surface elevation data
        'surface_elevation': ds["zlev"].values[-1],
        'trace_gases': {  # Reference profiles for trace gases
            'CO2': {'ref_profile': np.ones(nlay)},  # Dummy profile for CO2
            'CH4': {'ref_profile': np.ones(nlay)},  # Dummy profile for CH4
            'H2O': {'ref_profile': np.ones(nlay)}   # Dummy profile for H2O
        }
    }

    # Define Level 1B product data (latitude and longitude)
    geometry = Geometry.from_shape((nalt, nact))
    geometry.lat = ds["latitude"].values  # Latitude data from the dataset
    geometry.lon = ds["longitude"].values  # Longitude data from the dataset
    geometry.deg2rad()

    # Create the output file using the prepared data
    write_l2(filename, atm, l2, retrieval_init, geometry)

    print('=> l1bl2 finished successfully')

    return l2
