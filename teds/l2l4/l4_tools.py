from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
from scipy.optimize import curve_fit
import numpy as np
import xarray as xr
import numpy as np

def create_meshgrid(lat, lon, z):
    """
    Create meshgrid for 2D lat, 2D lon, and 1D z arrays.
    """
    # Number of layers corresponds to the size of z
    z_grid = np.repeat(z[np.newaxis, np.newaxis, :], lat.shape[0], axis=0)
    z_grid = np.repeat(z_grid, lat.shape[1], axis=1)
        
    lat_grid = np.repeat(lat[:, :, np.newaxis], len(z), axis=2)
    lon_grid = np.repeat(lon[:, :, np.newaxis], len(z), axis=2)
        
    return lat_grid, lon_grid, z_grid


def regrid_3d(u, lat, lon, z, lat_new, lon_new, z_new, method='nearest'):
    """
    Regrid a 3D matrix `u` defined on (lat, lon, z) to a new grid (lat_new, lon_new, z_new).

    Parameters:
    - u (numpy.ndarray): 3D matrix of shape (len(lat), len(lon), len(z)).
    - lat (numpy.ndarray): 1D array of original latitudes.
    - lon (numpy.ndarray): 1D array of original longitudes.
    - z (numpy.ndarray): 1D array of original altitudes.
    - lat_new (numpy.ndarray): 1D array of new latitudes.
    - lon_new (numpy.ndarray): 1D array of new longitudes.
    - z_new (numpy.ndarray): 1D array of new altitudes.
    - method (str): Interpolation method ('linear', 'nearest', etc.).

    Returns:
    - u_new (numpy.ndarray): Regridded 3D matrix of shape (len(lat_new), len(lon_new), len(z_new)).
    """

    # Create meshgrid for the original grid
    if len(lat.shape)==1:
        lat_grid, lon_grid, z_grid = np.meshgrid(lat, lon, z, indexing='ij')
    else:
        lat_grid, lon_grid, z_grid =create_meshgrid(lat, lon, z)

    # Create meshgrid for the new grid
    if len(lat_new.shape)==1:
        lat_new_grid, lon_new_grid, z_new_grid = np.meshgrid(lat_new, lon_new, z_new, indexing='ij')
    else:
        lat_new_grid, lon_new_grid, z_new_grid = create_meshgrid(lat_new, lon_new, z_new)

    # Reshape and stack to prepare points for interpolation
    points = np.column_stack((lat_grid.ravel(), lon_grid.ravel(), z_grid.ravel()))

    # Interpolator function for u
    interpolator = RegularGridInterpolator((lat, lon, z), u, method=method, bounds_error=False, fill_value=None)

    # Prepare the points for interpolation onto the new grid
    points_new = np.stack((lat_new_grid.ravel(), lon_new_grid.ravel(), z_new_grid.ravel()), axis=-1)

    # Interpolate u onto the new grid
    if len(lat_new.shape)==1:
        u_new = interpolator(points_new).reshape(len(lat_new), len(lon_new), len(z_new))
    else:
        u_new = interpolator(points_new).reshape(lat_new_grid.shape)

    return u_new

def calculate_divergence_test(Fx, Fy, dx, dy, method='fourth_order'):
    """
    Calculate the divergence of the vector field using the specified method.
    
    Args:
    - Fx: 2D numpy array of the x-component of the vector field.
    - Fy: 2D numpy array of the y-component of the vector field.
    - dx: Grid spacing in the x-direction (in meters).
    - dy: Grid spacing in the y-direction (in meters).
    - method: 'second_order' or 'fourth_order', method for calculating derivatives.
    
    Returns:
    - divergence: 2D numpy array of the divergence of the vector field.
    """
    if method == 'second_order':
        # Second-order central difference
        dFx_dx = (Fx[2:, 1:-1] - Fx[:-2, 1:-1]) / (2 * dx)
        dFy_dy = (Fy[1:-1, 2:] - Fy[1:-1, :-2]) / (2 * dy)
        
        # Ensure the dimensions match
        min_rows = min(dFx_dx.shape[0], dFy_dy.shape[0])
        min_cols = min(dFx_dx.shape[1], dFy_dy.shape[1])
        
        divergence = np.zeros((min_rows, min_cols))
        divergence = dFx_dx[:min_rows, :min_cols] + dFy_dy[:min_rows, :min_cols]
    
    elif method == 'fourth_order':
        # Fourth-order central difference
        if Fx.shape[0] < 4 or Fx.shape[1] < 4:
            raise ValueError("Input arrays must be at least 4x4 for fourth-order central differences")
        
        #dFx_dx = (Fx[2:, 2:-2] - 8 * Fx[1:-1, 2:-2] + 8 * Fx[:-2, 2:-2] - Fx[:-2, 2:-2]) / (12 * dx)
        #dFy_dy = (Fy[2:-2, 2:] - 8 * Fy[2:-2, 1:-1] + 8 * Fy[2:-2, :-2] - Fy[2:-2, :-2]) / (12 * dy)

        dFx_dx= (-Fx[2:-2, 4:] + 8*Fx[2:-2, 3:-1] - 8*Fx[2:-2, 1:-3] + Fx[2:-2, 0:-4])/(12*dx)
        dFy_dy= (-Fy[4:, 2:-2] + 8*Fy[3:-1, 2:-2] - 8*Fy[1:-3, 2:-2] + Fy[0:-4, 2:-2])/(12*dy)

        #dFx_dx=np.gradient(Fx,axis=1,edge_order=1)/12*dx
        #dFy_dy=np.gradient(Fy,axis=0,edge_order=1)/12*dy
        
        if False:
            from scipy.ndimage import convolve
            #kernel_x=np.array([[1,-2,1]])/(4*dx)
            #kernel_y=np.array([[1],[-2],[1]])/(4*dx)

            f=Fx*4+Fy*4

            kernel_x=np.array([[-1,16,-30,16,-1]])/(12*dx)
            kernel_y=np.array([[-1],[16],[-30],[16],[-1]])/(12*dx)
            dFx_dx=convolve(f,kernel_x,mode="reflect")
            dFy_dy=convolve(f,kernel_y,mode="reflect")


        # Ensure the dimensions match
        min_rows = min(dFx_dx.shape[0], dFy_dy.shape[0])
        min_cols = min(dFx_dx.shape[1], dFy_dy.shape[1])
        
        divergence = np.zeros((min_rows, min_cols))
        divergence = dFx_dx[:min_rows, :min_cols] + dFy_dy[:min_rows, :min_cols]
    
    else:
        raise ValueError("Method must be 'second_order' or 'fourth_order'")
    
    return divergence,dFx_dx[:min_rows, :min_cols], dFy_dy[:min_rows, :min_cols]

def calculate_divergence(Fx, Fy, dx, dy):
    """
    Calculate the divergence of the vector field using the specified method.
    
    Args:
    - Fx: 2D numpy array of the x-component of the vector field.
    - Fy: 2D numpy array of the y-component of the vector field.
    - dx: Grid spacing in the x-direction (in meters).
    - dy: Grid spacing in the y-direction (in meters).
     
    Returns:
    - divergence: 2D numpy array of the divergence of the vector field.
    """
      
    # Fourth-order central difference
    if Fx.shape[0] < 4 or Fx.shape[1] < 4:
        raise ValueError("Input arrays must be at least 4x4 for fourth-order central differences")
    
    dFx_dx= (-Fx[2:-2, 4:] + 8*Fx[2:-2, 3:-1] - 8*Fx[2:-2, 1:-3] + Fx[2:-2, 0:-4])/(12*dx)
    dFy_dy= (-Fy[4:, 2:-2] + 8*Fy[3:-1, 2:-2] - 8*Fy[1:-3, 2:-2] + Fy[0:-4, 2:-2])/(12*dy)

    
    # Ensure the dimensions match
    min_rows = min(dFx_dx.shape[0], dFy_dy.shape[0])
    min_cols = min(dFx_dx.shape[1], dFy_dy.shape[1])
    
    divergence = np.zeros((min_rows, min_cols))
    divergence = dFx_dx[:min_rows, :min_cols] + dFy_dy[:min_rows, :min_cols]
        
    return divergence,dFx_dx[:min_rows, :min_cols], dFy_dy[:min_rows, :min_cols]

def create_grouped_dataset(l4_filename, l4_product):
    """
    Create a NetCDF file grouped by methods, where each method contains variables 
    representing species emissions and their corresponding apriori, estimated, and error values.

    Parameters:
    l4_filename (str): The path where the new NetCDF file will be saved.
    l4_product (dict): A dictionary containing methods as keys and another dictionary 
                       as values representing species and their emission data.
                       Example structure:
                       {
                           'method_1': {
                               'XCO2_Proxy': (apriori, estimated, error),
                               ...
                           },
                           'method_2': {
                               'XCH4_Proxy': (apriori, estimated, error),
                               ...
                           }
                       }
    """
    data_vars = {}

    # Iterate over methods and their respective species results
    for method, species_dict in l4_product.items():
        for sp, values in species_dict.items():
            # Replace forward slashes with underscores in variable names
            safe_method = method.replace('/', '_')
            safe_species = sp.replace('/', '_')
            variable_name = f"{safe_method}_{safe_species}"
            # Assign the tuple (apriori, estimated, error) for each species under the method group
            data_vars[variable_name] = (['stat'], np.array(values))

    # Coordinate labels for the emission data
    stats_coord = ['apriori_emission', 'estimated_emission', 'error_std']

    # Create the xarray dataset
    dataset = xr.Dataset(
        data_vars=data_vars,
        coords={'stat': stats_coord},
        attrs={
            'units': 'kg/s',
            'description': 'Emission estimations grouped by method and species'
        }
    )

    # Write the dataset to a NetCDF file
    try:
        dataset.to_netcdf(l4_filename)
        print(f"Grouped NetCDF file '{l4_filename}' created successfully.")
    finally:
        dataset.close()


def microhh_fit(sgm_filename, sgm_ref_filename, l2_filename, ret_name, apr_name):
    """
    Perform emission estimation by fitting aprifrom scipy.optimize import curve_fitori and L2 data from NetCDF datasets.
    
    Parameters:
    sgm_filename (str): Path to the MicroHH NetCDF file.
    sgm_ref_filename (str): Path to the reference NetCDF file.
    l2_filename (str): Path to the L2 NetCDF file.
    ret_name (str): Name of the retrieval variable.
    apr_name (str): Name of the apriori variable.

    Returns:
    tuple: (apriori emission, retrieved emission, error of the retrieved emission)
    """

    def load_dataset(filename):
        try:
            return xr.open_dataset(filename)
        except Exception as e:
            print(f"Error loading dataset '{filename}': {e}")
            raise

    microhh = load_dataset(sgm_filename)
    apr = load_dataset(sgm_ref_filename)
    l2 = load_dataset(l2_filename)

    apr_emission = microhh[apr_name[1:] + " emission"].values[0]

    apr_data = np.array(apr[apr_name]) / apr_emission
    l2_data = np.array(l2[ret_name]) / 1e6

    # Background correction
    apr_data -= np.nanmedian(apr_data)
    l2_data -= np.nanmedian(l2_data)

    # Plume detection
    idx = apr_data > 0
    apr_data_flat = apr_data[idx].flatten()
    l2_data_flat = l2_data[idx].flatten()

    # Model for fitting
    def model(x, alpha):
        return alpha * x

    try:
        params, covariance = curve_fit(model, apr_data_flat, l2_data_flat, p0=[2.0])
        ret_emission = params[0]
        ret_emission_error = np.sqrt(covariance[0, 0])
        print(f"Fitted scaling factor (alpha, {ret_name}): {ret_emission:.3f} ± {ret_emission_error:.3f}")
    except Exception as e:
        print(f"Error in curve fitting: {e}")
        raise
    finally:
        apr.close()
        l2.close()
        microhh.close()

    return apr_emission, ret_emission, ret_emission_error