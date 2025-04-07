# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
from typing import Any
import configparser
import numpy as np
import os
import xarray as xr


def xy2latlon(origin: tuple[float, float],
              x1: np.ndarray,
              y1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert from x-y to lat-lon.

    Parameters
    ----------
    origin
        Origin in lat-lon coordinates
    x1
        x coordinate in km
    y1
        y coordinate in km

    Returns
    -------
    lat
        Latitude
    lon
        Longitude

    """
    phi0 = np.deg2rad(origin[0])
    ld0 = np.deg2rad(origin[1])
    s_p0 = np.sin(phi0)
    c_p0 = np.cos(phi0)
    earth_radius = 6371000.0
    x, y = x1 / earth_radius, y1 / earth_radius
    p = np.maximum(np.sqrt(x**2 + y**2), 1.0e-10)
    c = 2.0 * np.arcsin(p / 2.0)
    s_c = np.sin(c)
    c_c = np.cos(c)
    phi = np.arcsin(c_c * s_p0 + y * s_c * c_p0 / p)
    ld = ld0 + np.arctan2(
        x * s_c, (p * c_p0 * c_c - y * s_p0 * s_c))
    lat = np.rad2deg(phi)
    lon = np.rad2deg(ld)
    if isinstance(lat, np.ndarray):
        lat[lat > 90.0] -= 180.0
        lat[lat < -90.0] += 180.0
        lon[lon > 180.0] -= 360.0
        lon[lon < -180.0] += 360.0
    else:
        if abs(lat) > 90.0:
            if lat > 0:
                lat = lat - 180.0
            else:
                lat = lat + 180.0
        if abs(lon) > 180.0:
            if lon > 0:
                lon = lon - 360.0
            else:
                lon = lon + 360.0
    return lat, lon


def microhh_calculate_znodes(z_centers: np.ndarray) -> np.ndarray:
    """Calculate znodes from z center heights.

    Parameters
    ----------
    z_centers
        Array of z center heights.

    Returns
    -------
        Array of z node heights.

    """
    znodes = np.zeros(len(z_centers) + 1)
    # Midpoints between z levels
    znodes[1:-1] = 0.5 * (z_centers[:-1] + z_centers[1:])
    # Extrapolate lower bound
    znodes[0] = z_centers[0] - (z_centers[1] - z_centers[0]) / 2
    # Extrapolate upper bound
    znodes[-1] = z_centers[-1] + (z_centers[-1] - z_centers[-2]) / 2
    return znodes


def microhh_read_par(file_list: list[str],
                     par_list: list[str],
                     group: str | None = None) -> dict:
    """Read parameters from a list of NetCDF files.

    Parameters
    ----------
    file_list
        List of file paths.
    par_list
        List of parameter names to read.
    group
        Group within the NetCDF file. Defaults to None.

    Returns
    -------
        Dictionary of parameter values.

    """
    result = {}

    if len(file_list) == 1:
        with xr.open_dataset(file_list[0], group=group) as ds:
            for par in par_list:
                result[par] = ds[par].values
    else:
        for file, par in zip(file_list, par_list):
            with xr.open_dataset(file, group=group) as ds:
                result[par] = ds[par].values
    return result


def microhh_create_netcdf(
    filename: str,
    time_data: np.ndarray,
    x_data: np.ndarray,
    y_data: np.ndarray,
    z_data: np.ndarray,
    znodes_data: np.ndarray,
    lon_data: np.ndarray,
    lat_data: np.ndarray,
    data_list: list[np.ndarray],
    data_names: list[str],
    data_attrs: list[dict],
    title: str = "Dataset"
) -> None:
    """Create a NetCDF file with multiple fields and metadata.

    Parameters
    ----------
    filename
        Path to save the NetCDF file.
    time_data
        Time coordinate data.
    x_data
        X-coordinate data (cell centers in x-direction).
    y_data
        Y-coordinate data (cell centers in y-direction).
    z_data
        Z-coordinate data (cell centers in z-direction).
    znodes_data
        Z-nodes data (nodes in z-direction).
    lon_data
        Longitude data (2D array for y, x).
    lat_data
        Latitude data (2D array for y, x).
    data_list
        List of 4D data arrays (time, z, y, x).
    data_names
        List of names corresponding to the data arrays.
    data_attrs
        List of attribute dictionaries for each data array.
    title
        Title for the dataset. Default is "Dataset".

    """
    try:
        # Validate input dimensions
        if lon_data.shape != lat_data.shape:
            raise ValueError(
                "Longitude and latitude data must have the same shape.")
        for data in data_list:
            if data.shape != (
                    len(time_data), len(z_data), len(y_data), len(x_data)):
                raise ValueError(
                    "Each data array must have dimensions (time, z, y, x).")
        if len(data_list) != len(data_names) or len(data_list) != len(
                data_attrs):
            raise ValueError(
                "The number of data arrays, names, and attributes must match.")

        # Create the dataset
        variables: dict[Any, Any] = {
            "longitude": (("y", "x"), lon_data, {
                "long_name": "cell centers",
                "units": "degree",
                "grid_type": "3DCoRectMesh"}),
            "latitude": (("y", "x"), lat_data, {
                "long_name": "cell centers",
                "units": "degree",
                "grid_type": "3DCoRectMesh"}),
        }

        # Add each data array to the dataset with its attributes
        for data, name, attrs in zip(data_list, data_names, data_attrs):
            variables[name] = (("time", "z", "y", "x"), data, attrs)

        ds = xr.Dataset(
            variables,
            coords={
                "x": (("x"), x_data, {
                    "long_name": "cell centers in x-direction",
                    "units": "m",
                    "grid_type": "3DCoRectMesh"}),
                "y": (("y"), y_data, {
                    "long_name": "cell centers in y-direction",
                    "units": "m",
                    "grid_type": "3DCoRectMesh"}),
                "time": ("time", time_data, {
                    "units": "hours since simulation start",
                    "calendar": "standard"}),
                "z": ("z", z_data, {
                    "long_name": "cell centers in z-direction",
                    "units": "m",
                    "grid_type": "3DCoRectMesh"}),
                "znodes": ("znodes", znodes_data, {
                    "long_name": "Nodes in z-direction",
                    "units": "m",
                    "grid_type": "3DCoRectMesh"})
            },
            attrs={"title": title}
        )

        # Save to NetCDF
        ds.to_netcdf(filename)
        print(f"NetCDF file created successfully: {filename}")

    except Exception as e:
        print(f"Error creating NetCDF file: {e}")


def microhh_create_tedsnetcdf(config: dict,
                              return_data: str = "netcdf") -> None | dict:
    """Create a NetCDF file from MicroHH simulation data using the
    provided YAML configuration.

    Parameters
    ----------
    config
        Configuration dictionary parsed from YAML.
    return_data
        Whether to return data or create a NetCDF file. Default is
        "netcdf".

    Returns:
        Returns data dictionary if return_data is not "netcdf".

    """

    # Extract configuration values
    # Path to the MicroHH simulation data
    path = config["simulation"]["input_directory"]
    # Path and name of the output NetCDF file
    output_file = config["simulation"]["output_file"]
    # Rotation angle in degrees
    rotation_angle = config["simulation"]["rotation_angle"]
    # Latitude of the source
    lat0 = config["simulation"]["source_location"]["lat0"]
    # Longitude of the source
    lon0 = config["simulation"]["source_location"]["lon0"]
    # Default filename for auxiliary data
    default_filename = config["simulation"]["default_filename"]
    # Filename for plume.ini
    ini_filename = config["simulation"]["ini_filename"]
    # Ensure all names are uppercase
    data_names = np.array(
        [name.upper() for name in config["simulation"]["data_names"]])
    # Set emissions data
    data_set_emissons = config["simulation"]["data_set_emissions"]

    # --- Define Constants ---
    avogadro_number = 6.02214076e23  # molecules/mol

    # Molecular masses (kg/mol)
    molecular_masses = {
        "air": 0.0289647,
        "no": 0.030006,
        "no2": 0.0460055,
        "o3": 0.048,
        "co2": 0.044009,
        "ch4": 0.016042
    }

    # --- Read MicroHH Data ---
    gases = {}
    wind = {}

    for name in data_names:
        key = name.lower()  # Convert name to lowercase
        file_path = os.path.join(path, f"{key}.nc")  # Construct the file path

        if key in ["co2", "no", "no2", "o3"]:  # Check if it's a gas
            if os.path.exists(file_path):  # Ensure the file exists
                gases[key] = microhh_read_par([file_path], [key])[key]
            else:
                print(f"Warning: File for gas '{key}' not found at "
                      f"{file_path}. Skipping.")
        elif key in ["u", "v", "w"]:  # Check if it's a wind component
            if os.path.exists(file_path):  # Ensure the file exists
                wind[key] = microhh_read_par([file_path], [key])[key]
            else:
                print(f"Warning: File for wind component '{key}' not found "
                      f"at {file_path}. Skipping.")
        else:
            # CH4 is calculated later, so skip the warning for it
            if key != "ch4":
                print(f"Warning: '{name}' is not recognized as a gas or wind "
                      "component. Skipping.")

    # --- Read Grid and Auxiliary Data ---
    # Use "no.nc" as a reference file
    file_list_dim = [os.path.join(path, "no.nc")]
    par_list_dim = ["x", "y", "z", "time"]
    dim = microhh_read_par(file_list_dim, par_list_dim)
    dim["dx"] = np.abs(dim["x"][1] - dim["x"][0])
    dim["dy"] = np.abs(dim["y"][1] - dim["y"][0])
    dim["dz"] = np.abs(dim["z"][1] - dim["z"][0])

    file_list_aux = [os.path.join(path, default_filename)]
    par_list_aux = ["rhoref"]
    aux = microhh_read_par(file_list_aux, par_list_aux, group="thermo")
    aux["rho_kg"] = aux["rhoref"]
    # Convert to mol/m³
    aux["rho_mol"] = aux["rho_kg"] / molecular_masses["air"]

    # --- Read Source Information from plume.ini ---
    plume_ini_path = os.path.join(path, ini_filename)
    config_parser = configparser.ConfigParser()
    config_parser.read(plume_ini_path)

    source_info: dict[Any, Any] = {}
    if "source" in config_parser:
        source_info["sourcelist"] = config_parser.get(
            "source", "sourcelist").split(",")
        source_info["source_x0"] = [
            float(x) for
            x in config_parser.get("source", "source_x0").split(",")]
        source_info["source_y0"] = [
            float(x) for
            x in config_parser.get("source", "source_y0").split(",")]
        source_info["source_z0"] = [
            float(x) for
            x in config_parser.get("source", "source_z0").split(",")]
        source_info["strength"] = [
            float(x) for
            x in config_parser.get("source", "strength").split(",")]
    else:
        raise ValueError(
            "The [source] group is missing in the plume.ini file.")

    # --- Process Gases and Emissions ---
    sourcelist_array = np.array(source_info["sourcelist"])
    gas_strengths_kgps = {}

    for name in data_names:
        key = name.lower()
        if key in ["co2", "no", "no2"]:
            if key in sourcelist_array:
                idx = np.where(sourcelist_array == key)[0][0]
                strength_kmolps = source_info["strength"][idx]
                molecular_mass = molecular_masses[key]
                gas_strengths_kgps[key] = (
                    strength_kmolps * 1000 * molecular_mass)  # Convert to kg/s

                # Rescale the gas concentration to the emission
                # strength from the config file.
                sel = np.where(data_names == name)[0][0]
                if data_set_emissons[sel] != "None":
                    gases[key] = (gases[key]
                                  / gas_strengths_kgps[key]
                                  * data_set_emissons[sel])
                    gas_strengths_kgps[key] = data_set_emissons[sel]
            else:
                print(
                    f"Warning: Gas '{key}' not found in sourcelist. Skipping.")

    if "CH4" in data_names:
        if "CO2" not in data_names:
            raise ValueError("Error: CH4 can only be calculated when CO2 is "
                             "selected in data_names.")

        # Find the index of CH4 in data_names and select the
        # corresponding emission from data_set_emissions.
        # Get the index of CH4 in data_names
        ch4_index = np.where(data_names == "CH4")[0][0]
        # Select the emission value for CH4
        gas_strengths_kgps["ch4"] = data_set_emissons[ch4_index]
        gases["ch4"] = (gases["co2"]
                        / gas_strengths_kgps["co2"]
                        * gas_strengths_kgps["ch4"])

    for name in data_names:
        key = name.lower()
        if key in gases:
            gases[key] *= (aux["rho_mol"][:, np.newaxis, np.newaxis]
                           * avogadro_number * dim["dz"])

    # --- Rotate and Transform Coordinates ---
    x_vec = (dim["x"] - source_info["source_x0"][0])
    y_vec = (dim["y"] - source_info["source_y0"][0])
    x, y = np.meshgrid(x_vec, y_vec)

    theta = np.radians(rotation_angle)
    rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                                [np.sin(theta), np.cos(theta)]])
    coordinates = np.vstack([x.flatten(), y.flatten()])
    rotated_coordinates = np.dot(rotation_matrix, coordinates)
    x_rot = rotated_coordinates[0].reshape(x.shape)
    y_rot = rotated_coordinates[1].reshape(y.shape)

    lat, lon = xy2latlon((lat0, lon0), x_rot, y_rot)

    # --- Prepare Data for NetCDF Creation ---
    data_list = []
    data_attrs = []
    source_location = (source_info["source_y0"][0], lat0, lon0)

    for name in data_names:
        key = name.lower()
        if key in gases:
            data_list.append(gases[key])
            data_attrs.append({
                "units": "molecules/m^2",
                "long_name": f"{name} column density",
                "source": source_location,
                "emission_in_kgps": gas_strengths_kgps[key]
            })
        elif key in wind:
            data_list.append(wind[key])
            data_attrs.append({
                "units": "m/s",
                "long_name": f"{name} wind component"
            })
        else:
            print(f"Warning: '{name}' not found in gases or wind "
                  "dictionaries. Skipping.")

    if return_data == "netcdf":

        # --- Create NetCDF File ---
        microhh_create_netcdf(
            output_file,
            dim["time"] / 3600,  # Convert time to hours
            x_vec,
            y_vec,
            dim["z"],
            microhh_calculate_znodes(dim["z"]),
            lon,
            lat,
            data_list,
            list(data_names),
            data_attrs
        )
        return None
    else:
        return {
            "time": dim["time"] / 3600,
            "x": x_vec,
            "y": y_vec,
            "z": dim["z"],
            "znodes": microhh_calculate_znodes(dim["z"]),
            "lon": lon,
            "lat": lat,
            "gases": gases,
            "gas_strengths_kgps": gas_strengths_kgps,
            "wind": wind,
            "data_list": data_list,
            "data_names": data_names,
            "data_attrs": data_attrs
        }
