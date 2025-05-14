#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
from copy import deepcopy

import netCDF4 as nc
import numpy as np

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def _get_l1b(l1b_filename: str) -> dict:
    """
    Reads and extracts relevant L1B data from a NetCDF file. The function processes the input file,
    retrieving and organizing various observation and geolocation data, including radiance,
    wavelength, noise, solar angles, sensor angles, and spatial coordinates. Additionally,
    it derives a mask to identify valid data points within the radiance dataset.

    Parameters:
        l1b_filename: str
            The path to the L1B NetCDF file to be read.

    Returns:
        dict
            A dictionary containing the extracted L1B data. Keys include:
            'wavelength', 'radiance', 'noise', 'sza', 'saa', 'vza', 'vaa',
            'latitude', 'longitude', and 'mask', representing respective
            data arrays.
    """

    # getting l1b data from file
    logger.debug(f"Reading L1B data from {l1b_filename}")
    nc_l1b = nc.Dataset(l1b_filename, mode="r")

    l1b_data = {
        "wavelength": deepcopy(nc_l1b["observation_data/wavelength"][:]),
        "radiance": deepcopy(nc_l1b["observation_data/radiance"][:]),
        "noise": deepcopy(nc_l1b["observation_data/radiance_stdev"][:]),
        "sza": deepcopy(nc_l1b["geolocation_data/solar_zenith"][:]),
        "saa": deepcopy(nc_l1b["geolocation_data/solar_azimuth"][:]),
        "vza": deepcopy(nc_l1b["geolocation_data/sensor_zenith"][:]),
        "vaa": deepcopy(nc_l1b["geolocation_data/sensor_azimuth"][:]),
        "latitude": deepcopy(nc_l1b["geolocation_data/latitude"][:]),
        "longitude": deepcopy(nc_l1b["geolocation_data/longitude"][:]),
    }

    # Extract mask
    nc_var = nc_l1b["observation_data/radiance"]
    l1b_data["mask"] = ~nc_var[:].mask
    if not nc_var[:].mask.shape:
        l1b_data["mask"] = np.full(nc_var[:].shape, True)

    return l1b_data


def _sim_modified_output(filename: str, l1b_output: dict) -> None:
    """
    Creates a NetCDF output file containing geolocation and observation data based
    on the given inputs.

    This function generates a NetCDF file based on provided data arrays representing
    geolocation and observation data. It defines dimensions, groups, and variables
    within the file to store latitude, longitude, solar and sensor angles, as well
    as spectral radiance and noise levels.

    Arguments:
        filename (str): The path to the NetCDF file to be created.
        l1b_output (dict): A dictionary containing geolocation and observation
            data arrays. Expected keys are:
            - 'latitude': 2D numpy array with latitude values.
            - 'longitude': 2D numpy array with longitude values.
            - 'sza': 2D numpy array with solar zenith angles.
            - 'saa': 2D numpy array with solar azimuth angles.
            - 'vza': 2D numpy array with sensor zenith angles.
            - 'vaa': 2D numpy array with sensor azimuth angles.
            - 'wavelength': 2D numpy array with wavelength values.
            - 'radiance': 3D numpy array with spectral radiance values.
            - 'noise': 3D numpy array with radiance noise values.

    Raises:
        IOError: If there is an issue creating or writing to the file.

    Returns:
        None
    """
    output = nc.Dataset(filename, mode="w")
    output.title = "Tango Carbon level 1B data"

    nalt, nact, nwav = l1b_output["radiance"].shape
    output.createDimension("wavelength", nwav)
    output.createDimension("across_track_sample", nact)
    output.createDimension("along_track_sample", nalt)

    nc_grp = output.createGroup("geolocation_data")
    _dim2 = ("along_track_sample", "across_track_sample")

    nc_var = nc_grp.createVariable("latitude", "f4", _dim2, fill_value=-32767.0)
    nc_var.long_name = "latitude at bin locations"
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = "degrees_north"
    nc_var[:] = l1b_output["latitude"]

    nc_var = nc_grp.createVariable("longitude", "f4", _dim2, fill_value=-32767.0)
    nc_var.long_name = "longitude at bin locations"
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = "degrees_east"
    nc_var[:] = l1b_output["latitude"]

    nc_var = nc_grp.createVariable("solar_zenith", "f4", _dim2, fill_value=-32767.0)
    nc_var.long_name = "solar zenith angle at bin locations"
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = "degrees"
    nc_var[:] = l1b_output["sza"]

    nc_var = nc_grp.createVariable("solar_azimuth", "f4", _dim2, fill_value=-32767.0)
    nc_var.long_name = "solar azimuth angle at bin locations"
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = "degrees"
    nc_var[:] = l1b_output["saa"]

    nc_var = nc_grp.createVariable("sensor_zenith", "f4", _dim2, fill_value=-32767.0)
    nc_var.long_name = "sensor zenith angle at bin locations"
    nc_var.valid_min = -90.0
    nc_var.valid_max = 90.0
    nc_var.units = "degrees"
    nc_var[:] = l1b_output["vza"]

    nc_var = nc_grp.createVariable("sensor_azimuth", "f4", _dim2, fill_value=-32767.0)
    nc_var.long_name = "sensor azimuth angle at bin locations"
    nc_var.valid_min = -180.0
    nc_var.valid_max = 180.0
    nc_var.units = "degrees"
    nc_var[:] = l1b_output["vaa"]

    nc_grp = output.createGroup("observation_data")
    _dim3 = ("along_track_sample", "across_track_sample", "wavelength")
    _dim2 = ("across_track_sample", "wavelength")
    nc_var = nc_grp.createVariable("wavelength", "f8", _dim2, fill_value=-32767.0)
    nc_var.long_name = "wavelength"
    nc_var.valid_min = 0
    nc_var.valid_max = 8.0e3
    nc_var.units = "nm"
    nc_var[:] = l1b_output["wavelength"]

    nc_var = nc_grp.createVariable("radiance", "f8", _dim3, fill_value=-32767.0)
    nc_var.long_name = "spectral radiance"
    nc_var.valid_min = 0
    nc_var.valid_max = 1.0e28
    nc_var.units = "photons / (sr nm m2 s)"
    nc_var[:] = l1b_output["radiance"]

    nc_var = nc_grp.createVariable("radiance_stdev", "f8", _dim3, fill_value=-32767.0)
    nc_var.long_name = "relative radiance noise (1-sigma)"
    nc_var.units = "photons/(sr nm m2 s)"
    nc_var.valid_min = 0.0
    nc_var.valid_max = 1e24
    nc_var[:] = l1b_output["noise"]

    output.close()


def add_radiance_offset(filename_in: str, filename_out: str, rad_offset: float) -> None:
    logger.debug(f"Adding radiance offset of {rad_offset} to {filename_in}")
    """
    Adds a radiometric offset to radiance data and saves the modified data to a
    NetCDF file. The function processes input data by applying a scaled radiometric
    offset to the radiance field, while duplicating other data fields consistently
    in the output structure.

    Args:
        filename_in (str): The path to the input file containing radiance
        and associated data.
        filename_out (str): The path where the output file with the modified
        radiance data will be saved.
        rad_offset (float): The scaling factor used to calculate the radiometric
        offset applied to the radiance data.

    Raises:
        None

    Returns:
        None
    """
    l1b: dict = _get_l1b(filename_in)
    nalt, nact, nwav = l1b["radiance"].shape

    l1b_scaled: dict = {
        "sza": np.empty([nalt, nact]),
        "saa": np.empty([nalt, nact]),
        "vza": np.empty([nalt, nact]),
        "vaa": np.empty([nalt, nact]),
        "latitude": np.empty([nalt, nact]),
        "longitude": np.empty([nalt, nact]),
        "noise": np.empty([nalt, nact, nwav]),
        "radiance": np.empty([nalt, nact, nwav]),
        "wavelength": np.empty([nwav]),
    }

    for ialt in range(nalt):
        for iact in range(nact):
            offset = rad_offset * np.max(l1b["radiance"][0, iact, :])
            l1b_scaled["radiance"][ialt, iact, :] = l1b["radiance"][0, iact, :] + offset
        l1b_scaled["sza"][ialt, :] = l1b["sza"][0, :]
        l1b_scaled["saa"][ialt, :] = l1b["saa"][0, :]
        l1b_scaled["vza"][ialt, :] = l1b["vza"][0, :]
        l1b_scaled["vaa"][ialt, :] = l1b["vaa"][0, :]
        l1b_scaled["latitude"][ialt, :] = l1b["latitude"][0, :]
        l1b_scaled["longitude"][ialt, :] = l1b["longitude"][0, :]
        l1b_scaled["noise"][ialt, :, :] = l1b["noise"][0, :, :]
    l1b_scaled["wavelength"] = l1b["wavelength"]

    # output to netcdf file
    _sim_modified_output(filename_out, l1b_scaled)

    logger.info("Radiometric offset added successfully ")
