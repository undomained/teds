#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import netCDF4 as nc


def main():
    """
    Creates a NetCDF4 dataset and writes geolocation and observation data.

    This function generates a new NetCDF4 file named `sample_l1b.nc` with predefined
    dimensions and groups for geolocation and observation data. It initializes
    the data structures for geolocation variables and observations, and populates
    them with random or calculated values.
    """
    nalt, nact, nwav = 2, 3, 4
    with nc.Dataset("sample_l1b.nc", "w", format="NETCDF4") as ds:
        ds.createDimension("wavelength", nwav)
        ds.createDimension("across_track_sample", nact)
        ds.createDimension("along_track_sample", nalt)

        geo = ds.createGroup("geolocation_data")
        obs = ds.createGroup("observation_data")

        d2 = ("along_track_sample", "across_track_sample")
        d3 = ("along_track_sample", "across_track_sample", "wavelength")
        d2_wave = ("across_track_sample", "wavelength")

        for name in [
            "latitude",
            "longitude",
            "solar_zenith",
            "solar_azimuth",
            "sensor_zenith",
            "sensor_azimuth",
        ]:
            geo.createVariable(name, "f4", d2)[:] = np.random.uniform(-90, 90, (nalt, nact))

        obs.createVariable("wavelength", "f8", d2_wave)[:] = (
            np.linspace(1590, 1675, nwav).reshape(1, -1).repeat(nact, axis=0)
        )
        obs.createVariable("radiance", "f8", d3)[:] = np.random.rand(nalt, nact, nwav) * 1e5
        obs.createVariable("radiance_stdev", "f8", d3)[:] = np.random.rand(nalt, nact, nwav) * 1e2


if __name__ == "__main__":
    main()
