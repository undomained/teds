#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 12:54:22 2023

@author: manugv
"""

from numpy import median, sum
from .ModulePreProcessing import getendtoendsimdata, readyamlfile, microhhvelocityinterp
from .ModuleLSQ import emissionprecision
from .ModuleCFM import get_massflux


def emissions(filename):
    """Compute emissions based on different methods
    Parameters
    ----------
    filename : String
        File name of the yaml config file

    Examples
    --------
    from L4 import Compute
    Compute.emissions(filename)

    """

    # read parameters from yaml file
    params = readyamlfile(filename)

    if params["method"] == "LSQ":

        print("Estimating emissions using Least squares estimate (LSQ)")
        # here 405ppm is the background.
        data = getendtoendsimdata(params)
        print("   Data read successfully")
        # without averaging kernel
        back = median(data.actual_column)
        emission, precision = emissionprecision(data.actual_column - back,
                                                data.lvl2data - back,
                                                data.lvl2precision,
                                                params["plumethreshold"])
        emission = emission * params["emission"]
        print(f"LSQ estimated emission is{emission: .2f} kg/s")
        print(f"LSQ estimated level-4 precision is{precision: .2e} kg/s")
        # with averaging kernel correction
        corrected_column = 1e6 * (sum(data.dactual_column * data.avg_kernel, axis=2))/data.actual_column_air
        emission, precision = emissionprecision(corrected_column - back,
                                                data.lvl2data - back,
                                                data.lvl2precision,
                                                params["plumethreshold"])
        emission = emission * params["emission"]
        print(f"LSQ estimated emission after averaging kernel correction is{emission: .2f} kg/s")
        print(f"LSQ estimated level-4 precision after averaging kernel correction is{precision: .2e} kg/s")

    elif params["method"] == "CFM":

        print("Estimating emissions using Cross-sectional Flux Method (CFM)")
        data = getendtoendsimdata(params)
        interp_u, interp_v, microhhdata = microhhvelocityinterp(params)
        co2_conc_kg = data.lvl2data*data.ppm_to_kg_gas
        print("   Data read successfully")
        massflux, emission = get_massflux(co2_conc_kg, data.grid,
                                          params["plumethreshold"],
                                          interp_u, interp_v)
        if emission in None:
            print("CFM failed and emission is not estimated")
        else:
            print(f"CFM estimated emission is{emission: .2f} kg/s")
