#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 12:54:22 2023

@author: manugv
"""


def emissions(filename):
    """Compute emissions based on different methods
    `
    Parameters
    ----------
    filename : String
        File name of the yaml config file

    Examples
    --------
    from L4 import Compute
    Compute.emissions(filename)

    """

    import ModulePreProcessing as preprocess
    # read parameters from yaml file
    params = preprocess.readyamlfile(filename)

    if params["method"] == "LSQ":

        import ModuleLSQ as lsq
        from numpy import median

        print("Estimating emissions using Least squares estimate (LSQ)")
        # here 405ppm is the background.
        data = preprocess.getendtoendsimdata(params)
        print("   Data read successfully")
        back = median(data.actual_column)
        emission, precision = lsq.emissionprecision(data.actual_column - back,
                                                    data.lvl2data.data - back,
                                                    data.lvl2precision,
                                                    params["plumethreshold"])
        emission = emission * params["emission"]
        print(f"LSQ estimated emission is{emission: .2f} kg/s")
        print(f"LSQ estimated level-4 precision is{precision: .2e} kg/s")

    elif params["method"] == "CFM":

        import ModuleCFM as cfm

        print("Estimating emissions using Cross-sectional Flux Method (CFM)")
        data = preprocess.getendtoendsimdata(params)
        interp_u, interp_v, microhhdata = preprocess.microhhvelocityinterp(params)
        co2_conc_kg = data.lvl2data*data.ppm_to_kg_gas
        print("   Data read successfully")
        massflux, emission = cfm.get_massflux(co2_conc_kg, data.grid,
                                              params["plumethreshold"],
                                              interp_u, interp_v)
        if emission in None:
            print("CFM failed and emission is not estimated")
        else:
            print(f"CFM estimated emission is{emission: .2f} kg/s")
