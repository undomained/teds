#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 12:54:22 2023

@author: manugv
"""

from numpy import median, sum
from .ModulePreProcessing import getendtoendsimdata, microhhvelocityinterp
from .ModuleLevel2 import read_level2_product
from ..lib.libINV import lsq_fit
from .ModuleCFM import get_massflux
import numpy as np
import sys

def level2_to_level4_processor(config):
    """Compute emissions from level 2 product based on different methods
    Parameters
    ----------
    config : list 
        File name of the yaml config fileof configuration parameter
        ['method']  currently 'LSQ' and 'CFM' implemented
        []
    """
    if config["method"] == "LSQ":
        #Assume the transport is known, so the measurement can be written as 
        # ymeas = K x
        #with the measurement ymeas, the Jacobian K, and the emissions x

        print("Estimating emissions using Least squares estimate (LSQ)")

        model  = getendtoendsimdata(config["sgm_input"], config["gas"], config["lat_lon_src"])
        l2data = read_level2_product(config["l2_input"], config["gas"])

        print("   Data read successfully")

        # define plume mask
        background = median(l2data.Xgas[:].data)
        ymask = model.Xgas.flatten() < (1.+config["lsq_plumethreshold"])*background
        
        ymeas = l2data.Xgas.flatten()[ymask == False] 
        
        if(config['avg_kernel']):

            #calculate effective total model columns using the column avergaing kernel of the l2 product 
            if(config["gas"] == 'co2'):
                conv_fact = 1.E6  #ppm
            if(config["gas"] == 'ch4'):
                conv_fact = 1.E9  #ppb
                
            Acol_Xgas = conv_fact * (sum(model.dcol_gas * l2data.avg_kernel, axis=2))/model.column_air
            Kmat  = Acol_Xgas.flatten()[ymask == False]/config["emission"]
            
        else:
            Kmat  = model.Xgas.flatten()[ymask == False]/config["emission"]
            
        yprec = l2data.Xgas_precision.flatten()[ymask == False]
        Sy    = np.diag(yprec**2)
        
        flux, Sflux = lsq_fit(ymeas, Kmat, Sy)
        flux_prec = np.sqrt(Sflux)

        print(f"LSQ estimated emission is{flux: .2f} kg/s")
        print(f"LSQ estimated level-4 precision is{flux_prec: .3f} kg/s")
        
    elif config["method"] == "CFM":

        print("Estimating emissions using Cross-sectional Flux Method (CFM)")
        data = getendtoendsimdata(config)
        interp_u, interp_v, microhhdata = microhhvelocityinterp(config)
        co2_conc_kg = data.lvl2data*data.ppm_to_kg_gas
        print("   Data read successfully")
        massflux, emission = get_massflux(co2_conc_kg, data.grid,
                                          config["plumethreshold"],
                                          interp_u, interp_v)
        if emission in None:
            print("CFM failed and emission is not estimated")
        else:
            print(f"CFM estimated emission is{emission: .2f} kg/s")
