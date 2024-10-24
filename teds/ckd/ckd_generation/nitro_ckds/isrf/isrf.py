"""
Instrument Spectral Response Function
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import csv
from scipy.interpolate import CubicSpline
from teds.ckd.ckd_generation.nitro_ckds.spectral.wavelength import mirror_on_first_row

def generate(ncc): # Get the data and combine
    cfg = ncc.cfg
    filename = cfg['paths']['dir_external'] + 'Nominal_ISRF_dataset_17102024.xlsx'
    df_sim = pd.read_excel(filename, "SimulationLogic", skiprows = 3, nrows = 181)
    df_sim.columns = ['ix', 'wavelength', 'actpos']
    df_sim['wavelength'] *= 1000 # convert wavelength to nm
    df_wl = pd.read_excel(filename, "WavelengthDifferenceSampling", skiprows = 4, nrows = 1029, usecols= "A:FY")
    df_isrf = pd.read_excel(filename, "ISRF", skiprows = 5, nrows = 1029, usecols= "H:GF")

    sim_indices = df_sim.iloc[:,0]
    wavelengths = np.unique(df_sim.iloc[:,1])
    actpos = np.unique(df_sim.iloc[:,2])
    
    # fit supergaussian(x, x0, w, n)
    popt = np.zeros((len(sim_indices), 3))
    significant_indices = np.zeros((len(sim_indices)))
    for sim_ix in sim_indices:
        wldif = np.array(df_wl[sim_ix])  # nm
        isrf = np.array(df_isrf[sim_ix]) # AU
        
        p0 = [0, 0.6, 10]
        bmin = [-0.3, 0.05, 1]
        bmax = [0.3, 5, 100]
        popt[sim_ix], pcov = curve_fit(supergaussian, wldif, isrf, p0, bounds = (bmin, bmax))
        
        thr = 0.001
        mask = isrf > thr
        significant_indices[sim_ix] = np.max(np.abs(wldif[mask]))

        if False: # plot to check fit accuracy
            fig, ax = plt.subplots()
            fit = supergaussian(wldif, *popt[sim_ix])
            ax.plot(wldif, isrf)
            ax.plot(wldif, fit)
            fig.savefig("./ckd/ckd_generation/nitro_ckds/isrf/test.png")
    
    max_wl_diff_to_compute_isrf_at = np.max(significant_indices) # 0.61

    df_sim["x0"] = popt[:, 0]
    df_sim["w"] = popt[:, 1]
    df_sim["n"] = popt[:, 2]
    
    if False: # fit to check smoothness of parameters
        fig, ax = plt.subplots(2, 3, figsize = (12,8))
        colors = plt.cm.viridis(np.linspace(0, 1, len(wavelengths)))
        for i, wl in enumerate(wavelengths):
            mask = df_sim.iloc[:,1] == wl
            ax[1, 0].plot(df_sim.iloc[:,2][mask], df_sim["x0"][mask], c = colors[i])
            ax[1, 1].plot(df_sim.iloc[:,2][mask], df_sim["w"][mask], c = colors[i])
            ax[1, 2].plot(df_sim.iloc[:,2][mask], df_sim["n"][mask], c = colors[i])

        colors = plt.cm.viridis(np.linspace(0, 1, len(actpos)))
        for i, ap in enumerate(actpos):
            mask = df_sim.iloc[:,2] == ap
            ax[0, 0].plot(df_sim.iloc[:,1][mask], df_sim["x0"][mask], c = colors[i])
            ax[0, 1].plot(df_sim.iloc[:,1][mask], df_sim["w"][mask], c = colors[i])
            ax[0, 2].plot(df_sim.iloc[:,1][mask], df_sim["n"][mask], c = colors[i])

        ax[0, 0].set_title("x0")
        ax[0, 1].set_title("w")
        ax[0, 2].set_title("n")
        [a.set_xlabel("wavelength") for a in ax[0,:]]
        [a.set_xlabel("ACT position") for a in ax[1,:]]
        fig.suptitle("Tango Nitro - ISRF supergaussian fit parameters")
        plt.tight_layout()
        fig.savefig("./ckd/ckd_generation/nitro_ckds/isrf/test.png")

    # parameters are smooth so can be interpolated for ACT and wavelength dimensions of input spectra
    params_intermediate = np.zeros((len(actpos), cfg["dimensions"]["lbl_samples"], 3))
    params = np.zeros((cfg["dimensions"]["across_track"], cfg["dimensions"]["lbl_samples"], 3))
    # interpolate wavelengths
    lbl_wavelengths = np.linspace(cfg["lbl_min_wl"], cfg["lbl_max_wl"], cfg["dimensions"]["lbl_samples"])
    for i, ap in enumerate(actpos):
        mask = df_sim.iloc[:,2] == ap
        df_filt = df_sim[mask]
        for k, param in enumerate(["x0", "w", "n"]):
            wl = np.array(df_filt["wavelength"])
            p = np.array(df_filt[param])
            spl = CubicSpline(wl, p)
            params_intermediate[i, :, k] = spl(lbl_wavelengths)
    
    # interpolate act
    new_actpos = np.linspace(actpos[0], actpos[-1], cfg["dimensions"]["across_track"])
    for i, wl in enumerate(lbl_wavelengths):
        mask = df_sim.iloc[:,1] == wl
        for k, param in enumerate(["x0", "w", "n"]):
            spl = CubicSpline(actpos, params_intermediate[:, i, k])
            params[:,i,k] = spl(new_actpos)

    # Add fit parameters as variables
    dims = ['across_track', 'lbl_samples']

    supergaussian_str =  "Supergaussian = exp(-(|x-x0| / w)^2n)"
    attr = {"long_name": "x0, center position", "comment" : supergaussian_str}
    ncc.create_var('x0', dims, params[:,:,0], attr, 'f8')

    attr = {"long_name": "w, width", "comment" : supergaussian_str}
    ncc.create_var('w', dims, params[:,:,1], attr, 'f8')

    attr = {"long_name": "n, rank", "comment" : supergaussian_str}
    ncc.create_var('n', dims, params[:,:,2], attr, 'f8')

def supergaussian(x, x0, w, n):
    return np.exp(-(np.abs(x-x0)/w)**(2*n))