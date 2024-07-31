import numpy as np
import os
import sys
import yaml
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
import cartopy.crs as crs
import matplotlib as mpl
from scipy.stats import linregress
from pathlib import Path

from teds import log

def read_file(file):
    # read in netcdf file, return as dict
    # descends into groups (only 1 level) as saves vars to root of dict

    var = {}

    def read_vars(f):
        for key in f.variables.keys():
            var[key] = f[key][:]
            try:
                var[key].units = f[key].units
                var[key].long_name = f[key].long_name
            except:
                pass

    with nc.Dataset(file) as f:
        read_vars(f)

        for group in f.groups.keys():
            read_vars(f[group])

    return var

def plot_map(lat,lon,var,var_name,save_location):
    # map of var

    fig = plt.figure(figsize=(12,12))
    projection = crs.PlateCarree()
    ax = plt.axes(projection=projection)

    cs = ax.pcolormesh(lon,lat,var,transform=crs.PlateCarree())
        # vmin=args.var_min, vmax=args.var_max,alpha=args.opacity, cmap=args.colormap, zorder=3)
    plt.title(var_name)
    
    cax,kw = colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    cbar=fig.colorbar(cs,cax=cax,**kw)
    try:
        cbar.set_label(f'{var.long_name} [{var.units}]')
    except:
        pass

    ax.set_aspect('equal')

    savestring = var_name.lower().replace(" ", "_")
    plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')

    return

def plot_map_diff(lat,lon,var1,var2,var_name,plot_name,save_location):
    # map of difference (var1 - var2)

    assert lat.shape == lon.shape == var1.shape == var2.shape
    assert var1.units == var2.units


    diff = var1 - var2
    minmax = np.max(np.abs([np.min(diff),np.max(diff)]))

    fig = plt.figure(figsize=(12,12))
    projection = crs.PlateCarree()
    ax = plt.axes(projection=projection)

    cs = ax.pcolormesh(lon,lat,var1-var2,transform=crs.PlateCarree(),
        vmin=-minmax, vmax=minmax, cmap='bwr', zorder=3)
    plt.title(var_name)
    
    cax,kw = colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    cbar=fig.colorbar(cs,cax=cax,**kw)
    try:
        cbar.set_label(f'Difference {plot_name} [{var1.units}]')
    except:
        pass

    ax.set_aspect('equal')

    savestring = 'diff_'+var_name.lower().replace(" ", "_").replace("-","_")
    plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')


    return

def plot_scatter(var1,var2,var_name,var1_name,var2_name,save_location):
    # var 1: truth
    # var 2: retrieval

    assert var1.shape == var2.shape
    assert var1.units == var2.units

    # mask invalid and flatten
    var1_flat = np.ma.masked_invalid(var1).flatten()
    var2_flat = np.ma.masked_invalid(var2).flatten()

    # combine masks and remove masked data
    var1_mask = var1_flat.mask
    var2_mask = var2_flat.mask
    var12_mask = var1_mask | var2_mask

    var1_flat = var1_flat[~var12_mask].data
    var2_flat = var2_flat[~var12_mask].data

    # make 1D hist
    if (var1_flat == var1_flat[0]).all():
        plt.figure(figsize=(9,9))
        plt.hist(var2_flat, bins=25,label='L2')
        plt.axvline(x=var1_flat[0],linestyle='--',color='k',label='SGM')
        plt.ylabel('Count')
        plt.xlabel(f'{var1_name} [{var1.units}]')
        plt.legend()
        savestring = 'hist_'+var_name.lower().replace(" ", "_")
        plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')

    # make 2D hist
    # fit regressions line and calc stats
    else:
        slope, intercept, r_pearson, p_value, std_err = linregress(var1_flat, var2_flat)
        r2 = r_pearson*r_pearson

        bias = var1_flat-var2_flat
        sigma = np.std(bias)
        bias = np.mean(np.abs(bias))

        # plot

        fig,ax = plt.subplots(figsize=(9,9))
        h = plt.hist2d(var1_flat, var2_flat,bins=100, norm=mpl.colors.LogNorm())
        lims = np.array([np.min([ax.get_xlim(), ax.get_ylim()]),np.max([ax.get_xlim(), ax.get_ylim()])])
        ax.plot(lims, lims, 'k-', alpha=0.3, zorder=3,label='1:1')
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        plt.title('N = {}, R$^2$ = {:.3f}, mean $\sigma$ = {:.3E}, mean bias = {:.3E}'.format(var1_flat.size,r2, sigma, bias))
        plt.xlabel(f'{var1_name} [{var1.units}]')
        plt.ylabel(f'{var2_name} [{var2.units}]')
        plt.plot(lims,lims*slope+intercept,'k--',alpha=0.5,zorder=2,label='y={:.2f}x+{:.2E}'.format(slope, intercept))
        plt.legend()
        cax,kw = colorbar.make_axes(ax,location='right',pad=0.02,shrink=0.5)
        cbar=fig.colorbar(h[-1],cax=cax, extend='neither')
        cbar.set_label('Number of pixels')

        savestring = 'scatter_'+var_name.lower().replace(" ", "_")
        plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')

    return


def pam_nitro(cfg):


    log.info(f"Started PAM")

    # read SGM atm and L2 file

    sgm = read_file(cfg['io']['sgm_atm'])
    l2  = read_file(cfg['io']['l2'])
    basedir  = os.path.dirname(cfg['io']['l2'])

    savedir = os.path.join(basedir, 'figs' )
    Path(savedir).mkdir(exist_ok=True)
    
    plotvars = cfg['plot_list']

    log.info(f"Saving figures to: {savedir}")

    # loop over plotting vars

    for varname in plotvars:

        log.info(f'Plotting {varname}')

        plotvar = plotvars[varname]
        l2_var = l2[plotvar['l2_name']]
        sgm_var = sgm[plotvar['sgm_name']]

        # map of SGM var
        plot_map(sgm['lat'],sgm['lon'],sgm_var, f'SGM {varname}', savedir)

        # map of L2 var
        plot_map(l2['lat'],l2['lon'],l2_var, f'L2 {varname}', savedir)
        
        # diff map (L2 - SGM)
        plot_map_diff(l2['lat'],l2['lon'],l2_var,sgm_var, varname, f'(L2 - SGM) {varname}', savedir)

        # scatter plot SGM vs L2 var
        plot_scatter(sgm_var,l2_var,varname,f'SGM {varname}',f'L2 {varname}', savedir)

    log.info(f'Finished PAM')

    return

if __name__ == '__main__':

    # call with:
    # python pam.py pam_no2.yaml

    # reading yaml config
    cfg = yaml.safe_load(open(sys.argv[1]))

    pam_nitro(cfg)
