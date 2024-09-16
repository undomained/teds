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

def read_file(file, var_names):
    # read in netcdf file, return as dict
    # descends into groups (only 1 level) as saves vars to root of dict

    var = {}

    def read_vars(f):
        for key in f.variables.keys():
            if key in var_names:
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

def plot_map(lat, lon, var, f_name, var_name, save_location):
    # map of var

    fig = plt.figure(figsize=(12,12))
    # projection = crs.PlateCarree()
    projection = crs.Orthographic(central_longitude=np.mean(lon), central_latitude=np.mean(lat) )
    ax = plt.axes(projection=projection)

    cs = ax.pcolormesh(lon,lat,var,transform=crs.PlateCarree())
        # vmin=args.var_min, vmax=args.var_max,alpha=args.opacity, cmap=args.colormap, zorder=3)
    plt.title(f'{f_name} {var_name}')
    
    cax,kw = colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    cbar=fig.colorbar(cs,cax=cax,**kw)
    try:
        cbar.set_label(f'{var.long_name} [{var.units}]', labelpad = 25)
    except:
        pass

    ax.set_aspect('equal')

    savestring = f_name.lower().replace(" ", "_") +'_'+ var_name.lower().replace(" ", "_")
    plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')

    return

def plot_map_diff(lat, lon, var1, var2, f1_name, f2_name, var_name, save_location):
    # map of difference (var2 - var1)

    assert lat.shape == lon.shape == var1.shape == var2.shape
    assert var1.units == var2.units


    diff = var2 - var1
    minmax = np.max(np.abs([np.min(diff),np.max(diff)]))

    fig = plt.figure(figsize=(12,12))
    # projection = crs.PlateCarree()
    projection = crs.Orthographic(central_longitude=np.mean(lon), central_latitude=np.mean(lat) )
    ax = plt.axes(projection=projection)

    cs = ax.pcolormesh(lon,lat,var2-var1,transform=crs.PlateCarree(),
        vmin=-minmax, vmax=minmax, cmap='bwr', zorder=3)
    plt.title(f'({f2_name} - {f1_name}) {var_name}')
    
    cax,kw = colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    cbar=fig.colorbar(cs,cax=cax,**kw)
    try:
        cbar.set_label(f'Difference {var_name} [{var1.units}]', labelpad = 25)
    except:
        pass

    ax.set_aspect('equal')

    savestring = 'diff_'+var_name.lower().replace(" ", "_").replace("-","_")
    plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')


    return

def plot_scatter(var1, var2, f1_name, f2_name, var_name, save_location):
    # var 1: x
    # var 2: y

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
        plt.hist(var2_flat, bins=25,label=f2_name)
        plt.axvline(x=var1_flat[0],linestyle='--',color='k',label=f1_name)
        plt.ylabel('Count')
        plt.xlabel(f'{var_name} [{var1.units}]')
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
        plt.xlabel(f'{f1_name} {var_name} [{var1.units}]', labelpad=25)
        plt.ylabel(f'{f2_name} {var_name} [{var2.units}]')
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

    # set fontsize
    mpl.rcParams.update({'font.size': cfg['font_size']})

    # create folder for figs
    savedir = os.path.join(cfg['io']['base_dir'], 'figs' )
    Path(savedir).mkdir(exist_ok=True)
    log.info(f"Saving figures to: {savedir}")

    # loop over plotting vars
    plotvars = cfg['plot_list']

    for varname in plotvars:
        plotvar = plotvars[varname]
        log.info(f'Plotting {varname}')

        # read files
        f1 = read_file(cfg['io'][plotvar['f1']], [plotvar['f1_var'], 'lat', 'lon'])

        f1_var = f1[plotvar['f1_var']]
        f1_lat = f1['lat']
        f1_lon = f1['lon']
        f1_name = plotvar['f1_name']

        # map of f1 var
        plot_map(f1_lat, f1_lon, f1_var, f1_name, varname, savedir)

        # if two files are present, plot f2 map, diff and scatter
        if 'f2' in plotvar:
            f2 = read_file(cfg['io'][plotvar['f2']], [plotvar['f2_var'], 'lat', 'lon'])

            f2_lat = f2['lat']
            f2_lon = f2['lon']
            f2_var = f2[plotvar['f2_var']]
            f2_name = plotvar['f2_name']

            # map of f2 var
            plot_map(f2_lat, f2_lon, f2_var, f2_name, varname, savedir)
            
            # diff map (f2 - f1)
            plot_map_diff(f1_lat, f1_lon, f1_var, f2_var, f1_name, f2_name, varname, savedir)

            # scatter plot f1 vs f2 var
            plot_scatter(f1_var, f2_var, f1_name, f2_name, varname, savedir)

    log.info(f'Finished PAM')

    return

if __name__ == '__main__':

    # call with:
    # python pam.py pam_no2.yaml

    # reading yaml config
    cfg = yaml.safe_load(open(sys.argv[1]))

    pam_nitro(cfg)
