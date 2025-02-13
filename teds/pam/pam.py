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
import subprocess
import numpy.typing as npt
from matplotlib.lines import Line2D
from cartopy.feature import LAND, COASTLINE, RIVERS, LAKES
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable

from teds import log
from teds.lib.libNumTools import get_isrf


lat_lon_bb = {}  # (lon_low, lon_high, lat_low, lat_high)
lat_lon_bb['Matimba'] = (27.4, 27.8, -23.85, -23.5)
lat_lon_bb['Jaenschwalde'] = (14.15, 14.80, 51.65, 52.05)
lat_lon_bb['Belchatow'] = (18.9, 19.7, 51.05, 51.45)
lat_lon_bb['Lipetsk'] = (39.4, 40, 52.75, 52.35)
targetp = {}  # lon, lat
targetp['Matimba'] = (27.610838, -23.6688333)
targetp['Jaenschwalde'] = (14.46277, 51.83472)
targetp['Belchatow'] = (19.330556, 51.266389)
targetp['Lipetsk'] = (39.6945, 52.57123)

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

def plot_map(cfg, lat, lon, var, f_name, var_name, save_location):
    # map of var

    fig = plt.figure(figsize=(12,12))
    # projection = crs.PlateCarree()
    projection = crs.Orthographic(central_longitude=np.mean(lon), central_latitude=np.mean(lat) )
    ax = plt.axes(projection=projection)

    ax.add_feature(LAND)
    ax.add_feature(COASTLINE)
    ax.add_feature(RIVERS)
    ax.add_feature(LAKES)
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    cs = ax.pcolormesh(lon,lat,var,transform=crs.PlateCarree())
        # vmin=args.var_min, vmax=args.var_max,alpha=args.opacity, cmap=args.colormap, zorder=3)
    plt.title(f'{f_name} {var_name}')
    
    cax,kw = colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    cbar=fig.colorbar(cs,cax=cax,**kw)
    try:
        cbar.set_label(f'{var.long_name} [{var.units}]')
        cax.get_xaxis().get_offset_text().set_position((1.15,0))
        if cfg['font_size']:
            cax.get_xaxis().get_offset_text().set_fontsize(cfg['font_size'])
        cbar.ax.xaxis.OFFSETTEXTPAD = -15

    except:
        pass

    ax.set_aspect('equal')

    if save_location:
        savestring = f_name.lower().replace(" ", "_") +'_'+ var_name.lower().replace(" ", "_")
        plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')
    else:
        plt.show()

    return

def plot_map_diff(cfg, lat, lon, var1, var2, f1_name, f2_name, var_name, save_location):
    # map of difference (var2 - var1)

    assert lat.shape == lon.shape == var1.shape == var2.shape
    assert var1.units == var2.units

    # plot abs diff
    diff = var2 - var1
    minmax = np.max(np.abs([np.min(diff),np.max(diff)]))

    fig = plt.figure(figsize=(12,12))
    # projection = crs.PlateCarree()
    projection = crs.Orthographic(central_longitude=np.mean(lon), central_latitude=np.mean(lat) )
    ax = plt.axes(projection=projection)

    cs = ax.pcolormesh(lon,lat,diff,transform=crs.PlateCarree(),
        vmin=-minmax, vmax=minmax, cmap='bwr', zorder=3)
    plt.title(f'({f2_name} - {f1_name}) {var_name}')
    
    cax,kw = colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    cbar=fig.colorbar(cs,cax=cax,**kw)
    try:
        cbar.set_label(f'Difference {var_name} [{var1.units}]', labelpad = 25)
        cax.get_xaxis().get_offset_text().set_position((1.15,0))
        if cfg['font_size']:
            cax.get_xaxis().get_offset_text().set_fontsize(cfg['font_size'])
        cbar.ax.xaxis.OFFSETTEXTPAD = -15
    except:
        pass

    ax.set_aspect('equal')

    savestring = 'diff_'+var_name.lower().replace(" ", "_").replace("-","_")
    if save_location:
        plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')
    else:
        plt.show()

    # plot rel diff
    reldiff = (var2 - var1)/var1 * 100
    minmax = np.max(np.abs([np.min(reldiff),np.max(reldiff)]))

    fig = plt.figure(figsize=(12,12))
    # projection = crs.PlateCarree()
    projection = crs.Orthographic(central_longitude=np.mean(lon), central_latitude=np.mean(lat) )
    ax = plt.axes(projection=projection)
    cs = ax.pcolormesh(lon,lat,reldiff,transform=crs.PlateCarree(),
        vmin=-minmax, vmax=minmax, cmap='bwr', zorder=3)
    plt.title(f'({f2_name} - {f1_name})/{f1_name}*100 %')
    
    cax,kw = colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.5)
    cbar=fig.colorbar(cs,cax=cax,**kw)
    try:
        cbar.set_label(f'Difference {var_name} [%]', labelpad = 25)
        cax.get_xaxis().get_offset_text().set_position((1.15,0))
        if cfg['font_size']:
            cax.get_xaxis().get_offset_text().set_fontsize(cfg['font_size'])
        cbar.ax.xaxis.OFFSETTEXTPAD = -15
    except:
        pass

    ax.set_aspect('equal')

    if save_location:
        savestring = 'reldiff_'+var_name.lower().replace(" ", "_").replace("-","_")
        plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')
    else:
        plt.show()

    return

def plot_hist(cfg, var1, var2, f1_name, f2_name, var_name, save_location, req=None, req_name=None):
    # var1: x
    # var2: y
    # f1_name: var1 name for plot
    # f2_name: var2 name for plot
    # save_location: path to figs folder
    # req: requirement value
    # req_name: name of requirement for plot

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

    # make 1D hist of constants SGM vars (clouds for example)
    if (var1_flat != var1_flat[0]).all():
        plt.figure(figsize=(9,9))
        plt.hist(var2_flat, bins=25,label=f2_name)
        plt.axvline(x=var1_flat[0],linestyle='--',color='k',label=f1_name)
        plt.ylabel('Count')
        plt.xlabel(f'{var_name} [{var1.units}]')
        plt.legend()

        if save_location:
            savestring = 'hist_'+var_name.lower().replace(" ", "_")
            plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()

    # make 1D error hist
    else:
        diff = var2_flat- var1_flat

        stdev_diff = np.std(diff)
        mean_diff = np.mean(np.abs(diff))

        plt.figure(figsize=(9,9))
        plt.hist(diff, bins=100)
        plt.ylabel('Count')
        plt.xlabel(f'difference {var_name} [{var1.units}]')
        plt.title('stdev = {:.3E}, mean = {:.3E}'.format(stdev_diff, mean_diff))
        if req_name: # add requirement lines
            plt.axvline(x=req, color='gray',alpha=0.5, zorder = 1, label=req_name)
            plt.axvline(x=-req, color='gray',alpha=0.5, zorder = 1)
            plt.legend()
        savestring = 'hist_'+var_name.lower().replace(" ", "_")
        if save_location:
            plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()
    return


def plot_scatter(cfg, var1, var2, f1_name, f2_name, var_name, save_location, req=None, req_name=None):
    # var1: x
    # var2: y
    # f1_name: var1 name for plot
    # f2_name: var2 name for plot
    # save_location: path to figs folder
    # req: requirement value
    # req_name: name of requirement for plot

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

    # all values are equal in one or both arrays; scatterplot is useless.
    if (var1_flat == var1_flat[0]).all() or (var2_flat == var2_flat[0]).all():
        return
    
    # fit regressions line and calc stats
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
    plt.title('N = {}, R$^2$ = {:.3f}, stdev(error) = {:.3E}, mean(error) = {:.3E}'.format(var1_flat.size,r2, sigma, bias))
    plt.xlabel(f'{f1_name} {var_name} [{var1.units}]', labelpad=25)
    plt.ylabel(f'{f2_name} {var_name} [{var2.units}]')
    plt.plot(lims,lims*slope+intercept,'k--',alpha=0.5,zorder=2,label='y={:.2f}x+{:.2E}'.format(slope, intercept))
    
    if req_name == 'bias': # add requirement lines
        plt.plot(lims,lims+lims*req,'k-.',alpha=0.5,zorder=2,label=f'{req_name} requirement') 
        plt.plot(lims,lims-lims*req,'k-.',alpha=0.5,zorder=2)
    elif req_name == 'precision': 
        plt.plot(lims,lims*slope+intercept+req,'k-.',alpha=0.5,zorder=2,label=f'{req_name} requirement') 
        plt.plot(lims,lims*slope+intercept-req,'k-.',alpha=0.5,zorder=2)
    
    plt.legend()
    cax,kw = colorbar.make_axes(ax,location='right',pad=0.02,shrink=0.5)
    cbar=fig.colorbar(h[-1],cax=cax, extend='neither')
    cbar.set_label('Number of pixels')
    if cfg['font_size']:
        cax.get_xaxis().get_offset_text().set_fontsize(cfg['font_size'])
    
    if save_location:
        savestring = 'scatter_'+var_name.lower().replace(" ", "_")
        plt.savefig(f'{save_location}/{savestring}.png',format='png', dpi=1000, bbox_inches='tight')
    else:
        plt.show()
    return

def pam_l2(cfg, savedir):

    # loop over plotting vars
    plotvars = cfg['l2']['plot_list']

    for varname in plotvars:
        plotvar = plotvars[varname]
        log.debug(f'Plotting {varname}')

        # read files
        if plotvar['f1'] in cfg['io']:
            file1 = cfg['io'][plotvar['f1']]
        else:
            file1 = plotvar['f1']
        f1 = read_file(file1, [plotvar['f1_var'], 'latitude', 'longitude'])

        f1_var = f1[plotvar['f1_var']]
        f1_lat = f1['latitude']
        f1_lon = f1['longitude']
        f1_name = plotvar['f1_name']

        # map of f1 var
        plot_map(cfg, f1_lat, f1_lon, f1_var, f1_name, varname, savedir)

        # if two files are present, plot f2 map, diff and scatter
        if 'f2' in plotvar:

            if plotvar['f2'] in cfg['io']:
                file2 = cfg['io'][plotvar['f2']]
            else:
                file2 = plotvar['f2']
            f2 = read_file(file2, [plotvar['f2_var'], 'latitude', 'longitude'])

            f2_lat = f2['latitude']
            f2_lon = f2['longitude']
            f2_var = f2[plotvar['f2_var']]
            f2_name = plotvar['f2_name']

            # map of f2 var
            plot_map(cfg, f2_lat, f2_lon, f2_var, f2_name, varname, savedir)
            
            # diff map (f2 - f1)
            plot_map_diff(cfg, f1_lat, f1_lon, f1_var, f2_var, f1_name, f2_name, varname, savedir)
            
            # hist plot f1 vs f2 var
            if 'req' in plotvar: # plot requirement line 
                plot_hist(cfg, f1_var, f2_var, f1_name, f2_name, varname, savedir,
                            req=float(plotvar['req']),
                            req_name=plotvar['req_name'])
                plot_scatter(cfg, f1_var, f2_var, f1_name, f2_name, varname, savedir,
                            req=float(plotvar['req']),
                            req_name=plotvar['req_name'])
            else:
                plot_hist(cfg, f1_var, f2_var, f1_name, f2_name, varname, savedir)
                plot_scatter(cfg, f1_var, f2_var, f1_name, f2_name, varname, savedir)
    return


def geo_panel(ax: mpl.axes.Axes,
              lon: npt.NDArray[np.float64],
              lat: npt.NDArray[np.float64],
              data: npt.NDArray[np.float64],
              lon_lat_bbox: tuple[float, float, float, float],
              panel_title: str,
              colbar: mpl.colorbar,
              cbar_title: str = '',
              valmax: np.float64 | None = None,
              valmin: np.float64 | None = None) -> (
                  tuple[mpl.axes.Axes, mpl.colorbar]
                  | tuple[mpl.axes.Axes,
                          mpl.collections.QuadMesh]):

    project = crs.PlateCarree()
    ax.set_extent(lon_lat_bbox)
    ax.set_title(panel_title)
    ax.add_feature(LAND)
    ax.add_feature(COASTLINE)
    ax.add_feature(RIVERS)
    ax.add_feature(LAKES)
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    if valmin is None:
        vvmin = np.nanmin(data)
    else:
        vvmin = valmin
    if valmax is None:
        vvmax = np.nanmax(data)
    else:
        vvmax = valmax

    mesh = ax.pcolormesh(lon, lat, data, alpha=0.9, transform=project,
                         cmap='cividis', vmax=vvmax, vmin=vvmin)
    if colbar:
        cbar = plt.colorbar(mesh,
                            ax=ax,
                            orientation='vertical',
                            fraction=0.04,
                            pad=0.10)
        cbar.set_label(cbar_title)
        return (ax, cbar)
    else:
        return (ax, mesh)


def pam_gm(filen: str,
            station_name: str,
            plt_options: str,
            save_dir: str) -> None:
    """ simple plotting routine to visualize gm output
        input:  filen
                station_name ('Jaenschwalde', Belchotow, Lipetsk, Matimba)
    """
    gm_data = nc.Dataset(filen)

    # plt.rcParams.update({'font.size': 8})

    global lat_lon_bb, targetp

    lon_lat_bbox = lat_lon_bb[station_name]
    central_point = targetp[station_name]

    if plt_options == 'geoloc':
        fig, ax = plt.subplots(1,
                               1,
                               figsize=(10, 10),
                               dpi=100,
                               subplot_kw={
                                   'projection': crs.Orthographic(
                                       central_point[0],
                                       central_point[1])},)

        ax.set_extent(lon_lat_bbox)
        ax.add_feature(LAND)
        ax.add_feature(COASTLINE)
        ax.add_feature(RIVERS)
        ax.add_feature(LAKES)
        ax.set_title(station_name + " target")

        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gls = ax.gridlines(draw_labels=True)
        gls.top_labels = False   # suppress top labels
        gls.right_labels = False  # suppress right labels

        plt.plot(gm_data['longitude'][:],
                 gm_data['latitude'][:],
                 c='red',
                 marker='o',
                 alpha=0.8,
                 linestyle='None',
                 markersize=1,
                 transform=crs.PlateCarree())
        plt.plot(central_point[0],
                 central_point[1],
                 c='blue',
                 marker='o',
                 markersize=5,
                 alpha=0.8,
                 transform=crs.PlateCarree())

        point_obs = Line2D([0],
                           [0],
                           label='observation point',
                           marker='o',
                           markersize=1,
                           c='red',
                           linestyle='None',
                           alpha=0.8)
        point_tar = Line2D([0],
                           [0],
                           label='target point',
                           marker='o',
                           markersize=5,
                           c='blue',
                           linestyle='None',
                           alpha=0.8)
        plt.legend(handles=[point_obs, point_tar])
        if save_dir:
            plt.savefig(f'{save_dir}/geolocation.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()

    if plt_options == 'geometry':

        fig, axs = plt.subplots(2,
                                2,
                                figsize=(14, 10),
                                dpi=100,
                                subplot_kw={'projection': crs.Orthographic(
                                    central_point[0], central_point[1])},)
        fig.suptitle(station_name)#, fontsize=16)

        ax = axs[0, 0]
        ax, cbar = geo_panel(ax,
                             gm_data['longitude'],
                             gm_data['latitude'],
                             gm_data['solarzenithangle'],
                             lon_lat_bbox,
                             'Solar Zenith Angle',
                             True,
                             'SZA [degree]')

        ax = axs[0, 1]
        ax, cbar = geo_panel(ax,
                             gm_data['longitude'],
                             gm_data['latitude'],
                             gm_data['viewingzenithangle'],
                             lon_lat_bbox,
                             'Viewing Zenith Angle',
                             True,
                             'VZA [degree]')

        ax = axs[1, 0]
        ax, cbar = geo_panel(ax,
                             gm_data['longitude'],
                             gm_data['latitude'],
                             gm_data['solarazimuthangle'],
                             lon_lat_bbox,
                             'Solar Azimuth Angle',
                             True,
                             'SAA [degree]'
                             )

        ax = axs[1, 1]
        ax, cbar = geo_panel(ax,
                             gm_data['longitude'],
                             gm_data['latitude'],
                             gm_data['viewingazimuthangle'],
                             lon_lat_bbox,
                             'Viewing Azimuth Angle',
                             True,
                             'VAA [degree]'
                             )
        if save_dir:
            plt.savefig(f'{save_dir}/geometry.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()


def pam_sgm_gps(filen: str, station_name: str, plt_options: str, save_dir: str) -> None:
    """ simple plotting routine to visualize sgm-gps output
        input:  filen
                station_name ('Jaenschwalde', Belchotow, Lipetsk, Matimba)
                plt_options ('albedo', 'gas')
    """
    sgmgps_data = nc.Dataset(filen)

    # plt.rcParams.update({'font.size': 8})

    global lat_lon_bb, targetp

    lon_lat_bbox = lat_lon_bb[station_name]
    central_point = targetp[station_name]

    if 'albedo' in plt_options:

        key_list = list(sgmgps_data.variables.keys())
        S2_bands = [band for band in key_list if 'B' in band]
        if ((plt_options not in S2_bands)):
            raise Exception('plt_options not well chosen. For albedo, use an '
                            f'imput from the list{S2_bands}')

        fig, ax = plt.subplots(1,
                               1,
                               figsize=(14, 10),
                               dpi=100,
                               subplot_kw={'projection': crs.Orthographic(
                                   central_point[0], central_point[1])},)
        fig.suptitle(station_name)#, fontsize=16)

        ax, cbar = geo_panel(ax,
                             sgmgps_data['longitude'],
                             sgmgps_data['latitude'],
                             sgmgps_data[plt_options],
                             lon_lat_bbox,
                             plt_options,
                             True,
                             '$A_s$ [1]'
                             )
        if save_dir:
            plt.savefig(f'{save_dir}/sgm_raw_{sgmgps_data[plt_options].name}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()

    if 'col' in plt_options:
        key_list = list(sgmgps_data.variables.keys())
        Xgases = [Xgas for Xgas in key_list if 'col' in Xgas]

        if (plt_options not in Xgases):
            raise Exception('plt_options not well chosen. For gases use an '
                            'imput from the list', Xgases)

        fig, ax = plt.subplots(1, 1, figsize=(14, 10), dpi=100,
                               subplot_kw={'projection': crs.Orthographic(
                                   central_point[0], central_point[1])},)
        fig.suptitle(station_name)#, fontsize=16)

        cbar_title = plt_options + ' [' + sgmgps_data[plt_options].units + ']'
        panel_title = sgmgps_data[plt_options].long_name
        ax, cbar = geo_panel(ax,
                             sgmgps_data['longitude'],
                             sgmgps_data['latitude'],
                             sgmgps_data[plt_options],
                             lon_lat_bbox,
                             panel_title,
                             True,
                             cbar_title
                             )
        if save_dir:
            plt.savefig(f'{save_dir}/sgm_raw_{sgmgps_data[plt_options].name}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()

def pam_sgm_rad(file_rad: str,
                file_gm: str,
                station_name: str,
                wavel: float,
                ialt_iact: list[float],
                save_dir: str) -> None:
    """ simple plotting routine to visualize sgm-gps output
        input:  file_rad
                file_gm
                station_name ('Jaenschwalde', Belchotow, Lipetsk, Matimba)
                plt_options ('rad_map', 'rad_spec')
                wavel   single wavelength for radiation scene map
    """
    sgmrad_data = nc.Dataset(file_rad)
    gm_data = nc.Dataset(file_gm)

    idx = np.array(ialt_iact)
    nspec = np.size(idx, 0)
    if nspec > 4:
        raise Exception(
            'Too many spectra. Reduce the input to maximum 4 spectra')

    cols = ['blue', 'orange', 'green', 'hotpink']

    # plt.rcParams.update({'font.size': 8})

    global lat_lon_bb, targetp

    lon_lat_bbox = lat_lon_bb[station_name]
    central_point = targetp[station_name]

    fig, ax1 = plt.subplots(1,
                            1,
                            figsize=(13, 6),
                            dpi=100,
                            subplot_kw={'projection': crs.Orthographic(
                                central_point[0], central_point[1])},)
    fig.suptitle(station_name)#, fontsize=16)

    idxw = np.abs(sgmrad_data['wavelength'][:] - wavel).argmin()
    ax1, mesh = geo_panel(ax1,
                          gm_data['longitude'],
                          gm_data['latitude'],
                          sgmrad_data['radiance'][:, :, idxw],
                          lon_lat_bbox,
                          f"{sgmrad_data['wavelength'][idxw]:.2f}" + ' nm',
                          colbar=False,
                          )

    divider = make_axes_locatable(ax1)
    ax2 = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)

    fig.add_axes(ax2)
    plt.colorbar(mesh, cax=ax2)

    ax3 = divider.new_horizontal(size="100%", pad=1, axes_class=plt.Axes)
    fig.add_axes(ax3)
    fig.suptitle(station_name)#, fontsize=16)
    for ispec in range(nspec):
        ialt, iact = idx[ispec, :]

        ax1.plot(gm_data['longitude'][ialt, iact],
                 gm_data['latitude'][ialt, iact],
                 c=cols[ispec],
                 marker='o',
                 markersize=5,
                 alpha=0.8,
                 linestyle='None',
                 transform=crs.PlateCarree(),
                 label=f'(ialt, iact) =({ialt}.{iact})')

        ax3.plot(sgmrad_data['wavelength'][:],
                 sgmrad_data['radiance'][ialt, iact, :],
                 color=cols[ispec],
                 label=f'(ialt, iact) =({ialt}.{iact})',
                 alpha=0.6)
    ax3.set_xlabel('wavelength' + sgmrad_data['wavelength'].units)
    ax3.set_ylabel('I ' + sgmrad_data['radiance'].units)
    ax3.legend()
    ax1.legend()

    if save_dir:
        plt.savefig(f'{save_dir}/sgm_rad.png',format='png', dpi=1000, bbox_inches='tight')
    else:
        plt.show()


def pam_im(filen: str,
           image_no: int,
           data_max: float,
           save_dir: str,
           title: str = '') -> None:
    """
    input:  filen      filename of input data
            image_no   index pointing to along track readout
    """
    level0 = nc.Dataset(filen)

    # detector dimensions
    # plt.rcParams['font.size'] = 14

    # Read data
    science_data = level0['science_data']
    image = science_data['detector_image'][image_no,:,:]
    n_rows, n_cols = image.shape

    # print(np.max(image))
    # print(np.min(image))
    x_values = np.linspace(0, n_cols+1, n_cols+1)
    y_values = np.linspace(0, n_rows+1, n_rows+1)

    # Set up plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlabel('Detector column')
    ax.set_ylabel('Detector row')
    ax.set_title(title)

    # print(data_max)

    psm = ax.pcolormesh(
        x_values,
        y_values,
        image,
        cmap='rainbow')#,
        # vmax=data_max)

    cb = fig.colorbar(psm, ax=ax)
    cb.set_label('binary units [1]')

    if save_dir:
        plt.savefig(f'{save_dir}/im_{image_no}.png',format='png', dpi=1000, bbox_inches='tight')
    else:
        plt.show()

def pam_l1b(filen_l1b: str,
            filen_sgmrad: str,
            isrf_config: dict,
            plt_options: str,
            ialt: int,
            iact: int,
            spec_nominal: bool,
            save_dir: str) -> None:
    """
    input:  filen      filename of input data
            image_no   index pointing to along track readout
    """

    level1b = nc.Dataset(filen_l1b)
    sgmrad = nc.Dataset(filen_sgmrad)
    
    if plt_options == 'residuals':

        # define isrf function
        wave_lbl = sgmrad['wavelength'][:].data
        wave = level1b['observation_data']['wavelength'][iact, :].data
        isrf_convolution = get_isrf(wave, wave_lbl, isrf_config)
        sgmrad_conv = isrf_convolution(sgmrad['radiance'][ialt, iact, :].data)

        wave_min = wave_lbl.min()
        wave_max = wave_lbl.max()
        fig, axs = plt.subplots(2, 1, figsize=(10, 10), dpi=100,)
        fig.suptitle(
            f'Spectral residuals: Ialt = {ialt}, Iact = {iact}', fontsize=16)
        ax1 = axs[0]

        if spec_nominal:
            idx = np.where((wave >= 1590) & (wave <= 1675))
        else:
            idx = np.where((wave >= wave_min) & (wave <= wave_max))

        l1b = np.array(level1b['observation_data']['radiance'][ialt, iact, :])
        ax1.plot(wave[idx],
                 l1b[idx],
                 color='blue',
                 label='level 1B',
                 alpha=0.6)

        ax1.plot(wave,
                 sgmrad_conv,
                 color='darkblue',
                 linestyle=':',
                 label='convolved level 1B',
                 alpha=0.6)

        ax1.plot(sgmrad['wavelength'][:].data,
                 sgmrad['radiance'][ialt, iact, :].data,
                 color='blue',
                 label='sgm line-by-line',
                 alpha=0.2)

        ax1.set_ylabel('I ' + level1b['observation_data']['radiance'].units)
        ax1.set_xlim([wave_min, wave_max])
        ax1.legend()
        ax1.set_xlim([1590, 1675])

        ax2 = axs[1]
        error = np.array((
            level1b['observation_data']['radiance'][ialt, iact, :].data
            - sgmrad_conv) / sgmrad_conv * 100.)
        ax2.plot(wave[idx],
                 error[idx],
                 color='blue',
                 alpha=0.6,
                 label='level-1B - conv. sgm radiance')
        ax2.set_xlim([wave_min, wave_max])

        ax2.plot(
            level1b['observation_data']['wavelength'][iact, :].data,
            level1b['observation_data']['radiance_stdev'][ialt, iact, :].data
            / sgmrad_conv * 100.,
            color='orange',
            linestyle='-',
            label='$\\pm$ one-sigma',
            alpha=0.2)
        ax2.plot(
            level1b['observation_data']['wavelength'][iact, :].data,
            -level1b['observation_data']['radiance_stdev'][ialt, iact, :].data
            / sgmrad_conv * 100.,
            color='orange',
            linestyle='-',
            alpha=0.2)
        ax2.legend()
        ax2.set_xlabel(
            'wavelength ' + level1b['observation_data']['wavelength'].units)
        ax2.set_ylabel('$\\delta$ I [%]')

        ax2.set_xlim([wave_min, wave_max])
        ax1.set_xlim([wave_min, wave_max])

        if save_dir:
            plt.savefig(f'{save_dir}/l1b_residual_{ialt}_{iact}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()

        # print(wave_min, wave_max)

    if plt_options == 'histo':
        # Calculate a linear array of all errors normalized by the
        # spectral standard deviation.
        sgmrad_conv = np.empty(level1b['observation_data']['radiance'].shape)
        nalt = np.size(sgmrad_conv, axis=0)
        nact = np.size(sgmrad_conv, axis=1)
        for iact in range(nact):
            # define isrf function
            wave_lbl = sgmrad['wavelength'][:].data
            wave = level1b['observation_data']['wavelength'][iact, :].data
            isrf_convolution = get_isrf(wave, wave_lbl, isrf_config)
            for ialt in range(nalt):
                sgmrad_conv[ialt, iact, :] = isrf_convolution(
                    sgmrad['radiance'][ialt, iact, :].data)

        error = (
            (level1b['observation_data']['radiance'][:].data - sgmrad_conv)
            / level1b['observation_data']['radiance_stdev'][:].data)
        if spec_nominal:
            err = np.array([])
            for iact in range(nact):
                wave = level1b['observation_data']['wavelength'][iact, :].data
                idx = np.where((wave >= 1590) & (wave <= 1675))
                err = np.append(err, error[:, iact, idx].flatten())
        else:
            err = error.flatten()

        # print('==>>', err.size)
        num_bins = 201
        bindef = np.arange(num_bins)/5. - 20.
        error_mean = np.mean(err)
        error_std = np.std(err)

        label_txt = 'mean error: '+str("%.2E" % error_mean) + ' std dev: '+str(
            "%.4f" % error_std)

        # the histogram of the data
        plt.figure(figsize=(10, 6))
        n, bins, patches = plt.hist(err,
                                    bins=bindef,
                                    alpha=0.5,
                                    label=label_txt,
                                    color='blue')
        plt.xlabel('normalized radiance error [1]')
        plt.ylabel('frequency')
        plt.legend()
        if save_dir:
            plt.savefig(f'{save_dir}/l1b_hist.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()

    if plt_options == 'spectrum':
        rad = level1b['observation_data']['radiance'][ialt,iact,:]*1e-4
        wvl = level1b['observation_data']['wavelength'][iact,:]
        sgm_rad= sgmrad['radiance'][ialt, iact, :]*1e-4
        sgm_wvl = sgmrad['wavelength'][:]

        plt.figure(figsize=(15, 5))
        plt.plot(sgm_wvl,sgm_rad, label='SGM', alpha=0.5)
        plt.plot(wvl,rad, label='L1B')
        plt.xlabel('wavelength [nm]')
        plt.ylabel('Radiance [ph sr-1 cm-2 nm-2 s-1]')
        plt.legend()
        plt.title(f'Ialt = {ialt}, Iact = {iact}', fontsize=16)

        if save_dir:
            plt.savefig(f'{save_dir}/l1b_spectrum_{ialt}_{iact}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()

    if plt_options == 'snr' or plt_options == 'snr_req':
        noise = level1b['observation_data']['radiance_stdev'][ialt,iact,:]
        signal = level1b['observation_data']['radiance'][ialt,iact,:]
        wvl = level1b['observation_data']['wavelength'][iact,:]
        snr = signal/noise

        # if plt_options == 'snr_req':
        #     snr[np.where(signal<6.24E16)] = np.ma.masked

        plt.figure(figsize=(15, 5))
        plt.plot(wvl,snr)
        plt.xlabel('wavelength [nm]')
        plt.ylabel('SNR [-]')
        plt.title(f'SNR: Ialt = {ialt}, Iact = {iact}', fontsize=16)
        if plt_options == 'snr_req':
            plt.axhline(y=240, color='r', linestyle='--', label = 'SNR2 = 240' )
            plt.legend()

        if save_dir:
            plt.savefig(f'{save_dir}/l1b_snr_{ialt}_{iact}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()


        noise = level1b['observation_data']['radiance_stdev'][ialt,:,:]
        signal = level1b['observation_data']['radiance'][ialt,:,:]
        snr = signal/noise

        # scatter plot
        fig,ax = plt.subplots(figsize=(9,9))
        h = plt.hist2d(signal.flatten(), snr.flatten() ,bins=100, norm=mpl.colors.LogNorm())
        if plt_options == 'snr_req':
            snr_mean = np.ma.mean(snr)
            plt.title(f'Ialt = {ialt}. All ACT pixels. Mean SNR: {snr_mean:.0f}')
            plt.axhline(y=240, color='r', linestyle='--', label = 'SNR2 = 240' )
            plt.axvline(x=6.24E16, color='k', linestyle='--', label = 'L2 = 6.24E16' )
            plt.legend(loc='upper left')
        plt.xlabel('Radiance [ph sr-1 cm-2 nm-2 s-1]')
        plt.ylabel('SNR [-]')
        # plt.title('N = {}, R$^2$ = {:.3f}, stdev(error) = {:.3E}, mean(error) = {:.3E}'.format(var1_flat.size,r2, sigma, bias))
        # plt.plot(lims,lims*slope+intercept,'k--',alpha=0.5,zorder=2,label='y={:.2f}x+{:.2E}'.format(slope, intercept))
        
        cax,kw = colorbar.make_axes(ax,location='right',pad=0.02,shrink=0.5)
        cbar=fig.colorbar(h[-1],cax=cax, extend='neither')
        cbar.set_label('Number of pixels')

        if save_dir:
            plt.savefig(f'{save_dir}/l1b_snr_Scatter_{ialt}_{iact}.png',format='png', dpi=1000, bbox_inches='tight')
        else:
            plt.show()


def pam_nitro(cfg):

    log.info(f'Starting PAM')

    # set fontsize
    if 'font_size' in cfg:
        mpl.rcParams.update({'font.size': cfg['font_size']})

    if cfg['save']:
        # create folder for figs
        savedir = os.path.join(cfg['io']['base_dir'], 'figs' )
        Path(savedir).mkdir(exist_ok=True)
        log.info(f"Saving figures to: {savedir}")
    else:
        savedir = None
    
    if cfg['gm']['run']:
        log.info(f'Plotting GM')
        pam_gm(cfg['io']['gm'], cfg['scene'], 'geoloc', savedir)
        pam_gm(cfg['io']['gm'], cfg['scene'], 'geometry', savedir)

    if cfg['sgm']['run']:
        log.info(f'Plotting SGM geo')
        pam_sgm_gps(cfg['io']['sgm_atm_raw'], cfg['scene'], 'col_no2', savedir)
        pam_sgm_gps(cfg['io']['sgm_atm_raw'], cfg['scene'], 'albedo_B02', savedir)

        pam_sgm_rad(cfg['io']['sgm_rad'], cfg['io']['gm'], cfg['scene'], cfg['sgm']['wvl'], cfg['sgm']['ialt_iact'], savedir)

    if cfg['im']['run']:
        log.info(f'Plotting L1A')
        pam_im(cfg['io']['l1a'], cfg['im']['ialt'], 0, savedir, title=f'ALT image {cfg['im']['ialt']}')
    
    if cfg['l1b']['run']:
        log.info(f'Plotting L1B')
        isrf_config={}
        isrf_config['type'] = 'Gaussian'  #type of ISRF, currently only Gaussian or generalized_normal
        isrf_config['fwhm'] = cfg['l1b']['isrf_fwhm']  #fwhm [nm]

        for plot_option in cfg['l1b']['plot_options']:
            pam_l1b(cfg['io']['l1b'], cfg['io']['sgm_rad'], isrf_config, plot_option, cfg['l1b']['ialt'], cfg['l1b']['iact'], False, savedir)

    if cfg['l2']['run']:
        log.info(f'Plotting L2')
        pam_l2(cfg, savedir) 


    if cfg['create_pdf']:
        log.info(f'Creating PDF')
        subprocess.run(f'magick $(ls -tr {savedir}/*.png) {savedir}/pam_figures.pdf', shell=True) 
        subprocess.run(f'magick -density 300x300 -quality 80 -compress jpeg \
                        {savedir}/pam_figures.pdf {savedir}/pam_figures_compressed.pdf',
                        shell=True) 

    log.info(f'Finished PAM')

    return

if __name__ == '__main__':

    # call with:
    # python pam.py <yaml config file>

    # reading yaml config
    cfg = yaml.safe_load(open(sys.argv[1]))

    pam_nitro(cfg)
