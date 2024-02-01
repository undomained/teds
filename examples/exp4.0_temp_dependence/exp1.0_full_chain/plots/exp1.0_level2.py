#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 21:06:07 2023

@author: jochen
"""
import netCDF4 as nc
import numpy as np
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

path = '/home/jochen/TANGO_E2ES/EndtoEndProject/data/'

#true state from sgm file
file0 = path+'/interface_data/sgm/Tango_Carbon_sgm_atmosphere_exp1.0.nc'
#corrected
file1 = path+'/interface_data/level2/Tango_Carbon_l2_exp1.0.nc'

sgm_data = nc.Dataset(file0)
albedo = deepcopy(sgm_data['albedo'][:])
xco2_true = deepcopy(sgm_data['XCO2'][:])
#xch4_true = deepcopy(sgm_data['col_ch4'][:])

l2_data = nc.Dataset(file1)
lat = deepcopy(l2_data['lat'][:])
lon = deepcopy(l2_data['lon'][:])
xco2_proxy      = deepcopy(l2_data['XCO2 proxy'][:])
xch4_proxy      = deepcopy(l2_data['XCH4 proxy'][:])
xco2_ns         = deepcopy(l2_data['XCO2'][:])
xch4_ns         = deepcopy(l2_data['XCH4'][:])
prec_xco2_proxy = deepcopy(l2_data['precision XCO2 proxy'][:])
prec_xch4_proxy = deepcopy(l2_data['precision XCH4 proxy'][:])
prec_xco2_ns    = deepcopy(l2_data['precision XCO2'][:])
prec_xch4_ns    = deepcopy(l2_data['precision XCH4'][:])
l2_data.close()

# l2_data = nc.Dataset(file2)
# xco2_proxy_afgl = deepcopy(l2_data['XCO2 proxy'][:])
# l2_data.close()

xco2_plume = xco2_true
xco2_plume = np.where(xco2_plume<405.5, np.nan, xco2_plume)
cfig = '2D'

if(cfig=='2D'):
    # plotting the lat long grid
    projPC = ccrs.PlateCarree()
    # Read shape file and extract geometries [for only Germany]
    reader = shpreader.Reader(path+'cartopy/germany/natural.shp')
    # To plot individually
    # ['forest', 'park', 'riverbank', 'water'] four kinds of the data
    water_geometry     = [riv.geometry for riv in reader.records() if riv.attributes["type"] == "water"]
    forest_geometry    = [riv.geometry for riv in reader.records() if riv.attributes["type"] == "forest"]
    park_geometry      = [riv.geometry for riv in reader.records() if riv.attributes["type"] == "park"]
    riverbank_geometry = [riv.geometry for riv in reader.records() if riv.attributes["type"] == "riverbank"]

    central_lon, central_lat = 14.4, 58.1

    lat_low = 51.65
    lat_high = 52.00
    lon_low = 14.20
    lon_high = 14.70
    extent = [lon_low, lon_high, lat_low, lat_high]

    # Create four polar axes and access them through the returned array
    fig, axs = plt.subplots(1, 3, figsize=(16, 5), dpi=100,   subplot_kw={
        'projection': ccrs.Orthographic(central_lon, central_lat)},)
    plt.subplots_adjust(wspace=0.4)
    ax = axs[0]
    ax.set_title('scene')
    ax.set_extent(extent)
    ax.add_feature(cfeature.LAND)
    ax.add_geometries(forest_geometry, crs=projPC, facecolor='forestgreen', alpha=0.5, zorder=0)
    ax.add_geometries(water_geometry, crs=projPC, facecolor='blue', alpha=0.5, zorder=0)

    gls = ax.gridlines(draw_labels=True)
    gls.top_labels = False   # suppress top labels
    gls.right_labels = False  # suppress right labels
#    ax.plot(central_lon, central_lat, marker='o', transform=ccrs.PlateCarree(), color='blue', markersize=5)
    mesh1 = ax.pcolormesh(lon, lat, albedo, alpha=1.0, transform=ccrs.PlateCarree(),
                         cmap='cividis', vmax=0.4, vmin=0)
    mesh2 = ax.pcolormesh(lon, lat, xco2_plume, vmin=390., vmax=435., cmap = 'BuPu', transform=ccrs.PlateCarree())

    cbar1 = plt.colorbar(mesh1, ax=ax, orientation='vertical', fraction=0.04, pad=0.05)
    cbar1.set_label('albedo [1]')

#   ===========================================================================
    ax = axs[1]
    ax.set_title('level 2 proxy sgm prior')
    ax.set_extent(extent)
    ax.add_feature(cfeature.LAND)
    ax.add_geometries(forest_geometry, crs=projPC, facecolor='forestgreen', alpha=0.5, zorder=0)

    gls = ax.gridlines(draw_labels=True)
    gls.top_labels = False   # suppress top labels
    gls.right_labels = False  # suppress right labels

    ax.plot(central_lon, central_lat, marker='o', transform=ccrs.PlateCarree(), color='blue', markersize=5)
    mesh = ax.pcolormesh(lon, lat, xco2_proxy, alpha=1.0, transform=ccrs.PlateCarree(),
                         cmap='BuPu', vmax=435, vmin=390)

    cbar = plt.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.04, pad=0.05)
    cbar.set_label('XCO$_2$ [ppm]')

    plt.show()
#
    plt.savefig('plots/xco2_error.png',)

if(cfig == 'histo'):

    xco2_proxy_err = (np.array(xco2_proxy).flatten()-np.array(xco2_true).flatten())/(np.array(xco2_true).flatten()-405.)*100.
    xco2_proxy_err = xco2_proxy_err[xco2_true.flatten()>408.]

    xco2_proxy_afgl_err = (np.array(xco2_proxy_afgl).flatten()-np.array(xco2_true).flatten())/(np.array(xco2_true).flatten()-405.)*100.
    xco2_proxy_afgl_err = xco2_proxy_afgl_err[xco2_true.flatten()>408]

    num_bins = 120
    bindef = np.arange(num_bins)/1. - 30.
    prox_lab = 'XCO$_2^\mathrm{prox}$' + ' ( m ='+str("%.2f" % np.mean(xco2_proxy_err)) + \
        ' ppm $ \sigma = $' + str("%.2f" % np.std(xco2_proxy_err)) + ' ppm)'
    prox_afgl = 'XCO$_2^\mathrm{prox}$' + ' ( m ='+str("%.2f" % np.mean(xco2_proxy_afgl_err)) + \
        ' ppm $ \sigma = $' + str("%.2f" % np.std(xco2_proxy_afgl_err)) + ' ppm)'
#    prox_corr = 'XCO2$_2^\mathrm{corr}$' + ' ( m ='+str("%.2f" % np.mean(xco2_proxy_corr_err)) + \
#        ' ppm $ \sigma = $' + str("%.2f" % np.std(xco2_proxy_corr_err)) + ' ppm)'
    # the histogram of the data

    plt.figure(figsize=(10,6))
    n, bins, patches = plt.hist(xco2_proxy_err, bins=bindef, alpha=0.5, label = prox_lab,color = 'blue')
    n, bins, patches = plt.hist(xco2_proxy_afgl_err, bins=bindef, alpha=0.5, label = prox_afgl, color = 'orange')
#    n, bins, patches = plt.hist(xco2_proxy_corr_err, bins=bindef, alpha=0.5, label = prox_corr, color = 'green')

    plt.xlabel('XCO$_2$ enhancement error [%]')
    plt.ylabel('frequency')
    plt.xlim([-30,80])
    plt.legend()
    plt.show()
    plt.savefig('xco2_histo.png',)
