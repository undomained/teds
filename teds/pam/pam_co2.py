# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import numpy.typing as npt
from matplotlib.lines import Line2D
from cartopy.feature import LAND, COASTLINE, RIVERS, LAKES
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy

from teds.lib.libNumTools import get_isrf
from teds.gm import vincenty

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


def geo_panel(ax: matplotlib.axes.Axes,
              lon: npt.NDArray[np.float64],
              lat: npt.NDArray[np.float64],
              data: npt.NDArray[np.float64],
              lon_lat_bbox: tuple[float, float, float, float],
              panel_title: str,
              colbar: matplotlib.colorbar,
              cbar_title: str = '',
              valmax: np.float64 | None = None,
              valmin: np.float64 | None = None) -> (
                  tuple[matplotlib.axes.Axes, matplotlib.colorbar]
                  | tuple[matplotlib.axes.Axes,
                          matplotlib.collections.QuadMesh]):

    project = ccrs.PlateCarree()
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


def pam_gm_Tango_Carbon(filen: str,
                        station_name: str,
                        plt_options: str) -> None:
    """ simple plotting routine to visualize gm output
        input:  filen
                station_name ('Jaenschwalde', Belchotow, Lipetsk, Matimba)
    """
    gm_data = Dataset(filen)

    plt.rcParams.update({'font.size': 8})

    global lat_lon_bb, targetp

    lon_lat_bbox = lat_lon_bb[station_name]
    central_point = targetp[station_name]

    if plt_options == 'geoloc':
        fig, ax = plt.subplots(1,
                               1,
                               figsize=(10, 10),
                               dpi=100,
                               subplot_kw={
                                   'projection': ccrs.Orthographic(
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

        plt.plot(gm_data['lon'][:],
                 gm_data['lat'][:],
                 c='red',
                 marker='o',
                 alpha=0.8,
                 linestyle='None',
                 markersize=1,
                 transform=ccrs.PlateCarree())
        plt.plot(central_point[0],
                 central_point[1],
                 c='blue',
                 marker='o',
                 markersize=5,
                 alpha=0.8,
                 transform=ccrs.PlateCarree())

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
        plt.show()

    if plt_options == 'geometry':

        fig, axs = plt.subplots(2,
                                2,
                                figsize=(14, 10),
                                dpi=100,
                                subplot_kw={'projection': ccrs.Orthographic(
                                    central_point[0], central_point[1])},)
        fig.suptitle(station_name, fontsize=16)

        ax = axs[0, 0]
        ax, cbar = geo_panel(ax,
                             gm_data['lon'],
                             gm_data['lat'],
                             gm_data['sza'],
                             lon_lat_bbox,
                             'Solar Zenith Angle',
                             True,
                             'SZA [degree]')

        ax = axs[0, 1]
        ax, cbar = geo_panel(ax,
                             gm_data['lon'],
                             gm_data['lat'],
                             gm_data['vza'],
                             lon_lat_bbox,
                             'Viewing Zenith Angle',
                             True,
                             'VZA [degree]')

        ax = axs[1, 0]
        ax, cbar = geo_panel(ax,
                             gm_data['lon'],
                             gm_data['lat'],
                             gm_data['saa'],
                             lon_lat_bbox,
                             'Solar Azimuth Angle',
                             True,
                             'SAA [degree]'
                             )

        ax = axs[1, 1]
        ax, cbar = geo_panel(ax,
                             gm_data['lon'],
                             gm_data['lat'],
                             gm_data['vaa'],
                             lon_lat_bbox,
                             'Viewing Azimuth Angle',
                             True,
                             'VAA [degree]'
                             )

    plt.show()


def pam_sgm_gps(filen: str, station_name: str, plt_options: str) -> None:
    """ simple plotting routine to visualize sgm-gps output
        input:  filen
                station_name ('Jaenschwalde', Belchotow, Lipetsk, Matimba)
                plt_options ('albedo', 'gas')
    """
    sgmgps_data = Dataset(filen)

    plt.rcParams.update({'font.size': 8})

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
                               subplot_kw={'projection': ccrs.Orthographic(
                                   central_point[0], central_point[1])},)
        fig.suptitle(station_name, fontsize=16)

        ax, cbar = geo_panel(ax,
                             sgmgps_data['lon'],
                             sgmgps_data['lat'],
                             sgmgps_data[plt_options],
                             lon_lat_bbox,
                             plt_options,
                             True,
                             '$A_s$ [1]'
                             )

    if 'X' in plt_options:
        key_list = list(sgmgps_data.variables.keys())
        Xgases = [Xgas for Xgas in key_list if 'X' in Xgas]

        if (plt_options not in Xgases):
            raise Exception('plt_options not well chosen. For gases use an '
                            'imput from the list', Xgases)

        fig, ax = plt.subplots(1, 1, figsize=(14, 10), dpi=100,
                               subplot_kw={'projection': ccrs.Orthographic(
                                   central_point[0], central_point[1])},)
        fig.suptitle(station_name, fontsize=16)

        cbar_title = plt_options + ' [' + sgmgps_data[plt_options].units + ']'
        panel_title = sgmgps_data[plt_options].long_name
        ax, cbar = geo_panel(ax,
                             sgmgps_data['lon'],
                             sgmgps_data['lat'],
                             sgmgps_data[plt_options],
                             lon_lat_bbox,
                             panel_title,
                             True,
                             cbar_title
                             )


def pam_sgm_rad(filen: str,
                station_name: str,
                wavel: float,
                ialt_iact: list[float]) -> None:
    """ simple plotting routine to visualize sgm-gps output
        input:  filen
                station_name ('Jaenschwalde', Belchotow, Lipetsk, Matimba)
                plt_options ('rad_map', 'rad_spec')
                wavel   single wavelength for radiation scene map
    """
    sgmrad_data = Dataset(filen)

    idx = np.array(ialt_iact)
    nspec = np.size(idx, 0)
    if nspec > 4:
        raise Exception(
            'Too many spectra. Reduce the input to maximum 4 spectra')

    cols = ['blue', 'orange', 'green', 'hotpink']

    plt.rcParams.update({'font.size': 8})

    global lat_lon_bb, targetp

    lon_lat_bbox = lat_lon_bb[station_name]
    central_point = targetp[station_name]

    fig, ax1 = plt.subplots(1,
                            1,
                            figsize=(13, 6),
                            dpi=100,
                            subplot_kw={'projection': ccrs.Orthographic(
                                central_point[0], central_point[1])},)
    fig.suptitle(station_name, fontsize=16)

    idxw = np.abs(sgmrad_data['wavelength'][:] - wavel).argmin()
    ax1, mesh = geo_panel(ax1,
                          sgmrad_data['lon'],
                          sgmrad_data['lat'],
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
    fig.suptitle(station_name, fontsize=16)
    for ispec in range(nspec):
        ialt, iact = idx[ispec, :]

        ax1.plot(sgmrad_data['lon'][ialt, iact],
                 sgmrad_data['lat'][ialt, iact],
                 c=cols[ispec],
                 marker='o',
                 markersize=5,
                 alpha=0.8,
                 linestyle='None',
                 transform=ccrs.PlateCarree(),
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


def pam_im(filen: str,
           image_no: int,
           data_max: float,
           title: str = '') -> None:
    """
    input:  filen      filename of input data
            image_no   index pointing to along track readout
    """
    level0 = Dataset(filen)

    # detector dimensions
    plt.rcParams['font.size'] = 14
    n_bins = level0.dimensions['bin'].size
    n_cols = 640
    n_rows = np.int16(n_bins/n_cols)

    # Read data
    science_data = level0['science_data']
    image = science_data['detector_image'][image_no, :]
    image = image.reshape(n_rows, n_cols)

    print(np.max(image))
    print(np.min(image))
    x_values = np.linspace(0, n_cols+1, n_cols+1)
    y_values = np.linspace(0, n_rows+1, n_rows+1)

    # Set up plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlabel('Detector column')
    ax.set_ylabel('Detector row')
    ax.set_title(title)

    psm = ax.pcolormesh(
        x_values,
        y_values,
        image,
        cmap='rainbow',
        vmax=data_max)

    cb = fig.colorbar(psm, ax=ax)
    cb.set_label('binary units [1]')


def pam_l1b(filen_l1b: str,
            filen_sgmrad: str,
            isrf_config: dict,
            plt_options: str,
            ialt: int,
            iact: int,
            spec_nominal: bool) -> None:
    """
    input:  filen      filename of input data
            image_no   index pointing to along track readout
    """

    level1b = Dataset(filen_l1b)
    sgmrad = Dataset(filen_sgmrad)

    if plt_options == 'spectrum':

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
            'wavelength' + level1b['observation_data']['wavelength'].units)
        ax2.set_ylabel('$\\delta$ I [%]')

        ax2.set_xlim([wave_min, wave_max])
        ax1.set_xlim([wave_min, wave_max])

        print(wave_min, wave_max)

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

        print('==>>', err.size)
        num_bins = 201
        bindef = np.arange(num_bins)/20. - 5.
        error_mean = np.mean(err)
        error_std = np.std(err)

        label_txt = 'mean error: '+str("%.2f" % error_mean) + ' std dev: '+str(
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


def pam_l1b_geo(filen1: str,
                filen2: str,
                percentile: float) -> None:

    data1 = Dataset(filen1)
    data2 = Dataset(filen2)

    lat1 = deepcopy(data1['geolocation_data']['latitude'][:])
    lon1 = deepcopy(data1['geolocation_data']['longitude'][:])
    lat2 = deepcopy(data2['geolocation_data']['latitude'][:])
    lon2 = deepcopy(data2['geolocation_data']['longitude'][:])

    data1.close()
    data2.close()

    diff_lat = (lat2-lat1).flatten()
    diff_lon = (lon2-lon1).flatten()

    # the histogram of the data
    fig, axs = plt.subplots(1,
                            3,
                            figsize=(16, 5),
                            dpi=100,)
#    fig.suptitle(station_name, fontsize=16)

    ax = axs[0]
    num_bins = 201
    bindef = (np.arange(num_bins)/200. - 0.5)*6*np.std(diff_lat)
    textstr = f"mean: {np.mean(diff_lat):.2E}\n"
    textstr += f"std. dev.: {np.std(diff_lat):.2E}"
    n, bins, patches = ax.hist(diff_lat,
                               bins=bindef,
                               alpha=0.5,
                               label='latitude',
                               color='blue')
    ax.set_xlabel('$\\Delta$ latitude [degree]')
    ax.set_ylabel('frequency')
    ax.text(0.95,
            0.95,
            textstr,
            transform=ax.transAxes,
            va='top',
            ha='right',
            bbox=dict(boxstyle='round,pad=0.5',
                      fc='white',
                      ec='gray',
                      alpha=0.8))

    ax = axs[1]
    num_bins = 201
    bindef = (np.arange(num_bins)/200. - 0.5)*6*np.std(diff_lon)
    textstr = f"mean: {np.mean(diff_lon):.2E}\n"
    textstr += f"std. dev.: {np.std(diff_lon):.2E}"

    n, bins, patches = ax.hist(diff_lon,
                               bins=bindef,
                               alpha=0.5,
                               label='latitude',
                               color='blue')
    ax.set_xlabel('$\\Delta$ longitude [degree]')
    ax.set_ylabel('frequency')
    ax.text(0.95,
            0.95,
            textstr,
            transform=ax.transAxes,
            va='top',
            ha='right',
            bbox=dict(boxstyle='round,pad=0.5',
                      fc='white',
                      ec='gray',
                      alpha=0.8))

    lat1_rad = np.deg2rad(np.array(lat1).flatten())
    lat2_rad = np.deg2rad(np.array(lat2).flatten())
    lon1_rad = np.deg2rad(np.array(lon1).flatten())
    lon2_rad = np.deg2rad(np.array(lon2).flatten())

    distance = vincenty(lat1_rad, lat2_rad, lon1_rad,  lon2_rad)

    textstr = f"mean: {np.mean(distance):.2f}  m\n"
    textstr += (f" {percentile}%-percentile: "
                f"{np.percentile(distance, percentile):.2f} m")

    ax = axs[2]
    num_bins = 201
    bindef = (np.arange(num_bins)/200.)*6*np.std(distance)

    n, bins, patches = ax.hist(distance,
                               bins=bindef,
                               alpha=0.5,
                               label='latitude',
                               color='blue')
    ax.set_xlabel('$\\Delta r$ [m]')
    ax.set_ylabel('frequency')

    ax.text(0.95,
            0.95,
            textstr,
            transform=ax.transAxes,
            va='top',
            ha='right',
            bbox=dict(boxstyle='round,pad=0.5',
                      fc='white',
                      ec='gray',
                      alpha=0.8))


def pam_l2(filen: str,
           filen_ref: str,
           station_name: str,
           plt_options: str,
           vscale: tuple[float, float] | None = None) -> None:

    level2 = Dataset(filen)
    sgmgps_data = Dataset(filen_ref)

    XCO2 = np.array(level2['XCO2 proxy'][:].data)
    XCO2err = np.array(level2['precision XCO2 proxy'][:].data)
    XCO2sgm = np.array(sgmgps_data['XCO2'][:].data)

    plt.rcParams.update({'font.size': 8})

    if plt_options == 'map':

        global lat_lon_bb, targetp

        lon_lat_bbox = lat_lon_bb[station_name]
        central_point = targetp[station_name]

        fig, axs = plt.subplots(1,
                                2,
                                figsize=(16, 7),
                                dpi=100,
                                subplot_kw={'projection': ccrs.Orthographic(
                                    central_point[0], central_point[1])},)
        fig.suptitle(station_name, fontsize=16)

        if vscale is None:
            XCO2max = np.max(np.array([XCO2.max(), XCO2sgm.max()]))
            XCO2min = np.min(np.array([XCO2.min(), XCO2sgm.min()]))
            print(XCO2max, XCO2min)
        else:
            XCO2min = vscale[0]
            XCO2max = vscale[1]

        ax = axs[0]
        ax, cbar = geo_panel(ax,
                             level2['lon'],
                             level2['lat'],
                             XCO2,
                             lon_lat_bbox,
                             'retrieved XCO2',
                             True,
                             'XCO2 [ppm]',
                             valmax=XCO2max,
                             valmin=XCO2min)

        ax = axs[1]
        ax, cbar = geo_panel(ax,
                             sgmgps_data['lon'],
                             sgmgps_data['lat'],
                             sgmgps_data['XCO2'],
                             lon_lat_bbox,
                             'SGM-GEO XCO2',
                             True,
                             'XCO2 [ppm]',
                             valmax=XCO2max,
                             valmin=XCO2min)

    if plt_options == 'histo':
        # Calculate a linear array of all errors normalized by the
        # spectral standard deviation.
        # the histogram of the data
        fig, axs = plt.subplots(1,
                                2,
                                figsize=(16, 6),
                                dpi=100,
                                )
        ax = axs[0]
        error_norm = (XCO2.flatten() - XCO2sgm.flatten())/XCO2err.flatten()

        num_bins = 201
        bindef = np.arange(num_bins)/20. - 5.
        rmse = np.sqrt(np.mean(np.square(error_norm)))
        textstr = f"mean: {np.mean(error_norm):.2f}\n"
        textstr += f"std. dev.: {np.std(error_norm):.2f}"

        n, bins, patches = ax.hist(error_norm,
                                   bins=bindef,
                                   alpha=0.5,
                                   color='blue')
        ax.set_xlabel('normalized XCO2 error [1]')
        ax.set_ylabel('frequency')
        ax.title.set_text('(L1B-SGM)/$\\sigma$')

        ax.text(0.95,
                0.95,
                textstr,
                transform=ax.transAxes,
                va='top',
                ha='right',
                bbox=dict(boxstyle='round,pad=0.5',
                          fc='white',
                          ec='gray',
                          alpha=0.8))

        ax = axs[1]
        error = (XCO2.flatten() - XCO2sgm.flatten())

        num_bins = 201
        bindef = (np.arange(num_bins)/200. - 0.5)*6*np.std(error)

        rmse = np.sqrt(np.mean(np.square(error)))
        textstr = f"mean: {np.mean(error):.2f} ppm\n"
        textstr += f"std. dev.: {np.std(error):.2f} ppm\n"
        textstr += f"RMSE.: {rmse:.2f} ppm"

        n, bins, patches = ax.hist(error,
                                   bins=bindef,
                                   alpha=0.5,
                                   color='blue')
        ax.set_xlabel('XCO2 error [ppm]')
        ax.set_ylabel('frequency')
        ax.title.set_text('L1B-SGM')

        ax.text(0.95,
                0.95,
                textstr,
                transform=ax.transAxes,
                va='top',
                ha='right',
                bbox=dict(boxstyle='round,pad=0.5',
                          fc='white',
                          ec='gray',
                          alpha=0.8))
