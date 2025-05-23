from scipy.stats import linregress
import cartopy.crs as crs
import logging
import matplotlib as mpl
import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import numpy.typing as npt
import sys
import yaml


def read_file(file: str) -> dict[str, nc.Variable]:
    # Read in netcdf file, return as dict.
    # Descends into groups (only 1 level) as saves vars to root of dict
    var = {}

    def read_vars(f: nc.Dataset) -> None:
        for key in f.variables.keys():
            var[key] = f[key][:]
            try:
                var[key].units = f[key].units
                var[key].long_name = f[key].long_name
            except Exception:
                pass

    with nc.Dataset(file) as f:
        read_vars(f)

        for group in f.groups.keys():
            read_vars(f[group])

    return var


def plot_map(lat: npt.NDArray[np.float64],
             lon: npt.NDArray[np.float64],
             var: nc.Variable,
             var_name: str,
             save_location: str) -> None:
    # map of var

    fig = plt.figure(figsize=(12, 12))
    projection = crs.PlateCarree()
    ax = plt.axes(projection=projection)

    cs = ax.pcolormesh(lon, lat, var, transform=crs.PlateCarree())
    plt.title(var_name)

    cax, kw = colorbar.make_axes(ax, location='bottom', pad=0.05, shrink=0.5)
    cbar = fig.colorbar(cs, cax=cax, **kw)
    try:
        cbar.set_label(f'{var.long_name} [{var.units}]')
    except Exception:
        pass

    ax.set_aspect('equal')

    savestring = var_name.lower().replace(" ", "_")
    plt.savefig(f'{save_location}/{savestring}.png',
                format='png',
                dpi=1000,
                bbox_inches='tight')


def plot_map_diff(lat: npt.NDArray[np.float64],
                  lon: npt.NDArray[np.float64],
                  var1: nc.Variable,
                  var2: nc.Variable,
                  var_name: str,
                  plot_name: str,
                  save_location: str) -> None:
    # Map of difference (var1 - var2)

    assert lat.shape == lon.shape == var1.shape == var2.shape
    assert var1.units == var2.units

    diff = var1[:] - var2[:]
    minmax = np.max(np.abs([np.min(diff), np.max(diff)]))

    fig = plt.figure(figsize=(12, 12))
    projection = crs.PlateCarree()
    ax = plt.axes(projection=projection)

    cs = ax.pcolormesh(lon,
                       lat,
                       var1[:]-var2[:],
                       transform=crs.PlateCarree(),
                       vmin=-minmax,
                       vmax=minmax,
                       cmap='bwr',
                       zorder=3)
    plt.title(var_name)

    cax, kw = colorbar.make_axes(ax, location='bottom', pad=0.05, shrink=0.5)
    cbar = fig.colorbar(cs, cax=cax, **kw)
    try:
        cbar.set_label(f'Difference {plot_name} [{var1.units}]')
    except Exception:
        pass

    ax.set_aspect('equal')

    savestring = 'diff_'+var_name.lower().replace(" ", "_").replace("-", "_")
    plt.savefig(f'{save_location}/{savestring}.png',
                format='png',
                dpi=1000,
                bbox_inches='tight')


def plot_scatter(var1: nc.Variable,
                 var2: nc.Variable,
                 var_name: str,
                 var1_name: str,
                 var2_name: str,
                 save_location: str) -> None:
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

    # fit regressions line and calc stats

    slope, intercept, r_pearson, p_value, std_err = linregress(
        var1_flat, var2_flat)

    r2 = r_pearson*r_pearson
    bias = var1_flat-var2_flat
    sigma = np.std(bias)
    bias = np.mean(np.abs(bias))

    # plot

    fig, ax = plt.subplots(figsize=(9, 9))
    h = plt.hist2d(var1_flat, var2_flat, bins=100, norm=mpl.colors.LogNorm())
    lims = np.array([np.min([ax.get_xlim(), ax.get_ylim()]),
                     np.max([ax.get_xlim(), ax.get_ylim()])])
    ax.plot(lims, lims, 'k-', alpha=0.3, zorder=3, label='1:1')
    ax.set_aspect('equal')
    ax.set_xlim(lims[0], lims[1])
    ax.set_ylim(lims[0], lims[1])
    plt.title('N = {}, R$^2$ = {:.3f}, mean $\\sigma$ = {:.3E}, '
              'mean bias = {:.3E}'.format(var1_flat.size, r2, sigma, bias))
    plt.xlabel(f'{var1_name} [{var1.units}]')
    plt.ylabel(f'{var2_name} [{var2.units}]')
    plt.plot(lims,
             lims*slope+intercept,
             'k--',
             alpha=0.5,
             zorder=2,
             label='y={:.2f}x+{:.2E}'.format(slope, intercept))
    plt.legend()
    cax, kw = colorbar.make_axes(ax, location='right', pad=0.02, shrink=0.5)
    cbar = fig.colorbar(h[-1], cax=cax, extend='neither')
    cbar.set_label('Number of pixels')

    savestring = 'scatter_'+var_name.lower().replace(" ", "_")
    plt.savefig(f'{save_location}/{savestring}.png',
                format='png',
                dpi=1000,
                bbox_inches='tight')


def pam_nitro(logger: logging.Logger, cfg: dict) -> None:
    logger.info("Started PAM")

    # Read SGM atm and L2 file
    sgm = read_file(cfg['sgm_atm_file'])
    l2 = read_file(cfg['l2_file'])

    plotvars = cfg['pam']['plot_list']

    logger.info(f"Saving figures to: {cfg['pam']['figure_dir']}")

    # loop over plotting vars

    for varname in plotvars:

        logger.info(f'Plotting {varname}')

        plotvar = plotvars[varname]
        l2_var = l2[plotvar['l2_name']]
        sgm_var = sgm[plotvar['sgm_name']]

        # map of SGM var
        plot_map(sgm['lat'][:],
                 sgm['lon'][:],
                 sgm_var,
                 f'SGM {varname}',
                 cfg['pam']['figure_dir'])

        # Map of L2 var
        plot_map(l2['lat'][:],
                 l2['lon'][:],
                 l2_var,
                 f'L2 {varname}',
                 cfg['pam']['figure_dir'])

        # Diff map (L2 - SGM)
        plot_map_diff(l2['lat'][:],
                      l2['lon'][:],
                      l2_var,
                      sgm_var,
                      varname,
                      f'(L2 - SGM) {varname}',
                      cfg['pam']['figure_dir'])

        # scatter plot SGM vs L2 var
        plot_scatter(sgm_var[:],
                     l2_var[:],
                     varname,
                     f'SGM {varname}',
                     f'L2 {varname}',
                     cfg['pam']['figure_dir'])

    logger.info('Finished PAM')


if __name__ == '__main__':

    # call with:
    # python pam.py pam_no2.yaml

    # or with logging to file:
    # python pam.py pam_no2.yaml pam_no2.log

    # Reading yaml config
    cfg = yaml.safe_load(open(sys.argv[1]))

    # cfg = cfg['pam']

    loglevel = logging.INFO

    # setup the logging to screen and to file
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')

    if len(sys.argv) > 2:
        fh = logging.FileHandler(sys.argv[2], mode='w')
        fh.setLevel(logging.ERROR)
        fh.setFormatter(formatter)

        ch = logging.StreamHandler()
        ch.setLevel(loglevel)
        ch.setFormatter(formatter)

        logging.basicConfig(level=loglevel, handlers=[ch, fh])

        logging.info(f'Logging to file: {sys.argv[2]}')

    else:
        ch = logging.StreamHandler()
        ch.setLevel(loglevel)
        ch.setFormatter(formatter)
        logging.basicConfig(level=loglevel, handlers=[ch])

    logger = logging.getLogger()

    pam_nitro(logger, cfg)
