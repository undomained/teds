# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

from netCDF4 import Dataset
from shapely.geometry import shape
from shapely.ops import unary_union
from shapely.prepared import prep
import logging
import numpy as np
import shapely.geometry as sgeom

from .collocation_algorithm import collocate_cloudmask
from .write_module import write_pixel_parameter_to_nc

_logger = logging.getLogger(__name__)


class CloudmaskS2(object):

    def __init__(self, path, aux_path, cartopy_path_directory, clear=False):
        self.path = path
        self.aux_path = aux_path
        self.cartopy_path_directory = cartopy_path_directory
        self.filename = "gauss_conv_300m_berlin_cloudmask.nc"
        self.source = (
            'CO2M prep from Sentinel level2 stitched tiles over Berlin area.')
        self.description = 'Cloud mask from Sentinel2 level1 data over Berlin'

        self.n_pixels = None
        self.lat_orbit = None
        self.lon_orbit = None

        self.clear = clear

        self.cloud_mask_in = None
        self.lats = None
        self.lons = None
        self.cloudmask_orbit = None

        self.cloud_reff_orbit = None
        self.cloud_optical_thick_orbit = None
        self.cloud_top_height_asl_orbit = None
        self.cloud_phase_orbit = None
        self.cloud_fraction_orbit = None
        self.original_cloud_fraction = None
        self.watermask_orbit = None

        self.land = None
        self.lake = None
        self.river = None

        self.varname_cloud = ['cloud_reff',
                              'cloud_optical_thickness',
                              'cloud_top_height_asl',
                              'cloud_phase',
                              'cloud_fraction']

    def read_file(self):

        _logger.info("Reading cloud file over Berlin area {}.".format(
            self.path + "/" + self.filename))
        # Read one in for shape and dimensions.
        root_grp = Dataset(self.path + "/" + self.filename)

        self.lats = root_grp.variables["lat"][:]
        self.lons = root_grp.variables["lon"][:]
        self.cloud_mask_in = root_grp.variables["cloud"][:]

        root_grp.close()

    def collocate(self, julday_orbit, lat, lon, n_pixels, use_watermask=True):
        """Given time and lat-lon box of the groundpixel, return the
        collocated cloud.

        """

        _logger.info("Collocation of cloud file {}".format(self.filename))

        self.lat_orbit = lat
        self.lon_orbit = lon
        self.n_pixels = n_pixels

        self.cloud_reff_orbit = np.zeros((self.n_pixels))
        self.cloud_optical_thick_orbit = np.zeros((self.n_pixels))
        self.cloud_top_height_asl_orbit = np.zeros((self.n_pixels))
        self.cloud_phase_orbit = np.zeros((self.n_pixels))
        self.cloud_fraction_orbit = np.zeros((self.n_pixels))
        self.original_cloud_fraction = np.zeros((self.n_pixels))
        if self.clear:
            self.cloudmask_orbit = np.zeros((self.n_pixels))
        else:
            # Fill regional ensemble with cloudy pixels if there is no data.
            self.cloudmask_orbit = np.ones((self.n_pixels))

        if use_watermask:
            self.prepare_shapes(self.aux_path)
            self.watermask_orbit = np.zeros((self.n_pixels)).astype(int)

            for i in range(0, len(self.lat_orbit)):
                if (
                        (self.is_land(self.lon_orbit[i],
                                      self.lat_orbit[i]) is False) or
                        (self.is_lake(self.lon_orbit[i],
                                      self.lat_orbit[i]) is True) or
                        (self.is_river(self.lon_orbit[i],
                                       self.lat_orbit[i]) is True)):
                    self.watermask_orbit[i] = 1
                else:
                    self.watermask_orbit[i] = 0

        if self.clear:
            return

        if self.cloud_mask_in is None:
            self.read_file()
        (self.cloudmask_orbit,
         idx_cloud,
         idx_clear,
         self.original_cloud_fraction) = collocate_cloudmask(
             self.lat_orbit,
             self.lon_orbit,
             n_pixels,
             self.lats,
             self.lons,
             self.cloud_mask_in,
             fraction_minimum=0.05,
             lonlat_gridded=False)

        # Use contants where there is no info.
        self.cloud_top_height_asl_orbit[idx_cloud] = 2000
        self.cloud_reff_orbit[idx_cloud] = 10
        self.cloud_optical_thick_orbit[idx_cloud] = 10  # To 10 percent.
        self.cloud_phase_orbit[idx_cloud] = 100

        # Use fraction instead of mask.
        self.cloud_fraction_orbit[idx_clear] = 0.1
        self.cloud_fraction_orbit[idx_cloud] = 0.9

        # Cut off fraction below fraction minimum, keep the rest.
        self.original_cloud_fraction[idx_clear] = 0.0

    def post_process(self, albedo_data):
        # Calculate optical thickness using the albedo.
        _logger.info("Adjusting the cloud optical thickness with the albedo")

        albedo_data_band560 = albedo_data[2, :]

        tot_albedo = (
            self.original_cloud_fraction
            * ((10 * (0.14)) / ((10 * (0.14)) + 2))
            + (1 - self.original_cloud_fraction) * albedo_data_band560)

        self.cloud_optical_thick_orbit[:] = (
            (-2 * tot_albedo) / ((tot_albedo - 1) * 0.14))

    def write_output(self, output_dir):

        output_file = output_dir + 'sentinel_cloudmask.nc'

        write_pixel_parameter_to_nc(self.cloudmask_orbit,
                                    'cloudmask',
                                    self.lat_orbit,
                                    self.lon_orbit,
                                    self.npixels,
                                    output_file)

    def write_to_group(self, group, start=0, end=None):

        _logger.info(f"Writing cloudmask to {group.name}")
        nvar_cloud = len(self.varname_cloud)

        data = -999 * np.ones((nvar_cloud, self.n_pixels))
        data[0, :] = self.cloud_reff_orbit
        data[1, :] = self.cloud_optical_thick_orbit
        data[2, :] = self.cloud_top_height_asl_orbit
        data[3, :] = self.cloud_phase_orbit
        data[4, :] = self.cloud_fraction_orbit

        ###############
        for ivar in range(nvar_cloud):
            var = self.varname_cloud[ivar]
            group.variables[var][start:end] = data[ivar, :]

    def write_variables_to_group(self,
                                 group,
                                 var_selection,
                                 start=0,
                                 end=None):

        if 'cloud_top_height_asl' in var_selection:
            group.variables['cloud_top_height_asl'][start:end] = (
                self.cloud_top_height_asl_orbit[:])

        if 'land_fraction' in var_selection:
            land_fraction = np.zeros(self.n_pixels)
            idx_land = np.where(self.watermask_orbit < 1)[0]
            idx_water = np.where(self.watermask_orbit > 0)[0]
            land_fraction[idx_land] = 0.9
            land_fraction[idx_water] = 0.1
            group.variables['land_fraction'][start:end] = land_fraction

        if 'cloud_fraction' in var_selection:

            group.variables['cloud_fraction'][start:end] = (
                self.cloud_fraction_orbit[:])

        if 'original_cloud_fraction' in var_selection:
            group.variables['original_cloud_fraction'][start:end] = (
                self.original_cloud_fraction[:])

    def prepare_shapes(self, path):
        try:
            import shapefile
            shapefiles_path = self.cartopy_path_directory + '/shapefiles'
            land_shp_fname = (
                shapefiles_path + '/natural_earth/physical/ne_50m_land.shp')
            lake_shp_fname = (
                shapefiles_path + '/natural_earth/physical/ne_110m_lakes.shp')
            river_shp_fname = (
                shapefiles_path
                + '/natural_earth/physical/ne_110m_rivers_lake_centerlines.'
                'shp')
            land_sf = shapefile.Reader(land_shp_fname)
            lake_sf = shapefile.Reader(lake_shp_fname)
            river_sf = shapefile.Reader(river_shp_fname)

            land_geom = unary_union(
                [shape(land_sf.shape(rec.oid).__geo_interface__)
                 for rec in land_sf.iterRecords()
                 if rec["featurecla"] != "Null island"])
            lake_geom = unary_union(
                [shape(lake_sf.shape(rec.oid).__geo_interface__)
                 for rec in lake_sf.iterRecords()
                 if "Lake" in rec["featurecla"]])
            river_geom = unary_union(
                [shape(river_sf.shape(rec.oid).__geo_interface__)
                 for rec in river_sf.iterRecords()
                 if "River" in rec["featurecla"]])

        except ImportError:
            import cartopy
            import cartopy.io.shapereader as shpreader
            cartopy.config['pre_existing_data_dir'] = (
                self.cartopy_path_directory)

            land_shp_fname = shpreader.natural_earth(
                resolution='50m', category='physical', name='land')

            land_geom = unary_union(
                [record.geometry
                 for record in shpreader.Reader(land_shp_fname).records()
                 if record.attributes.get('featurecla') != "Null island"])

            lake_shp_fname = shpreader.natural_earth(
                resolution='110m', category='physical', name='lakes')

            lake_geom = unary_union(
                [record.geometry
                 for record in shpreader.Reader(lake_shp_fname).records()
                 if "Lake" in record.attributes.get('featurecla')])

            river_shp_fname = shpreader.natural_earth(
                resolution='110m',
                category='physical',
                name='rivers_lake_centerlines')

            river_geom = unary_union(
                [record.geometry
                 for record in shpreader.Reader(river_shp_fname).records()
                 if "River" in record.attributes.get('featurecla')])

        self.land = prep(land_geom)
        self.lake = prep(lake_geom)
        self.river = prep(river_geom)

    def is_land(self, x, y):
        if not self.land:
            self.prepare_shapes(self.aux_path)
        return self.land.contains(sgeom.Point(x, y))

    def is_lake(self, x, y):
        if not self.lake:
            self.prepare_shapes(self.aux_path)
        return self.lake.contains(sgeom.Point(x, y))

    def is_river(self, x, y):
        if not self.river:
            self.prepare_shapes(self.aux_path)
        return self.river.contains(sgeom.Point(x, y))
