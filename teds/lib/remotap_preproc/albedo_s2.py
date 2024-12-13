from netCDF4 import Dataset
import numpy as np
from .collocation_algorithm import interpolate_to_orbit
import logging

_logger = logging.getLogger(__name__)


class AlbedoS2(object):

    def __init__(self, path):
        self.path = path
        self.filename_template = "gauss_conv_300m_berlin_albedo_band{}.nc"
        self.name_bands = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B8A", "B09", "B11", "B12"]
        self.wavelength_bands = [443, 490, 560, 665, 705, 740, 783, 865, 945, 1610, 2190]
        self.source = 'CO2M prep from Sentinel level2 stitched tiles over Berlin area.'
        self.description = 'Average surface albedo from Sentinel2 level2 data over Berlin'

        self.n_pixels = None
        self.lat_orbit = None
        self.lon_orbit = None

        self.albedo_in = None
        self.lats = None
        self.lons = None
        self.albedo = None

        # Band 1 - Coastal aerosol	0.443	60
        # Band 2 - Blue	0.490	10
        # Band 3 - Green	0.560	10
        # Band 4 - Red	0.665	10
        # Band 5 - Vegetation Red Edge	0.705	20
        # Band 6 - Vegetation Red Edge	0.740	20
        # Band 7 - Vegetation Red Edge	0.783	20
        # Band 8 - NIR	0.842	10
        # Band 8A - Vegetation Red Edge	0.865	20
        # Band 9 - Water vapour	0.945	60
        # Band 10 - SWIR - Cirrus	1.375	60
        # Band 11 - SWIR	1.610	20
        # Band 12 - SWIR	2.190
        
    def get_nbands(self):
        return len(self.wavelength_bands)

    def read_file(self):

        _logger.info("Reading albedo file")

        # Read one in for shape and dimensions.
        root_grp = Dataset(self.path + "/" + self.filename_template.format(self.name_bands[0]))
        nlat = root_grp.dimensions['nlat'].size
        nlon = root_grp.dimensions['nlon'].size

        self.lats = root_grp.variables["lat"][:]
        self.lons = root_grp.variables["lon"][:]

        nlon = self.lons.size
        nlat = self.lats.size

        root_grp.close()

        # Read all.
        self.albedo_in = np.zeros((len(self.wavelength_bands), nlat, nlon))
        for i, band in enumerate(self.name_bands):
            _logger.debug("Reading in band {}".format(band))
            root_grp = Dataset(self.path + "/" + self.filename_template.format(band))
            self.albedo_in[i, :, :] = root_grp.variables['albedo'][:]
            root_grp.close()

    def collocate(self, julday_orbit, lat, lon, n_pixels):
        """Given time and lat-lon box of the groundpixel, return the collocated albedo"""

        if self.albedo_in is None:
            self.read_file()

        _logger.info("Collocation of albedo file.")
        self.lat_orbit = lat
        self.lon_orbit = lon
        self.n_pixels = n_pixels

        self.albedo = interpolate_to_orbit(self.lat_orbit, self.lon_orbit, n_pixels, self.lats,
                                           self.lons, self.albedo_in, len(self.wavelength_bands), method='nearest')

    def post_process(self, fill_data, bands):

        # The albedo file contains data on a equally spaced lon/lat grid.
        _logger.info("Filling missing albedo data after collocation with xbrdf from gome2 and modis used in the orbit ensemble.")
        indices = np.where((self.lat_orbit < np.min(self.lats)) | (self.lat_orbit > np.max(self.lats)) |
                           (self.lon_orbit < np.min(self.lons)) | (self.lon_orbit > np.max(self.lons)))
        for i, band in enumerate(self.wavelength_bands):

            # for each band, find the closest related band frequency and replace nan with the fill_data.
            diff_array = np.abs(bands - self.wavelength_bands[i])
            related_band = diff_array.argmin()
            self.albedo[i, indices] = fill_data[related_band, indices]

    def write_output(self, output_dir):

        output_file = output_dir + 'albedo_sentinel2.nc'
        da = Dataset(output_file, 'w', format='NETCDF4')

        da.createDimension('Nbands', len(self.wavelength_bands))
        da.createDimension('Npixels_orbit', self.n_pixels)

        da.createVariable('Bands', 'f4', ('Nbands'))
        da.variables['Bands'][:] = self.wavelength_bands

        da.createVariable('albedo', 'f4', ('Nbands', 'Npixels_orbit'))

        da.variables['albedo'][:] = self.albedo

        da.createVariable('lon', 'f4', ('Npixels_orbit'))
        da.createVariable('lat', 'f4', ('Npixels_orbit'))
        da.variables['lat'][:] = self.lat_orbit
        da.variables['lon'][:] = self.lon_orbit

    def write_to_group(self, group, start=0, end=None):

        _logger.info("Writing Xbrdf information to {}".format(group.name))

        group.variables['Bands_combine'][:] = self.wavelength_bands
        group.variables['Xbrdf_wave_combine'][:, start:end] = self.albedo
