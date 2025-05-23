# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions dealing with input/output.

These are functions for formatting the output to stdout, for
reading/writing input/output files, and tests for things like whether
a file is readable/writable.

"""
from datetime import datetime
from datetime import timezone
from netCDF4 import Dataset
from pathlib import Path
from pyquaternion import Quaternion
import numpy as np
import numpy.typing as npt
import platform
import xarray as xr
import yaml

from .interpolate_navigation_data import interpolate_navigation_data
from .types import BinningTable
from .types import CKD
from .types import CKDDark
from .types import CKDNoise
from .types import CKDNonlin
from .types import CKDPRNU
from .types import CKDRadiometric
from .types import CKDSpectral
from .types import CKDStray
from .types import CKDSwath
from .types import L1
from .types import ProcLevel
from teds import log
from teds.gm.io import nc_write_geometry
from teds.gm.types import Geometry
from teds.gm.types import Navigation
from teds.lib.io import get_git_commit_hash


def read_binning_table(binning_file: str,
                       binning_table_id: int,
                       n_detector_rows: int,
                       n_detector_cols: int) -> BinningTable:
    """Read the binning pattern.

    Parameters
    ----------
    binning_file
        Path of file with binning patterns.
    binning_table_id
        Identifier of binning pattern, 0 is unbinned.
    n_detector_rows
        Number of detector rows.
    n_detector_cols
        Number of detector columns.

    Returns
    -------
        Binning table consisting of a two-dimensional (unbinned)
        detector map of bin indices and one-dimensional list of number
        of detector pixels per bin.

    """
    if binning_table_id == 0:
        bin_indices = np.arange(n_detector_rows * n_detector_cols,
                                dtype=np.int32).reshape(
            n_detector_rows, n_detector_cols)
        count_table: npt.NDArray[np.int32] = np.ones(
            n_detector_rows * n_detector_cols, np.int32)
    else:
        if not Path(binning_file).is_file():
            _text = "binning file" if binning_file is None else binning_file
            raise SystemExit(f"ERROR: {_text} not found")
        binning_table = xr.load_dataset(binning_file,
                                        group=f'Table_{binning_table_id}')
        bin_indices = binning_table.binning_table.values
        uniq_bin_indices, tmp_count_table = np.unique(bin_indices,
                                                      return_counts=True)
        count_table = np.zeros(np.max(uniq_bin_indices) + 1, np.int32)
        count_table[uniq_bin_indices] = tmp_count_table
    return BinningTable(bin_indices, count_table)


def read_proc_level(filename: str) -> ProcLevel:
    """Read the processing level of a L1 file.

    Parameters
    ----------
    filename
        Path of an L1A file, L1B file, or anything in between.

    Returns
    -------
        Processing (calibration) level of the given file.

    """
    ds = Dataset(filename)
    if ds.product_type in ('L1A', 'L1B', 'SGM'):
        return ProcLevel[ds.product_type.lower()]
    return ProcLevel(ds.l1x_level)


def read_ckd(filename: str) -> CKD:
    """Read a CKD file.

    All available data are read, so unused data are not skipped and
    missing data are not checked.

    Parameters
    ----------
    filename
        Path of CKD file.

    Returns
    -------
        Calibration key data (CKD).

    """
    nc = Dataset(filename)
    # Dark CKD
    ckd_dark = CKDDark(
        offset=nc['dark/offset'][:].data.ravel(),  # counts
        current=nc['dark/current'][:].data.ravel()  # counts/s
    )
    # Noise. conversion_gain can be a detector map or 1 value,
    # optionally dependent on temperature. Best practice is to use the
    # same value for all pixels and move a dependence to the PRNU CKD.
    ckd_noise = CKDNoise(
        conversion_gain=1 / nc['noise/g'][:].data.ravel(),  # e/counts
        read_noise=nc['noise/n'][:].data.ravel()  # counts
    )
    # Nonlinearity
    ckd_nonlin = CKDNonlin(
        observed=nc['nonlinearity/observed'][:].data,  # counts
        expected=nc['nonlinearity/expected'][:].data  # counts
    )
    if (
            not (np.diff(ckd_nonlin.expected) > 0).all()
            or not (np.diff(ckd_nonlin.observed) > 0).all()):
        raise SystemExit(f"ERROR: nonlinearity data in {filename} are not "
                         "strictly increasing")
    # PRNU including quantum efficiency
    ckd_prnu = CKDPRNU(nc['prnu/prnu'][:].data.ravel())
    # Stray light
    n_kernels = nc['stray'].dimensions['kernel'].size
    kernel_n_rows = nc['stray/kernel_rows'][:].data
    kernel_n_cols = nc['stray/kernel_cols'][:].data
    kernel_fft_sizes = nc['stray/kernel_fft_sizes'][:].data
    kernels_fft = nc['stray/kernels_fft'][:].data
    ckd_kernels_fft: list[npt.NDArray[np.complex128]] = []
    # Rearrange memory layout of kernel FFTs if packed
    for i_kernel in range(n_kernels):
        global_beg = 2 * sum(kernel_fft_sizes[0:i_kernel])
        global_end = global_beg + 2 * kernel_fft_sizes[i_kernel]
        kernel_fft = kernels_fft[global_beg:global_end]
        # If not packed
        if kernel_fft_sizes[i_kernel] == (kernel_n_rows[i_kernel]
                                          * kernel_n_cols[i_kernel]):
            n_comp = 2
            kernel_fft_shaped = kernel_fft.reshape(
                n_comp, kernel_n_rows[i_kernel], kernel_n_cols[i_kernel])
            re = kernel_fft_shaped[0, :, :]
            im = kernel_fft_shaped[1, :, :]
            fft_full = re + 1j * im
        else:
            fft_re_packed = np.reshape(
                kernel_fft[::2],
                (kernel_n_rows[i_kernel] // 2 + 1, kernel_n_cols[i_kernel]))
            fft_im_packed = np.reshape(kernel_fft[1::2], fft_re_packed.shape)
            fft_packed = fft_re_packed + 1j * fft_im_packed
            fft_full = np.zeros((kernel_n_rows[i_kernel],
                                 kernel_n_cols[i_kernel]),
                                dtype='complex128')
            mid = kernel_n_rows[i_kernel] // 2 + 1
            fft_full[:mid, :] = fft_packed
            fft_full[mid::, 0] = np.conjugate(fft_packed[mid - 2:0:-1, 0])
            fft_full[mid::, 1:] = np.conjugate(
                fft_packed[mid - 2:0:-1, kernel_n_cols[i_kernel]:0:-1])
        ckd_kernels_fft.append(fft_full)
    ckd_stray = CKDStray(ckd_kernels_fft,
                         nc['stray/eta'][:].data,
                         nc['stray/weights'][:].data,
                         nc['stray/edges'][:].data)
    # Swath
    ckd_swath = CKDSwath(nc['swath/act_angle'][:].data,
                         nc['swath/act_map'][:].data,
                         nc['swath/wavelength_map'][:].data,
                         nc['swath/row_map'][:].data,
                         nc['swath/col_map'][:].data,
                         nc['swath/line_of_sight'][:].data)
    # Smile (spectral distortion) description:
    # wavelength at each detector pixel
    ckd_spectral = CKDSpectral(nc['spectral/wavelength'][:].data)  # nm
    # Radiometric, nm-1 sr-1 m-2 e count-1
    ckd_radiometric = CKDRadiometric(nc['radiometric/radiometric'][:].data)
    return CKD(
        nc.dimensions['detector_row'].size,
        nc.dimensions['detector_column'].size,
        nc['pixel_mask'][:].astype(bool).ravel(),  # False good, True bad
        ckd_dark,
        ckd_noise,
        ckd_nonlin,
        ckd_prnu,
        ckd_stray,
        ckd_swath,
        ckd_spectral,
        ckd_radiometric)


def read_geometry(l1_product: L1, config: dict) -> Geometry:
    """Read geometry from an original geometry file or an L1 file.

    Geometry is not used for calculations in L1 processing. The
    dictionaries are only used to read some settings.

    Parameters
    ----------
    l1_product
        L1 product (signal and detector settings).
    config
        Configuration parameters.

    Returns
    -------
        Coordinates and angles.

    """
    filename = config['io']['geometry']
    alt_beg = config['alt_beg']
    alt_end = None if config['alt_end'] is None else config['alt_end']+1
    with Dataset(filename) as root:
        groups = list(root.groups)
    if 'geolocation_data' in groups:  # L1 file
        # Assumption: the input file currently processed or a file
        # with the same geometry data
        nc = Dataset(filename)
        grp = nc['geolocation_data']
        geo = Geometry(grp['latitude'][alt_beg:alt_end, :],
                       grp['longitude'][alt_beg:alt_end, :],
                       grp['height'][alt_beg:alt_end, :],
                       grp['solar_zenith'][alt_beg:alt_end, :],
                       grp['solar_azimuth'][alt_beg:alt_end, :],
                       grp['sensor_zenith'][alt_beg:alt_end, :],
                       grp['sensor_azimuth'][alt_beg:alt_end, :])
    else:
        nc = Dataset(filename)
        if 'height' in groups:
            _height = nc['height'][:]
        else:
            _height = np.zeros_like(nc['lat'][alt_beg:alt_end, :])
        geo = Geometry(nc['lat'][alt_beg:alt_end, :],
                       nc['lon'][alt_beg:alt_end, :],
                       _height,
                       nc['sza'][alt_beg:alt_end, :],
                       nc['saa'][alt_beg:alt_end, :],
                       nc['vza'][alt_beg:alt_end, :],
                       nc['vaa'][alt_beg:alt_end, :])
    return geo


def copy_geometry(l1a_filename: str,
                  geo_filename: str,
                  i_alt_beg: int,
                  l1_product: L1) -> None:
    """Copy geolocation data from the geometry file directly to the
    L1B product.

    This is a placeholder function until geolocation is properly
    implemented.

    Parameters
    ----------
    l1a_filename
        Path to the input L1A file. Only used for reading the
        configuration used in the IM which contains the original
        alt_beg.
    geo_filename
        Path to a file containing the geometry to be copied.
    i_alt_beg
        First along track position of geometry to be copied. This is
        in reference to the current list of L1 products. In reference
        to the geometry file, this value is added to alt_beg read
        from the L1A file.
    l1_product
        L1 product (signal and detector settings). Used for
        getting the current size of the L1 product (number of along
        track positions) and for storing the geometry.

    """
    nc_l1a = Dataset(l1a_filename)
    configuration = yaml.safe_load(nc_l1a['configuration'][:])
    alt_beg = configuration['alt_beg']
    _beg = i_alt_beg + alt_beg
    _end = _beg + max(l1_product.signal.shape[0], l1_product.spectra.shape[0])
    nc_geo = Dataset(geo_filename)
    if 'height' in nc_geo.variables:
        _height = nc_geo['height'][:]
    else:
        _height = np.zeros_like(nc_geo['lat'][:])
    l1_product.geometry = Geometry(nc_geo['latitude'][_beg:_end, :],
                                   nc_geo['longitude'][_beg:_end, :],
                                   _height[_beg:_end, :],
                                   nc_geo['solar_zenith'][_beg:_end, :],
                                   nc_geo['solar_azimuth'][_beg:_end, :],
                                   nc_geo['sensor_zenith'][_beg:_end, :],
                                   nc_geo['sensor_azimuth'][_beg:_end, :])


def copy_navigation_data(navigation_filename: str, l1a_filename: str) -> None:
    """Copy navigation data from one NetCDF file to another.

    Parameters
    ----------
    navigation_filename
        Path to file containing orbit positions and attitude
        quaternions.
    l1a_filename
        Path to L1A (or L1X) product to which navigation data should
        be appended.

    """
    nc_nav = Dataset(navigation_filename)
    nc_l1a = Dataset(l1a_filename, 'a')
    grp = nc_l1a.createGroup('navigation_data')
    for name, dim in nc_nav.dimensions.items():
        grp.createDimension(name, dim.size)
    for name, var in nc_nav.variables.items():
        grp.createVariable(name, var.datatype, var.dimensions)
        grp[name].setncatts(var.__dict__)
        grp[name][:] = var[:].data


def read_l1(filename: str,
            alt_beg: int,
            alt_end: int | None,
            spectra_in_memory: bool = True,
            geometry_filename: str = '') -> L1:
    """Read a level 1 file.

    Read a list of L1 products from a NetCDF file. The input data
    level may be L1A, L1B, or anything in between. alt_beg/end
    specify a subrange to process.

    Parameters
    ----------
    filename
        Path of L1A file, L1B file, or anything in between.
    alt_beg
        First (zero-based) frame to read, default 0.
    alt_end
        Last frame to read, default last frame in data.
    spectra_in_memory
        Whether to load all spectra into memory. This will run faster
        but might be too memory intensive if the input file has SGM
        line-by-line radiances.
    geometry_filename
        This is where the image timestamps are read.

    Returns
    -------
        L1 product (signal and detector settings).

    """
    if alt_end is not None:
        alt_end += 1
    l1_product = L1.from_empty()
    l1_product.proc_level = read_proc_level(filename)
    nc = Dataset(filename)
    # Attempt to read variables that are present in the NetCDF file
    if 'science_data' in nc.groups:
        grp = nc['science_data']
        l1_product.signal = (
            grp['detector_image'][alt_beg:alt_end, :].astype('f8').data)
        if 'detector_stdev' in grp.variables:
            l1_product.noise = grp['detector_stdev'][alt_beg:alt_end, :].data
        else:
            l1_product.noise = np.ones(l1_product.signal.shape)
        n_alt = l1_product.signal.shape[0]
    elif 'observation_data' in nc.groups:
        grp = nc['observation_data']
        l1_product.wavelengths = grp['wavelength'][:].data
        l1_product.spectra = grp['radiance'][alt_beg:alt_end, :, :].data
        if 'radiance_stdev' in grp.variables:
            l1_product.spectra_noise = grp['radiance_stdev'][
                alt_beg:alt_end, :, :].data
        else:
            l1_product.spectra_noise = np.ones(l1_product.spectra.shape)
        n_alt = l1_product.spectra.shape[0]
    elif 'radiance' in nc.variables:
        l1_product.wavelengths = nc['wavelength'][:].data
        if spectra_in_memory:
            l1_product.spectra = nc['radiance'][alt_beg:alt_end, :, :].data
        else:
            # Dummy slice just to extract the ALT and ACT dimension
            l1_product.spectra = nc['radiance'][alt_beg:alt_end, :, :0].data
        n_alt = l1_product.spectra.shape[0]
    if 'solar_irradiance' in nc.variables:
        l1_product.solar_irradiance = nc['solar_irradiance'][:].data
    if 'image_attributes' in nc.groups:
        grp = nc['image_attributes']
        l1_product.navigation = Navigation.from_shape((n_alt,))
        l1_product.navigation.time = grp['time'][alt_beg:alt_end].data
        l1_product.time_units = grp['time'].units
        l1_product.tai_seconds = grp['tai_seconds'][alt_beg:alt_end].data
        l1_product.tai_subsec = (
            grp['tai_subsec'][alt_beg:alt_end].data / 65535.0)
        l1_product.binning_table_id = grp['binning_table'][:]
        l1_product.coad_factor = grp['nr_coadditions'][:]
        l1_product.exposure_time = grp['exposure_time'][:]
    elif geometry_filename:
        nc_geo = Dataset(geometry_filename)
        l1_product.navigation = Navigation.from_shape((n_alt,))
        l1_product.navigation.time = nc_geo['time'][alt_beg:alt_end].data
        l1_product.time_units = nc_geo['time'].units
        l1_product.tai_seconds = nc_geo['tai_seconds'][alt_beg:alt_end].data
        l1_product.tai_subsec = (
            nc_geo['tai_subsec'][alt_beg:alt_end].data / 65535.0)
    elif l1_product.proc_level < ProcLevel.l1b:
        log.warning('L1 product has no image timestamps because no geometry '
                    'file given')
    if 'navigation_data' in nc.groups:
        grp = nc['navigation_data']
        n_times = grp.dimensions['time'].size
        platform_navigation = Navigation.from_shape((n_times,))
        platform_navigation.time = grp['time'][:].data
        platform_navigation.orb_pos = grp['orb_pos'][:].data
        for i_alt in range(n_times):
            q = grp['att_quat'][i_alt, :].data
            platform_navigation.att_quat[i_alt] = Quaternion(
                q[3], *q[:3]).normalised
        l1_product.navigation = interpolate_navigation_data(
            platform_navigation, l1_product.navigation.time)
    # Reading of main variables done. Now perform a few sanity checks.
    if n_alt == 0:
        raise SystemExit(f"error: data slice [{alt_beg}:{alt_end}] is empty")
    return l1_product


def write_l1(filename: str,
             config: dict,
             l1_product: L1) -> None:
    """Write a L1 file.

    Parameters
    ----------
    filename
        Path of an L1A file, L1B file, or anything in between.
    config
        Configuration parameters.
    l1_product
        L1 product (signal and detector settings).
    geometry
        Include geometry data if available and data consists of spectra.

    """
    default_fill_value = -32767
    out = Dataset(filename, 'w')
    out.Conventions = "CF-1.11"
    out.project = 'TANGO'
    out.instrument = config['instrument']
    out.creator_name = 'SRON/Earth Science'
    out.creator_url = 'https://earth.sron.nl/project/tango'
    out.processing_version = config['processing_version']
    out.product_name = Path(filename).name
    out.date_created = datetime.now(timezone.utc).isoformat(timespec='seconds')
    if not platform.system() == 'Windows':
        out.git_commit = get_git_commit_hash()
    # out.history = ""
    # In C++ code for 'title', config['project'] is fixed to
    # 'Tango' and config['instrument'] is fixed to 'Carbon'
    if l1_product.proc_level in (
            ProcLevel.l1a, ProcLevel.l1b, ProcLevel.sgm):
        if l1_product.proc_level == ProcLevel.l1a:
            out.title = f'TANGO {config["instrument"]} level 1A data'
        elif l1_product.proc_level == ProcLevel.l1b:
            out.title = f'TANGO {config["instrument"]} level 1B data'
        elif l1_product.proc_level == ProcLevel.sgm:
            out.title = f'TANGO {config["instrument"]} SGM radiometric scene'
        out.product_type = str(l1_product.proc_level)
    else:
        out.title = f'TANGO {config["instrument"]} level 1X data'
        out.product_type = "L1X"
        out.l1x_level = l1_product.proc_level.value
        # not in C++ code
        out.l1x_level_name = str(l1_product.proc_level)
    var_config = out.createVariable('configuration', str)
    config_text = (
        yaml.dump(config).replace(' true', ' yes').replace(' false', ' no'))
    var_config[:] = np.array([config_text], dtype='object')
    var_config.comment = 'configuration parameters used to produce this file'
    if l1_product.proc_level <= ProcLevel.stray:
        grp_data = out.createGroup('science_data')
        dim_alt = out.createDimension(
            'along_track_sample', l1_product.signal.shape[0])
        dim_bins = out.createDimension('bin', l1_product.signal.shape[1])
    else:
        dim_alt = out.createDimension('along_track_sample',
                                      l1_product.spectra.shape[0])
        dim_act = out.createDimension('across_track_sample',
                                      l1_product.spectra.shape[1])
        dim_waves = out.createDimension('wavelength',
                                        l1_product.spectra.shape[2])
        if l1_product.proc_level < ProcLevel.sgm:
            grp_data = out.createGroup('observation_data')
    if l1_product.proc_level < ProcLevel.l1b:
        grp_attr = out.createGroup('image_attributes')
        var_timestamps = grp_attr.createVariable(
            'time', 'f8', (dim_alt,), fill_value=default_fill_value)
        var_timestamps.long_name = 'detector image time'
        var_timestamps.units = l1_product.time_units
        var_timestamps.valid_min = 0.0
        var_timestamps.valid_max = 172800.0  # 2 x day
        var_timestamps[:] = l1_product.navigation.time
        var_tai_seconds = grp_attr.createVariable(
            'tai_seconds', 'u4', (dim_alt,), fill_value=0)
        var_tai_seconds.long_name = 'detector image TAI time (seconds)'
        var_tai_seconds.units = 'seconds since 1958-01-01 00:00:00 TAI'
        var_tai_seconds.valid_min = np.uint(1956528000)
        var_tai_seconds.valid_max = np.uint(2493072000)
        var_tai_seconds[:] = l1_product.tai_seconds
        var_tai_subsec = grp_attr.createVariable(
            'tai_subsec', 'f8', (dim_alt,), fill_value=default_fill_value)
        var_tai_subsec.long_name = 'detector image TAI time (subseconds)'
        var_tai_subsec.units = '1/65536'
        var_tai_subsec.valid_min = 0
        var_tai_subsec.valid_max = 65535
        var_tai_subsec[:] = l1_product.tai_subsec * 65535
        var_binning_table = grp_attr.createVariable('binning_table', 'i1')
        var_binning_table[:] = l1_product.binning_table_id
        var_binning_table.long_name = 'binning table ID'
        var_coad_factor = grp_attr.createVariable('nr_coadditions', 'u2')
        var_coad_factor[:] = l1_product.coad_factor
        var_coad_factor.long_name = 'coaddition factor'
        var_coad_factor.comment = 'number of detector read-outs summed'
        var_exposure_time = grp_attr.createVariable(
            'exposure_time', 'f8', fill_value=default_fill_value)
        var_exposure_time.long_name = 'exposure time'
        var_exposure_time.units = 's'
        var_exposure_time.comment = 'exposure time per detector read-out'
        var_exposure_time[:] = l1_product.exposure_time
    if l1_product.proc_level == ProcLevel.l1a:
        var_signal = grp_data.createVariable('detector_image',
                                             'i4',
                                             (dim_alt, dim_bins),
                                             compression='zlib',
                                             complevel=5,
                                             chunksizes=[1, len(dim_bins)],
                                             fill_value=default_fill_value)
        var_signal.long_name = 'image'
        var_signal.units = 'counts'
        var_signal.valid_min = 0
        var_signal.valid_max = 60000
        var_signal[:] = l1_product.signal.astype(int)
    elif l1_product.proc_level <= ProcLevel.stray:
        var_signal = grp_data.createVariable('detector_image',
                                             'f8',
                                             (dim_alt, dim_bins),
                                             compression='zlib',
                                             fill_value=default_fill_value)
        var_signal.long_name = 'signal of each detector pixel'
        var_signal.units = 'counts'
        var_signal.valid_min = -1e100
        var_signal.valid_max = 1e100
        var_signal[:] = l1_product.signal
        var_noise = grp_data.createVariable('detector_stdev',
                                            'f8',
                                            (dim_alt, dim_bins),
                                            compression='zlib',
                                            fill_value=default_fill_value)
        var_noise.long_name = 'standard deviation of detector bin'
        var_noise.units = 'counts'
        var_noise.valid_min = 0.0
        var_noise.valid_max = 1e100
        if l1_product.noise.size == 0:
            l1_product.noise = np.ones(l1_product.signal.shape)
        var_noise[:] = l1_product.noise
    else:
        var_waves = grp_data.createVariable('wavelength', 'f8',
                                            dim_waves,
                                            compression='zlib',
                                            fill_value=default_fill_value)
        var_waves.long_name = 'wavelength'
        var_waves.units = 'nm'
        var_waves.valid_min = 0.0
        var_waves.valid_max = 8000.0
        var_waves[:] = l1_product.wavelengths
        var_spectra = grp_data.createVariable(
            'radiance',
            'f8',
            (dim_alt, dim_act, dim_waves),
            compression='zlib',
            fill_value=default_fill_value)
        var_spectra[:] = l1_product.spectra
        if l1_product.proc_level > ProcLevel.swath:
            var_spectra.long_name = 'spectral photon radiance'
            var_spectra.units = 'nm-1 s-1 sr-1 m-2'
        else:
            var_spectra.long_name = 'signal of each detector pixel'
            var_spectra.units = 'counts'
        var_spectra.valid_min = 0.0
        var_spectra.valid_max = 1e20
        var_noise = grp_data.createVariable(
            'radiance_stdev',
            'f8',
            (dim_alt, dim_act, dim_waves),
            compression='zlib',
            fill_value=default_fill_value)
        var_noise.long_name = 'standard deviation of radiance in bin'
        var_noise.units = 'nm-1 s-1 sr-1 m-2'
        var_noise.valid_min = 0.0
        var_noise.valid_max = 1e20
        if l1_product.spectra_noise.size == 0:
            l1_product.spectra_noise = np.ones(l1_product.spectra.shape)
        var_noise[:] = l1_product.spectra_noise
    if l1_product.geometry.lat.size > 0:
        grp_geo = out.createGroup('geolocation_data')
        # If doing geolocation but not calibrating up to L1B
        if 'across_track_sample' not in out.dimensions:
            dim_act = out.createDimension('across_track_sample',
                                          l1_product.geometry.lat.shape[1])
        nc_write_geometry(grp_geo, l1_product.geometry)
    out.close()
