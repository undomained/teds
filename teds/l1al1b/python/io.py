# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Functions dealing with input/output.

These are functions for formatting the output to stdout, for
reading/writing input/output files, and tests for things like whether
a file is readable/writable.

"""
from datetime import datetime
from datetime import timezone
from pathlib import Path
import os

from netCDF4 import Dataset
from netCDF4 import default_fillvals
import numpy as np
import numpy.typing as npt
import subprocess
import xarray as xr
import yaml

from .types import CKD
from .types import Geometry
from .types import L1
from .types import ProcLevel


def print_heading(heading: str, empty_line: bool = True) -> None:
    """Print the name of a processing section.

    For example, reading the CKD could start with
    ######################
    # CKD initialization #
    ######################

    This is meant to be called directly as opposed to using the
    logger. Things like timestamps are not necessary for just printing
    a section heading.

    Parameters
    ----------
    heading
        Title of the section to be displayed.
    empty_line
        Whether to print an empty line before printing the heading.

    """
    if empty_line:
        print()
    print('#' * (len(heading) + 4))
    print(f'# {heading} #')
    print('#' * (len(heading) + 4))


def get_git_commit_hash() -> str:
    """Return short git hash if .git found"""
    git_hash = subprocess.run(
        ['git', 'rev-parse', 'HEAD'],
        shell=False,
        capture_output=True,
        cwd=os.path.dirname(__file__)).stdout.decode('utf-8')
    if git_hash:
        return git_hash[:8]
    return ''


def print_system_info() -> None:
    """Print information about the host system and some runtime
    options.

    """
    # Project version
    _version_file = os.path.join(os.path.dirname(__file__),
                                 '../../../project_version.txt')
    with open(_version_file) as f:
        for line in f.readlines():
            if not line.startswith('#'):
                version = line.rstrip()
    print('Version                 :', version)
    # Short git hash (only if .git found)
    git_hash = get_git_commit_hash()
    if git_hash:
        print('Commit hash             :', git_hash)
    # Datetime and contacts
    print('Date and timezone       :', datetime.now().strftime('%Y %B %d %a'))
    print('Contacts                : raullaasner@gmail.com')
    print('                          '
          'bitbucket.org/sron_earth/teds/issues (request permission)')
    # Platform
    host_system = subprocess.run(
        ['uname', '-sr'], shell=False, capture_output=True).stdout
    print('Host system             :',
          host_system.decode('utf-8').split('\n')[0])


def read_binning_pattern(binning_file: str,
                         binning_table_id: int,
                         n_detrows: int,
                         n_detcols: int) -> tuple[npt.NDArray[np.int32],
                                                  npt.NDArray[np.int32]]:
    """Read a binning pattern.

    Parameters
    ----------
    binning_file
        Path of file with binning patterns.
    binning_table_id
        Identifier of binning pattern, 0 is unbinned.
    n_detrows
        Number of detector rows.
    n_detcols
        Number of detector columns.

    Returns
    -------
        Two-dimensional (unbinned) detector map of bin indices and
        one-dimensional list of number of detector pixels per bin.

    """
    if binning_table_id == 0:
        bin_indices = np.arange(n_detrows * n_detcols, dtype=np.int32).reshape(
            n_detrows, n_detcols)
        count_table = np.ones(n_detrows * n_detcols, np.int32)
    else:
        if not Path(binning_file).is_file():
            binning_text = ("binning file" if binning_file is None
                            else binning_file)
            raise SystemExit(f"ERROR: {binning_text} not found")
        ckd_binning = xr.load_dataset(binning_file,
                                      group=f'Table_{binning_table_id}')
        bin_indices = ckd_binning['binning_table'].values
        uniq_bin_indices, tmp_count_table = np.unique(bin_indices,
                                                      return_counts=True)
        count_table = np.zeros(np.max(uniq_bin_indices) + 1, np.int32)
        count_table[uniq_bin_indices] = tmp_count_table
    return bin_indices, count_table


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
    ckd: CKD = {}
    ckd['n_detrows'] = nc.dimensions['detector_row'].size
    ckd['n_detcols'] = nc.dimensions['detector_column'].size
    # False good, True bad
    ckd['pixel_mask'] = nc['pixel_mask'][:].astype(bool)
    if 'dark' in nc.groups:
        ckd['dark'] = {
            'offset': nc['dark/offset'][:],  # counts
            'current': nc['dark/current'][:],  # counts/s
        }
    if 'noise' in nc.groups:
        # Conversion_gain can be a detector map or 1 value, optionally
        # dependent on temperature. Best practice is to use the same
        # value for all pixels and move a dependence to the PRNU CKD.
        ckd['noise'] = {
            'conversion_gain': 1 / nc['noise/g'][:],  # e/counts
            'read_noise': nc['noise/n'][:],  # counts
        }
    if 'nonlinearity' in nc.groups:
        ckd['nonlin'] = {
            'expected': nc['nonlinearity/y'][:],  # counts
            'observed': nc['nonlinearity/knots'][:],  # counts
        }
        if (
                not (np.diff(ckd['nonlin']['expected']) > 0).all()
                or not (np.diff(ckd['nonlin']['observed']) > 0).all()):
            raise SystemExit(f"ERROR: nonlinearity data in {filename} are not "
                             "strictly increasing")
    if 'prnu' in nc.groups:
        # PRNU including quantum efficiency
        ckd['prnu'] = {
            'prnu_qe': nc['prnu/prnu'][:]
        }
        # Not in C++ code
        ckd['prnu']['prnu_qe'][ckd['prnu']['prnu_qe'] <= 0] = np.nan
    if 'stray' in nc.groups:
        ckd['stray'] = {}
        ckd['stray']['kernels_fft'] = []
        n_kernels = nc['stray'].dimensions['kernel'].size
        kernel_n_rows = nc['stray/kernel_rows'][:]
        kernel_n_cols = nc['stray/kernel_cols'][:]
        kernel_fft_sizes = nc['stray/kernel_fft_sizes'][:]
        kernels_fft = nc['stray/kernels_fft'][:]
        for i_kernel in range(n_kernels):
            global_beg = 2 * sum(kernel_fft_sizes[0:i_kernel])
            global_end = global_beg + 2 * kernel_fft_sizes[i_kernel]
            kernel_fft = kernels_fft[global_beg:global_end]
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
            ckd['stray']['kernels_fft'].append(fft_full)
        # Internal scattering factor
        ckd['stray']['eta'] = nc['stray/eta'][:]
        # Margins between kernels and detector array: bottom, top,
        # left, right of kernel
        ckd['stray']['edges'] = nc['stray/edges'][:].astype('i4')
        ckd['stray']['weights'] = nc['stray/weights'][:]
    if 'radiometric' in nc.groups:
        # nm-1 sr-1 m-2 e count-1
        ckd['radiometric'] = {
            'rad_corr': nc['radiometric/radiometric'][:],
        }
    if 'spectral' in nc.groups:
        # Smile (spectral distortion) description:
        # wavelength at each detector pixel
        ckd['spectral'] = {
            'wavelengths': nc['spectral/wavelength'][:],  # nm
        }
    if 'swath' in nc.groups:
        # not yet used
        ckd['swath'] = {
            'line_of_sights': nc['swath/line_of_sight'][:],
            # keystone (spatial distortion) description: for each ACT
            # ground pixel (spectrum), the row index (float) at each
            # column is given.
            'spectrum_rows': nc['swath/row_index'][:],
        }
        if not monotonic(ckd['swath']['spectrum_rows'], axis=0):
            raise SystemExit(f"ERROR: swath/row_index in {filename} is not "
                             "monotonic")
    return ckd


def monotonic(x: npt.NDArray[np.float64], axis: int = -1) -> np.bool:
    """Test whether an array is strictly increasing/decreasing or not."""
    dx = np.diff(x, axis=axis)
    return np.all(dx < 0) or np.all(dx > 0)


def read_geometry(l1_products: L1, config: dict) -> Geometry:
    """Read geometry from an original geometry file or an L1 file.

    Geometry is not used for calculations in L1 processing. The
    dictionaries are only used to read some settings.

    Parameters
    ----------
    l1_products
        L1 products (signal and detector settings).
    config
        Configuration parameters.

    Returns
    -------
        Coordinates and angles.

    """
    filename = config['io']['geometry']
    original_image_start = l1_products['original_image_start']
    original_image_end = l1_products['original_image_end']+1
    image_start = config['image_start']
    image_end = None if config['image_end'] is None else config['image_end']+1
    with Dataset(filename) as root:
        groups = list(root.groups)
        variables = list(root.variables)
    geo: Geometry = {}
    if 'geolocation_data' in groups:  # L1 file
        # Assumption: the input file currently processed or a file
        # with the same geometry data
        nc = Dataset(filename)
        grp = nc['geolocation_data']
        geo['latitude'] = grp['latitude'][image_start:image_end, :]
        geo['longitude'] = grp['longitude'][image_start:image_end, :]
        geo['height'] = grp['height'][image_start:image_end, :]
        geo['saa'] = grp['solar_azimuth'][image_start:image_end, :]
        geo['sza'] = grp['solar_zenith'][image_start:image_end, :]
        geo['vaa'] = grp['sensor_azimuth'][image_start:image_end, :]
        geo['vza'] = grp['sensor_zenith'][image_start:image_end, :]
    elif all(x in variables for x
             in ['lat', 'lon', 'saa', 'sza', 'vaa', 'vza']):  # geometry file
        # Assumption: the geometry file is consistent with the data in
        # l1_data['original_file']
        nc = Dataset(filename)
        _beg = original_image_start
        _end = original_image_end
        geo['latitude'] = nc['lat'][_beg:_end, :]
        geo['longitude'] = nc['lon'][_beg:_end, :]
        geo['saa'] = nc['saa'][_beg:_end, :]
        geo['sza'] = nc['sza'][_beg:_end, :]
        geo['vaa'] = nc['vaa'][_beg:_end, :]
        geo['vza'] = nc['vza'][_beg:_end, :]
        if 'height' in groups:
            geo['height'] = nc['height'][_beg:_end, :]
        else:
            geo['height'] = np.zeros_like(geo['latitude'])
    return geo


def read_l1(filename: str,
            image_start: int,
            image_end: int | None = None) -> L1:
    """Read a level 1 file.

    Read a list of L1 products from a NetCDF file. The input data
    level may be L1A, L1B, or anything in between. image_start/end
    specify a subrange to process.

    Parameters
    ----------
    filename
        Path of L1A file, L1B file, or anything in between.
    image_start
        First (zero-based) frame to read, default 0.
    image_end
        Last frame to read, default last frame in data.

    Returns
    -------
        L1 products (signal and detector settings).

    """
    image_end = None if image_end is None else image_end+1
    proc_level = read_proc_level(filename)
    l1_products: L1 = {}
    l1_products['proc_level'] = proc_level
    nc = Dataset(filename)
    if proc_level <= ProcLevel.stray:
        grp = nc['science_data']
        l1_products['signal'] = (
            grp['detector_image'][image_start:image_end, :].astype('f8'))
        if proc_level > ProcLevel.l1a:
            l1_products['noise'] = (
                grp['detector_stdev'][image_start:image_end, :])
        n_data = l1_products['signal'].shape[0]
    elif proc_level < ProcLevel.sgm:
        grp = nc['observation_data']
        l1_products['spectra'] = grp['radiance'][image_start:image_end, :, :]
        l1_products['spectra_noise'] = grp['radiance_stdev'][
            image_start:image_end, :, :]
        l1_products['wavelengths'] = grp['wavelength'][:]
        n_data = l1_products['spectra'].shape[0]
        if not monotonic(l1_products['wavelengths']):
            raise SystemExit("ERROR: observation_data/wavelength in "
                             f"{filename} is not monotonic")
    else:  # SGM data
        l1_products['spectra'] = nc['radiance'][image_start:image_end, :, :]
        l1_products['wavelengths'] = nc['wavelength'][:]
        try:
            l1_products['solar_irradiance'] = nc['solar_irradiance'][:]
        except KeyError:
            l1_products['solar_irradiance'] = np.array(0, dtype=np.float64)
        # C++ code fills with 1
        l1_products['spectra_noise'] = np.full_like(l1_products['spectra'],
                                                    np.nan)
        n_data = l1_products['spectra'].shape[0]
        if not monotonic(l1_products['wavelengths']):
            raise SystemExit("ERROR: observation_data/wavelength in "
                             f"{filename} is not monotonic")
    if n_data == 0:
        raise SystemExit(f"ERROR: data slice [{image_start}:{image_end}] is "
                         "empty")
    try:
        grp = nc['image_attributes']
        l1_products['timestamps'] = grp['image_time'][image_start:image_end]
        l1_products['binning_table_ids'] = (
            grp['binning_table'][image_start:image_end])
        l1_products['coad_factors'] = (
            grp['nr_coadditions'][image_start:image_end])
        l1_products['exptimes'] = grp['exposure_time'][image_start:image_end]
    except IndexError:  # SGM/L1B data
        # Apparently not relevant for L1B
        l1_products['timestamps'] = np.zeros(n_data)
        # Assume data is unbinned if not specified
        l1_products['binning_table_ids'] = np.zeros(n_data)
        l1_products['coad_factors'] = np.ones(n_data)
        l1_products['exptimes'] = np.full(n_data, np.nan)
    if len(np.unique(l1_products['binning_table_ids'])) > 1:
        raise SystemExit("ERROR: different binning tables in data set, use "
                         "subset")
    # Some attributes that are not in C++ code. The first and last
    # frame of the original file used during the current processing.
    try:
        # Currently processed data is not the original data, but has
        # been written by this code before.
        l1_products['original_file'] = nc.original_file
        l1_products['original_image_start'] = (
            nc.original_image_start + image_start)
        l1_products['original_image_end'] = (
            nc.original_image_start + image_start + n_data - 1)
    except AttributeError:
        if 'configuration' in list(nc.variables):
            # Currently processed data is not the original data, but
            # has been written by the C++ code. If the original data
            # was subsetted more than 1 processing earlier, the values
            # are wrong.
            previous_config = yaml.safe_load(str(nc['configuration'][:]))
            l1_products['original_file'] = ''
            l1_products['original_image_start'] = (
                previous_config['image_start']+image_start)
            l1_products['original_image_end'] = (
                previous_config['image_start']+image_start+n_data-1)
        else:
            # currently processed data is the original data
            l1_products['original_file'] = filename
            l1_products['original_image_start'] = image_start
            l1_products['original_image_end'] = image_start+n_data-1
    return l1_products


def write_l1(filename: str,
             config: dict,
             l1_products: L1,
             geometry: bool = False) -> str:
    """Write a L1 file.

    Parameters
    ----------
    filename
        Path of an L1A file, L1B file, or anything in between.
    config
        Configuration parameters.
    l1_products
        L1 products (signal and detector settings).
    geometry
        Include geometry data if available and data consists of spectra.

    Returns
    -------
        Name of the NetCDF created

    """
    # Fill values are set explicitly to the netCDF default where the
    # C++ code always uses -32767.
    out = Dataset(filename, 'w')
    out.Conventions = "CF-1.11"
    out.project = 'TANGO'
    out.instrument = config['instrument']
    out.creator_name = 'SRON/Earth Science'
    out.creator_url = 'https://earth.sron.nl/project/tango'
    out.processing_version = config['processing_version']
    out.product_name = Path(filename).name
    out.date_created = datetime.now(timezone.utc).isoformat(timespec='seconds')
    # Not in C++ code
    out.original_file = l1_products['original_file']
    # Not in C++ code
    out.original_image_start = l1_products['original_image_start']
    # Not in C++ code
    out.original_image_end = l1_products['original_image_end']
    out.git_commit = get_git_commit_hash()
    # out.history = ""
    # In C++ code for 'title', config['project'] is fixed to
    # 'Tango' and config['instrument'] is fixed to 'Carbon'
    if l1_products['proc_level'] in (
            ProcLevel.l1a, ProcLevel.l1b, ProcLevel.sgm):
        if l1_products['proc_level'] == ProcLevel.l1a:
            out.title = f'TANGO {config["instrument"]} level 1A data'
        elif l1_products['proc_level'] == ProcLevel.l1b:
            out.title = f'TANGO {config["instrument"]} level 1B data'
        elif l1_products['proc_level'] == ProcLevel.sgm:
            out.title = f'TANGO {config["instrument"]} SGM radiometric scene'
        out.product_type = str(l1_products['proc_level'])
    else:
        out.title = f'TANGO {config["instrument"]} level 1X data'
        out.product_type = "L1X"
        out.l1x_level = l1_products['proc_level'].value
        # not in C++ code
        out.l1x_level_name = str(l1_products['proc_level'])
    var_config = out.createVariable('configuration', str)
    config_text = (
        yaml.dump(config).replace(' true', ' yes').replace(' false', ' no'))
    var_config[:] = np.array([config_text], dtype='object')
    var_config.comment = 'configuration parameters used to produce this file'
    # A variable should only be given the same name as a dimension in
    # a netCDF file when it is to be used as a coordinate variable.
    if l1_products['proc_level'] <= ProcLevel.stray:
        grp_data = out.createGroup('science_data')
        dim_data = out.createDimension('image', l1_products['signal'].shape[0])
        dim_bins = out.createDimension('bin', l1_products['signal'].shape[1])
    else:
        dim_data = out.createDimension('along_track_sample',
                                       l1_products['spectra'].shape[0])
        dim_act = out.createDimension('across_track_sample',
                                      l1_products['spectra'].shape[1])
        # 'wavelength' is also a variable, suggest 'spectral_sample'
        dim_waves = out.createDimension('wavelength',
                                        l1_products['spectra'].shape[2])
        if l1_products['proc_level'] < ProcLevel.sgm:
            grp_data = out.createGroup('observation_data')
    if l1_products['proc_level'] < ProcLevel.l1b:
        grp_attr = out.createGroup('image_attributes')
        var_timestamps = grp_attr.createVariable(
            'image_time',
            'f8',
            (dim_data,),
            compression='zlib',
            fill_value=default_fillvals['f8'])
        var_timestamps[:] = np.nan_to_num(l1_products['timestamps'],
                                          nan=default_fillvals['f8'])
        var_timestamps.long_name = 'time'  # C++ code has 'image time'
        var_timestamps.units = 'seconds since 2022-03-21'  # weird date
        var_timestamps.valid_min = 0.0
        var_timestamps.valid_max = 92304.0  # weird value
        var_binning_table_ids = grp_attr.createVariable(
            'binning_table',
            'i1',
            (dim_data,),
            compression='zlib',
            fill_value=default_fillvals['i1'])
        var_binning_table_ids[:] = np.nan_to_num(
            l1_products['binning_table_ids'], nan=default_fillvals['i1'])
        var_binning_table_ids.long_name = 'binning table ID'
        var_coad_factors = grp_attr.createVariable(
            'nr_coadditions',
            'u2',
            (dim_data,),
            compression='zlib',
            fill_value=default_fillvals['u2'])
        var_coad_factors[:] = np.nan_to_num(
            l1_products['coad_factors'], nan=default_fillvals['u2'])
        var_coad_factors.long_name = 'coaddition factor'
        var_coad_factors.comment = 'number of detector read-outs summed'
        var_exptimes = grp_attr.createVariable(
            'exposure_time',
            'f8',
            (dim_data,),
            compression='zlib',
            fill_value=default_fillvals['f8'])
        var_exptimes[:] = np.nan_to_num(
            l1_products['exptimes'], nan=default_fillvals['f8'])
        var_exptimes.long_name = 'exposure time'
        var_exptimes.units = 's'
        var_exptimes.comment = 'exposure time per detector read-out'
    if l1_products['proc_level'] == ProcLevel.l1a:
        # In C++ code only L1A signal data are compressed
        var_signal = grp_data.createVariable('detector_image',
                                             'i4',
                                             (dim_data, dim_bins),
                                             compression='zlib',
                                             complevel=5,
                                             chunksizes=[1, len(dim_bins)],
                                             fill_value=default_fillvals['i4'])
        var_signal[:] = np.nan_to_num(
            l1_products['signal'], nan=default_fillvals['i4']).astype(int)
        var_signal.long_name = 'signal'
        var_signal.units = 'counts'
        var_signal.valid_min = 0
        var_signal.valid_max = 60000
    elif l1_products['proc_level'] <= ProcLevel.stray:
        # C++ code does not use fill values
        var_signal = grp_data.createVariable('detector_image',
                                             'f8',
                                             (dim_data, dim_bins),
                                             compression='zlib',
                                             fill_value=default_fillvals['f8'])
        var_signal[:] = np.nan_to_num(
            l1_products['signal'], nan=default_fillvals['f8'])
        var_signal.long_name = 'signal'  # C++ code has 'detector images'
        var_signal.units = 'counts'
        var_signal.valid_min = -1e100
        var_signal.valid_max = 1e100
        var_noise = grp_data.createVariable('detector_stdev',
                                            'f8',
                                            (dim_data, dim_bins),
                                            compression='zlib',
                                            fill_value=default_fillvals['f8'])
        var_noise[:] = np.nan_to_num(
            l1_products['noise'], nan=default_fillvals['f8'])
        if l1_products['proc_level'] != ProcLevel.raw:
            # C++ code has 'standard deviation of detector bin'
            var_noise.long_name = 'noise'
            var_noise.units = 'counts'
        else:  # not in C++ code
            var_noise.long_name = 'binning factor'
        var_noise.valid_min = 0.0
        var_noise.valid_max = 1e100
    elif l1_products['proc_level'] < ProcLevel.sgm:
        var_waves = grp_data.createVariable('wavelength', 'f4',
                                            (dim_act, dim_waves),
                                            compression='zlib',
                                            fill_value=default_fillvals['f4'])
        var_waves[:] = np.nan_to_num(
            l1_products['wavelengths'], nan=default_fillvals['f4'])
        var_waves.long_name = 'wavelength'
        var_waves.units = 'nm'
        var_waves.valid_min = 0.0
        var_waves.valid_max = 999.0
        is_l1b = l1_products['proc_level'] >= ProcLevel.l1b
        ftype = 'f4' if is_l1b else 'f8'
        var_spectra = grp_data.createVariable(
            'radiance',
            ftype,
            (dim_data, dim_act, dim_waves),
            compression='zlib',
            fill_value=default_fillvals[ftype])
        var_spectra[:] = np.nan_to_num(
            l1_products['spectra'], nan=default_fillvals[ftype])
        if l1_products['proc_level'] > ProcLevel.swath:
            var_spectra.long_name = 'spectral photon radiance'
            var_spectra.units = 'nm-1 s-1 sr-1 m-2'
        else:  # not in C++ code
            var_spectra.long_name = 'signal'
            var_spectra.units = 'counts'
        var_spectra.valid_min = (
            np.array(0, ftype) if is_l1b else np.array(-1e20, ftype))
        var_spectra.valid_max = np.array(1e20, ftype)
        var_noise = grp_data.createVariable(
            'radiance_stdev',
            ftype,
            (dim_data, dim_act, dim_waves),
            compression='zlib',
            fill_value=default_fillvals[ftype])
        var_noise[:] = np.nan_to_num(
            # C++ code writes 1s
            l1_products['spectra_noise'], nan=default_fillvals[ftype])
        # C++ code has 'standard deviation of radiance in bin'
        var_noise.long_name = 'spectral photon radiance noise'
        var_noise.units = 'nm-1 s-1 sr-1 m-2'
        var_noise.valid_min = (
            np.array(0, ftype) if is_l1b else np.array(-1e20, ftype))
        var_noise.valid_max = np.array(1e20, ftype)
    else:  # SGM data (subset), C++ does not give this option
        var_waves = out.createVariable(
            'wavelength',
            'f8',
            dim_waves,
            compression='zlib',
            fill_value=default_fillvals['f8'])
        var_waves[:] = np.nan_to_num(
            l1_products['wavelengths'], nan=default_fillvals['f8'])
        var_waves.long_name = 'wavelength'
        var_waves.units = 'nm'
        var_waves.valid_min = 0.0
        var_waves.valid_max = 8000.0
        ftype = 'f8'
        var_spectra = out.createVariable(
            'radiance',
            ftype,
            (dim_data, dim_act, dim_waves),
            compression='zlib',
            fill_value=default_fillvals[ftype])
        var_spectra[:] = np.nan_to_num(
            l1_products['spectra'], nan=default_fillvals[ftype])
        # C++ code has 'line-by-line radiance'
        var_spectra.long_name = 'spectral photon radiance'
        # C++ code has 'ph nm-1 s-1 sr-1 m-2'
        var_spectra.units = 'nm-1 s-1 sr-1 m-2'
        var_spectra.valid_min = np.array(0, ftype)
        var_spectra.valid_max = np.array(1e28, ftype)
        if l1_products['solar_irradiance']:
            var_irr = out.createVariable(
                'solar_irradiance',
                ftype,
                dim_waves,
                compression='zlib',
                fill_value=default_fillvals[ftype])
            var_irr[:] = np.nan_to_num(l1_products['solar_irradiance'],
                                       nan=default_fillvals[ftype])
            # C++ code has 'line-by-line solar irradiance'
            var_irr.long_name = 'spectral photon irradiance'
            # C++ code has 'ph nm-1 s-1 m-2'
            var_irr.units = 'nm-1 s-1 m-2'
            var_irr.valid_min = np.array(0, ftype)
            var_irr.valid_max = np.array(1e30, ftype)
    if geometry and (l1_products['proc_level'] > ProcLevel.stray):
        # C++ code writes 0s in data towards L1A.
        geo: Geometry = read_geometry(l1_products, config)
        grp_geo = out.createGroup('geolocation_data')
        var_latitude = grp_geo.createVariable(
            'latitude',
            'f4',
            (dim_data, dim_act),
            compression='zlib',
            fill_value=default_fillvals['f4'])
        var_latitude[:] = np.nan_to_num(geo['latitude'],
                                        nan=default_fillvals['f4'])
        # C++ code has 'latitude at bin locations'
        var_latitude.long_name = 'latitude'
        var_latitude.units = 'degrees_north'
        var_latitude.valid_min = -90.0
        var_latitude.valid_max = 90.0
        var_longitude = grp_geo.createVariable(
            'longitude',
            'f4',
            (dim_data, dim_act),
            compression='zlib',
            fill_value=default_fillvals['f4'])
        var_longitude[:] = np.nan_to_num(geo['longitude'],
                                         nan=default_fillvals['f4'])
        # C++ code has 'longitude at bin locations'
        var_longitude.long_name = 'longitude'
        var_longitude.units = 'degrees_east'
        var_longitude.valid_min = -180.0
        var_longitude.valid_max = 180.0
        var_height = grp_geo.createVariable('height',
                                            'f4',
                                            (dim_data, dim_act),
                                            compression='zlib',
                                            fill_value=default_fillvals['f4'])
        var_height[:] = np.nan_to_num(geo['height'],
                                      nan=default_fillvals['f4'])
        # C++ code has 'height at bin locations'
        var_height.long_name = 'height'
        var_height.units = 'm'
        var_height.valid_min = -1000.0
        var_height.valid_max = 10000.0
        var_saa = grp_geo.createVariable('solar_azimuth',
                                         'f4',
                                         (dim_data, dim_act),
                                         compression='zlib',
                                         fill_value=default_fillvals['f4'])
        var_saa[:] = np.nan_to_num(geo['saa'], nan=default_fillvals['f4'])
        # C++ code has 'solar azimuth angle at bin locations'
        var_saa.long_name = 'solar azimuth angle'
        var_saa.units = 'degrees'
        var_saa.valid_min = -180.0
        var_saa.valid_max = 180.0
        var_sza = grp_geo.createVariable('solar_zenith',
                                         'f4',
                                         (dim_data, dim_act),
                                         compression='zlib',
                                         fill_value=default_fillvals['f4'])
        var_sza[:] = np.nan_to_num(geo['sza'], nan=default_fillvals['f4'])
        # C++ code has 'solar zenith angle at bin locations'
        var_sza.long_name = 'solar zenith angle'
        var_sza.units = 'degrees'
        var_sza.valid_min = -90.0
        var_sza.valid_max = 90.0
        var_vaa = grp_geo.createVariable('sensor_azimuth',
                                         'f4',
                                         (dim_data, dim_act),
                                         compression='zlib',
                                         fill_value=default_fillvals['f4'])
        var_vaa[:] = np.nan_to_num(geo['vaa'], nan=default_fillvals['f4'])
        # C++ code has 'sensor azimuth angle at bin locations'
        var_vaa.long_name = 'sensor azimuth angle'
        var_vaa.units = 'degrees'
        var_vaa.valid_min = -180.0
        var_vaa.valid_max = 180.0
        var_vza = grp_geo.createVariable('sensor_zenith',
                                         'f4',
                                         (dim_data, dim_act),
                                         compression='zlib',
                                         fill_value=default_fillvals['f4'])
        var_vza[:] = np.nan_to_num(geo['vza'], nan=default_fillvals['f4'])
        # C++ code has 'sensor zenith angle at bin locations'
        var_vza.long_name = 'sensor zenith angle'
        var_vza.units = 'degrees'
        var_vza.valid_min = -90.0
        var_vza.valid_max = 90.0
    return filename
