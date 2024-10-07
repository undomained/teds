# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Binning and unbinning related operations."""
from scipy.interpolate import RBFInterpolator
import numpy as np
import numpy.typing as npt


def bin_data(data: npt.NDArray[np.float64],
             bin_indices: npt.NDArray[np.int32],
             count_table: npt.NDArray[np.int32]) -> npt.NDArray[np.float64]:
    """Bin each frame in a set of frames and take the sum per bin.

    The sum over all data remains the same. If the set of frames
    contain booleans, an output bin is True if at least one bin value
    is True. One binned frame is a 1D array. An input frame and the
    bin indices can have any number of dimensions, as long as they
    have the same number of elements.

    Parameters
    ----------
    data
        Unbinned frames.
    bin_indices
        Bin index of each pixel in an unbinned frame.
    count_table
        Number of pixels in each bin of a binned frame.

    Returns
    -------
        Binned frames.

    """
    if data.dtype == bool:
        binned_data = np.full((len(data), len(count_table)), False)
        for i, s in zip(bin_indices.ravel(), data.reshape(len(data), -1).T):
            binned_data[:, i] = binned_data[:, i] | s
    else:
        binned_data = np.zeros((len(data), len(count_table)))
        for i, s in zip(bin_indices.ravel(), data.reshape(len(data), -1).T):
            binned_data[:, i] += s
    return binned_data


def unbin_data(binned_data: npt.NDArray[np.float64],
               bin_indices: npt.NDArray[np.int32],
               count_table: npt.NDArray[np.int32],
               method: str = 'nearest') -> npt.NDArray[np.float64]:
    """Unbin each frame in a set of binned frames.

    Bin values are spread over all positions in an unbinned frame
    where the array of bin indices has the corresponding bin index,
    with a tweak depending on the method.

    Parameters
    ----------
    binned_data
        Binned frames.
    bin_indices
        Bin index of each pixel in an unbinned frame.
    count_table
        Number of pixels in each bin of a binned frame.
    method
        Method for unbinning, 'nearest' splits each bin value equally between
        all corresponding pixels (keeping the sum over all data the same),
        'linear' or 'cubic' puts each ratio of bin value and pixel count at the
        centroid of the corresponding pixels and interpolates.

    Returns
    -------
        Unbinned frames.

    """
    # Actually, method can also be 'thin_plate_spline' and 'quintic'.
    # binned_data has dimensions (n_frames, n_bins), even if it
    # contains unbinned data. bin_indices has dimensions (n_rows,
    # n_cols). C++ code makes bad signals 0 but can also copy values
    # from neighbouring pixels when method is linear or cubic.
    scaled_binned_data = binned_data/count_table
    if (
            (bin_indices.size == scaled_binned_data.shape[1])
            and (bin_indices.ravel()
                 == np.arange(scaled_binned_data.shape[1])).all()):
        # data are already unbinned
        data = scaled_binned_data
    elif method == 'nearest':
        data = scaled_binned_data[:, bin_indices]
        data = data.reshape(data.shape[0], -1)
    else:
        # Assuming only groups of directly neighbouring pixels are
        # binned, so a centroid makes sense. C++ code also assumes
        # binning is only along columns with the same binning pattern
        # in all columns. Here a 2D interpolation is used, so
        # absorption lines may get smoothed out.
        # Guess how many neighbouring pixels should be used.
        area = {'linear': 9, 'cubic': 25, 'quintic': 49}.get(method, 25)
        idx_dict: dict = {}
        for pos in np.ndindex(bin_indices.shape):
            bi = bin_indices[pos]
            if np.isfinite(scaled_binned_data[:, bi]).all():
                idx_dict[bi] = idx_dict.get(bi, []) + [pos]
        # unique values in bin_indices or relevant indices of binned data
        bin_idx = np.array(list(idx_dict.keys()))
        # Average detector positions of bins
        bin_centroids = np.array(
            [np.mean(idx_dict[bi], axis=0) for bi in bin_idx])
        det_grid = np.mgrid[:bin_indices.shape[0],
                            :bin_indices.shape[1]].reshape(2, -1).T
        data_list: list = []
        for binned_frame in scaled_binned_data:
            # https://stackoverflow.com/questions/37872171
            rbf = RBFInterpolator(bin_centroids,
                                  binned_frame[bin_idx],
                                  smoothing=0,
                                  neighbors=area,
                                  kernel=method)
            data_list.append(rbf(det_grid))
        data = np.array(data_list)
        # bad pixels are not set to NaN here
        data = data.reshape(data.shape[0], -1)
    return data
