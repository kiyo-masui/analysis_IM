# TODO: consider moving azimuthalAverage and _nd version to the attic
import numpy as np
import math
from core import algebra


def radius_array(input_array):
    """Find the Euclidian distance of all the points in an array using their
    axis meta-data. e.g. x_axis[0]^2 + y_axis[0]^2. (n-dim)
    """
    index_array = np.indices(input_array.shape)
    scale_array = np.zeros(index_array.shape)

    for axis_index in range(input_array.ndim):
        axis_name = input_array.axes[axis_index]
        axis = input_array.get_axis(axis_name)
        scale_array[axis_index, ...] = axis[index_array[axis_index, ...]]

    scale_array = np.rollaxis(scale_array, 0, scale_array.ndim)

    return np.sum(scale_array ** 2., axis=-1) ** 0.5


def suggest_bins(input_array, truncate=True, nbins=40, logbins=True,
                 radius_arr=None):
    """Bin the points in an array by radius

    Parameters
    ----------
    input_array: np.ndarray
        array over which to bin
    truncate: boolean
        maximum radius is the smallest dimension of the array
    nbins: scalar
        number of bins to use if not given some in advance
    logbins: boolean
        generate a log-spaced binning
    radius_arr: np.ndarray
        optional array of |k| to avoid recalculation
    TODO: throw out bins where the counts/bin are too low
    """
    if radius_arr is None:
        radius_arr = radius_array(input_array)

    radius_sorted = np.sort(radius_arr.flat)

    if truncate:
        axis_range = np.zeros((input_array.ndim))
        for axis_index in range(input_array.ndim):
            axis_name = input_array.axes[axis_index]
            axis = input_array.get_axis(axis_name)
            axis_range[axis_index] = axis.max()

        max_r = axis_range.min()
    else:
        max_r = radius_sorted[-1]

    # ignore the bin at k=0 (index 0 in sorted)
    if logbins:
        bins = np.logspace(math.log10(radius_sorted[1]),
                           math.log10(max_r),
                           num=(nbins + 1), endpoint=True)
    else:
        bins = np.linspace(radius_sorted[1], max_r,
                           num=(nbins + 1), endpoint=True)

    print "%d bins from %10.15g to %10.15g" % (nbins, radius_sorted[1], max_r)

    return bins


def bin_edges(bins, log=False):
    """report the bin edges and centers using the same convention as
    np.histogram. `log` reports the log-center
    """
    bin_left, bin_right = bins[:-1], bins[1:]
    if log:
        bin_center = 10 ** (0.5 * (np.log10(bin_left) +
                               np.log10(bin_right)))
    else:
        bin_center = 0.5 * (bin_left + bin_right)

    return bin_left, bin_center, bin_right


def bin_an_array(input_array, bins, radius_arr=None):
    """Bin the points in an array by radius (n-dim)

    Parameters
    ----------
    input_array: np.ndarray
        array over which to bin
    bins: np.ndarray
        the bins
    radius_arr: np.ndarray
        optional array of |k| to avoid recalculation
    """
    if radius_arr is None:
        radius_arr = radius_array(input_array)

    radius_flat = radius_arr.flat
    arr_flat = input_array.flat

    counts_histo = np.histogram(radius_flat, bins)[0]
    binsum_histo = np.histogram(radius_flat, bins,
                                weights=arr_flat)[0]

    binavg = binsum_histo / counts_histo.astype(float)

    return counts_histo, binavg


def azimuthal_average(image, center=None, bw = 3):
    """
    Calculate the azimuthally averaged radial profile.

    Parameters
    ----------
    image : np.ndarray
        The 2D image.
    center : array_like, optional
        The [x,y] pixel coordinates used as the center. The default is
        None, which then uses the center of the image (including
        fractional pixels).

    Returns
    -------
    bl : np.ndarray
        The lower limit of the bin radius.
    radial_prof : np.ndarray
        The radial averaged profile.

    http://www.astrobetter.com/wiki/tiki-index.php?page=python_radial_profiles
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    maxr = np.array([image.shape[0] - center[0], image.shape[1] - center[1]]).min()

    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    maxind = np.searchsorted(r_sorted, maxr+0.00001)
    r_sorted = r_sorted[:maxind]
    i_sorted = i_sorted[:maxind]

    numbins = int(maxr / bw)
    maxr = numbins * bw

    bn = np.linspace(0.0, maxr, numbins+1)

    bc = 0.5*(bn[1:] + bn[:-1])
    bl = bn[:-1]

    # Get the integer part of the radii (bin size = 1)
    #r_int = r_sorted.astype(int)
    r_int = np.digitize(r_sorted, bn)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    #pdb.set_trace()
    rind = np.where(deltar)[0]       # location of changed radius
    rind = np.insert(rind, 0, -1)
    nr = rind[1:] - rind[:-1]        # number of radius bin

    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=image.dtype)
    csim = np.insert(csim, 0, 0.0)
    tbin = csim[rind[1:]+1] - csim[rind[:-1]+1]

    radial_prof = tbin / nr

    return bl, radial_prof


def azimuthal_average_nd(image, center=None, bw = 3):
    """
    Calculate the azimuthally averaged radial profile.

    Parameters
    ----------
    image : np.ndarray
        The 2D image.
    center : array_like, optional
        The [x,y] pixel coordinates used as the center. The default is
        None, which then uses the center of the image (including
        fractional pixels).

    Returns
    -------
    bl : np.ndarray
        The lower limit of the bin radius.
    radial_prof : np.ndarray
        The radial averaged profile.
    """
    # Calculate the indices from the image
    ia = np.indices(image.shape)
    ia = np.rollaxis(ia, 0, ia.ndim)

    if not center:
        center = (np.array(image.shape) - 1.0) / 2.0

    r = np.sum((ia - center)**2, axis=-1)**0.5

    # Get sorted radii
    #maxr = np.array([image.shape[0] - center[0], image.shape[1] - center[1]]).min()
    maxr = (np.array(image.shape) - center).min()

    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    maxind = np.searchsorted(r_sorted, maxr+0.00001)
    r_sorted = r_sorted[:maxind]
    i_sorted = i_sorted[:maxind]

    numbins = int(maxr / bw)
    maxr = numbins * bw

    bn = np.linspace(0.0, maxr, numbins+1)

    bc = 0.5*(bn[1:] + bn[:-1])
    bl = bn[:-1]

    # Get the integer part of the radii (bin size = 1)
    #r_int = r_sorted.astype(int)
    r_int = np.digitize(r_sorted, bn)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    #pdb.set_trace()
    rind = np.where(deltar)[0]       # location of changed radius
    rind = np.insert(rind, 0, -1)
    nr = rind[1:] - rind[:-1]        # number of radius bin

    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=image.dtype)
    csim = np.insert(csim, 0, 0.0)
    tbin = csim[rind[1:]+1] - csim[rind[:-1]+1]

    radial_prof = tbin / nr

    return bl, radial_prof
