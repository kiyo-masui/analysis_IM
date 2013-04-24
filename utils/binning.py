# TODO: consider moving azimuthalAverage and _nd version to the attic
# TODO: rerun unit test on histogram3d, develop unit tests for all funcs
import numpy as np
import math
from core import algebra
import unittest


def radius_array(input_array, zero_axes=[]):
    """Find the Euclidian distance of all the points in an array using their
    axis meta-data. e.g. x_axis[0]^2 + y_axis[0]^2. (n-dim)
    optionally give the list of axes to set to zero -- for example, in 3D:
    zero_axes = [0] leaves the x^2 out of x^2+y^2+z^2
    zero_axes = [1,2] leave the y^2 and z^2 out of x^2+y^2+z^2 (e.g. x^2)
    """
    index_array = np.indices(input_array.shape)
    scale_array = np.zeros(index_array.shape)

    for axis_index in range(input_array.ndim):
        axis_name = input_array.axes[axis_index]
        axis = input_array.get_axis(axis_name)
        if axis_index in zero_axes:
            axis = np.zeros_like(axis)

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


def bin_an_array_2d(input_array, radius_array_x, radius_array_y,
                    bins_x, bins_y):
    """Bin the points in an array by radius (n-dim)

    Parameters
    ----------
    input_array: np.ndarray
        array over which to bin
    radius_array_x: np.ndarray
        array of the x-vectors at each point in the series
    radius_array_y: np.ndarray
        array of the y-vectors at each point in the series
    bins_x and bins_y: np.ndarray
        the bins in x and y
    """
    radius_flat_x = radius_array_x.flatten()
    radius_flat_y = radius_array_y.flatten()
    arr_flat = input_array.flatten()

    counts_histo = np.histogram2d(radius_flat_x, radius_flat_y,
                                  bins=(bins_x, bins_y))[0]
    binsum_histo = np.histogram2d(radius_flat_x, radius_flat_y,
                                  bins=(bins_x, bins_y),
                                  weights=arr_flat)[0]

    binavg = binsum_histo / counts_histo.astype(float)

    return counts_histo, binavg


def find_edges(axis, delta=None):
    """
    service function for bin_catalog_data which
    finds the bin edges for the histogram
    """
    if delta is None:
        delta = axis[1] - axis[0]

    edges = np.array(axis) - delta / 2.
    return np.append(edges, edges[-1] + delta)


def print_edges(sample, edges, name):
    """print bin edges for a catalog
    """
    print "Binning %s from range (%5.3g, %5.3g) into (%5.3g, %5.3g)" % (
           name, min(sample), max(sample), min(edges), max(edges))


def histogram3d(sample, xedges, yedges, zedges):
    """Make a 3D histogram from the sample and edge specification
    indices in the sample: 0=x, 1=y, 2=z;
    histogramdd was having problems with the galaxy catalogs
    """
    numcatalog = sample.size
    x_size = xedges.size - 1
    y_size = yedges.size - 1
    z_size = zedges.size - 1
    box_index = np.zeros(numcatalog)
    count_array = np.zeros((x_size + 1) * (y_size + 1) * (z_size + 1))
    # the final array to return is the value within the bin
    count_cube = np.zeros((x_size, y_size, z_size))

    # find which bin each galaxies lies in
    x_index = np.digitize(sample[:, 0], xedges)
    y_index = np.digitize(sample[:, 1], yedges)
    z_index = np.digitize(sample[:, 2], zedges)

    # digitize puts values outside of the bins either to 0 or len(bins)
    x_out = np.logical_or((x_index == 0), (x_index == (x_size + 1)))
    y_out = np.logical_or((y_index == 0), (y_index == (y_size + 1)))
    z_out = np.logical_or((z_index == 0), (z_index == (z_size + 1)))
    # now flag all those point which are inside the region
    box_in = np.logical_not(np.logical_or(np.logical_or(x_out, y_out), z_out))

    # the 0th bin center is recorded in the digitized index=1, so shift
    # also throw out points that are not in the volume
    x_index = x_index[box_in] - 1
    y_index = y_index[box_in] - 1
    z_index = z_index[box_in] - 1

    box_index = x_index + y_index * x_size + z_index * x_size * y_size

    # note that bincount will only count up to the largest object in the list,
    # which may be smaller than the dimension of the full count cube
    try:
        count_array[0:max(box_index) + 1] = np.bincount(box_index)

        # make the x y and z axes which index the bincount output
        count_index = np.arange(x_size * y_size * z_size)
        zind = count_index / (x_size * y_size)
        yind = (count_index - x_size * y_size * zind) / x_size
        xind = count_index - x_size * y_size * zind - x_size * yind

        #count_cube[xind, yind, zind] = count_array[xind + yind * x_size +
        #                                           zind * x_size * y_size]
        count_cube[xind, yind, zind] = count_array[count_index]
        #split_indices = cartesian((np.arange(z_size),
        #                           np.arange(y_size),
        #                           np.arange(x_size)))
        #count_cube[split_indices] = count_array[count_index]
    except MemoryError:
        print "histogram3d: all points out of the volume"

    return count_cube


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


class CatalogGriddingTest(unittest.TestCase):
    """Unit test class for catalog gridding
    """

    def test_simple(self):
        """bin a simple 3x3x3 array
        """

        parent_axis = np.array([0.25, 0.75, 1.25])
        edges = find_edges(parent_axis)
        self.assertTrue(np.array_equal(edges, [0., 0.5, 1., 1.5]))

        # test a sample (with some outliers)
        sample = np.array([[0., 0., 0.],
                           [0.75, 0., 0.],
                           [1.25, 0., 0.],
                           [1.75, 0., 0.],
                           [0., 0., 0.],
                           [0.75, 0.75, 0.75],
                           [1.25, 1.25, 1.25]])

        result = histogram3d(sample, edges, edges, edges)
        alternate, histo_edges = np.histogramdd(sample,
                                                bins=[edges, edges, edges])

        answer = np.array([[[2,  0,  0],
                            [0,  0,  0],
                            [0,  0,  0]],
                           [[1,  0,  0],
                            [0,  1,  0],
                            [0,  0,  0]],
                           [[1,  0,  0],
                            [0,  0,  0],
                            [0,  0,  1]]])

        self.assertTrue(np.array_equal(answer, result))
        self.assertTrue(np.array_equal(alternate, result))

        # test the case where no points are in the volume
        sample2 = np.array([[-1., -1., -1.]])
        result2 = histogram3d(sample2, edges, edges, edges)
        alternate2, histo_edges = np.histogramdd(sample2,
                                                bins=[edges, edges, edges])

        answer2 = np.zeros((3, 3, 3), dtype=int)
        self.assertTrue(np.array_equal(answer2, result2))
        self.assertTrue(np.array_equal(alternate2, result2))

    def test_timing(self):
        """compare the timing of histogram3d and histogramdd"""
        # TODO: compare sum of two histogram methods;
        # edge cases do not seem to match
        # TODO: speed up histogram3d class
        edges = np.array([0., 0.25, 0.75, 1.])
        sample = np.random.rand(1e7, 3)

        # profiling tools do not seem to work well with numpy
        start = time.clock()
        result = histogram3d(sample, edges, edges, edges)
        end = time.clock()
        print (end - start) / 1000.
        alternate, histo_edges = np.histogramdd(sample,
                                                bins=[edges, edges, edges])
        endalt = time.clock()
        print (endalt - end) / 1000.
        print result - alternate
