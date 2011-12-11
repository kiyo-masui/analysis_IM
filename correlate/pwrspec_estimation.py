"""
    TODO:
        test that the array, weight shapes and axes are compatible
        calculate the full window function instead of just the diagonal
            i, j -> delta k -> nearest index k
        make plots of the real data's 3D power
"""
import numpy as np
import math
from core import algebra
from utils import data_paths
from simulations import ps_estimation
from simulations import corr21cm
import copy


def radius_array(input_array):
    """Find the Euclidian distance of all the points in an array using their
    axis meta-data. e.g. x_axis[0]^2 + y_axis[0]^2
    """
    index_array = np.indices(input_array.shape)
    scale_array = np.zeros(index_array.shape)

    for axis_index in range(input_array.ndim):
        axis_name = input_array.axes[axis_index]
        axis = input_array.get_axis(axis_name)
        scale_array[axis_index, ...] = axis[index_array[axis_index, ...]]

    scale_array = np.rollaxis(scale_array, 0, scale_array.ndim)

    return np.sum(scale_array ** 2., axis=-1) ** 0.5


def crossps(arr1, arr2, weight1, weight2, window=True):
    """Calculate the radially average cross-power spectrum of a two nD fields.

    The arrays must be identical and have the same length (physically
    and in pixel number) along each axis.

    Parameters
    ----------
    arr1, arr2: np.ndarray
        The cubes to calculate the cross-power spectrum of. If arr2 is
        None then return the standard (auto-)power spectrum.
    window: boolean
        Apply an additional Blackman window
    """
    if window:
        blackman_window = ps_estimation.blackman_nd(arr1.shape)
        weight1 *= blackman_window
        weight2 *= blackman_window

    arr1 *= weight1
    arr2 *= weight2
    ndim = arr1.ndim

    fft_arr1 = np.fft.fftshift(np.fft.fftn(arr1))
    fft_arr2 = np.fft.fftshift(np.fft.fftn(arr2))
    xspec = fft_arr1 * fft_arr2.conj()
    xspec = xspec.real

    # correct for the weighting
    product_weight = weight1 * weight2
    xspec /= np.sum(product_weight)

    # make the axes
    k_axes = tuple(["k_" + axis_name for axis_name in arr1.axes])
    xspec_arr = algebra.make_vect(xspec, axis_names=k_axes)

    info = {'axes': k_axes, 'type': 'vect'}
    width = np.zeros(ndim)
    for axis_index in range(ndim):
        n_axis = arr1.shape[axis_index]
        axis_name = arr1.axes[axis_index]
        axis_vector = arr1.get_axis(axis_name)
        delta_axis = abs(axis_vector[1] - axis_vector[0])
        width[axis_index] = delta_axis

        k_axis = np.fft.fftshift(np.fft.fftfreq(n_axis, d=delta_axis))
        k_axis *= 2. * math.pi
        delta_k_axis = abs(k_axis[1] - k_axis[0])

        k_name = k_axes[axis_index]
        info[k_name + "_delta"] = delta_k_axis
        info[k_name + "_centre"] = 0.
        #print k_axis
        #print k_name, n_axis, delta_axis

    xspec_arr.info = info
    #print xspec_arr.get_axis("k_dec")

    xspec_arr *= width.prod()

    return xspec_arr


def make_unitless(xspec_arr, radius_arr=None):
    """multiply by  surface area in ndim sphere / 2pi^ndim * |k|^D
    (e.g. k^3 / 2 pi^2 in 3D)
    """
    if radius_arr is None:
        radius_arr = radius_array(xspec_arr)

    ndim = xspec_arr.ndim
    factor = 2. * math.pi ** (ndim / 2.) / math.gamma(ndim / 2.)
    factor /= (2. * math.pi) ** ndim
    return xspec_arr * radius_arr ** ndim * factor


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


def binpwrspec(input_array, bins, radius_arr=None):
    """Bin the points in an array by radius

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


def calculate_xspec(cube1, cube2, weight1, weight2,
                    window=True, truncate=False, nbins=40,
                    unitless=True, logbins=True):

    print "finding the signal power spectrum"
    pwrspec3d_signal = crossps(cube1, cube2, weight1, weight2,
                               window=window)
    radius_arr = radius_array(pwrspec3d_signal)

    if unitless:
        print "making the power spectrum unitless"
        pwrspec3d_signal = make_unitless(pwrspec3d_signal,
                                         radius_arr=radius_arr)

    print "binning the power spectrum"
    bins = suggest_bins(pwrspec3d_signal, truncate=truncate, logbins=logbins,
                        nbins=nbins, radius_arr=radius_arr)

    counts_histo, binavg = binpwrspec(pwrspec3d_signal, bins,
                                      radius_arr=radius_arr)

    bin_left, bin_center, bin_right = bin_edges(bins, log=logbins)

    return bin_left, bin_center, bin_right, counts_histo, binavg


def test_with_simulation(unitless=True):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()
    filename = datapath_db.fetch('simideal_15hr_physical', intend_read=True,
                                 pick='1')
    zfilename = datapath_db.fetch('simideal_15hr', intend_read=True,
                                 pick='1')
    print filename
    cube1 = algebra.make_vect(algebra.load(filename))
    cube2 = algebra.make_vect(algebra.load(filename))
    zspace_cube = algebra.make_vect(algebra.load(zfilename))

    weight1 = algebra.ones_like(cube1)
    weight2 = algebra.ones_like(cube2)

    bin_left, bin_center, bin_right, counts_histo, binavg = \
                    calculate_xspec(cube1, cube2, weight1, weight2,
                                    window=True, truncate=False, nbins=40,
                                    unitless=unitless, logbins=True)

    #volume = 1.
    #for axis_name in cube1.axes:
    #    axis_vector = cube1.get_axis(axis_name)
    #    volume *= abs(axis_vector[1]-axis_vector[0])

    #k_vec = np.logspace(math.log10(1.e-2),
    #                    math.log10(5.),
    #                    num=55, endpoint=True)
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, binavg,
                        pwrspec_input):
        print ("%10.15g " * 6) % specdata


def test_with_random(unitless=True):
    """Test the power spectral estimator using a random noise cube"""

    nside = 256
    delta = 0.23

    cube1 = algebra.make_vect(np.random.normal(0, 1,
                                size=(nside, nside, nside)))

    info = {'axes': ["freq", "ra", "dec"], 'type': 'vect',
            'freq_delta': delta, 'freq_centre': 0.,
            'ra_delta': delta, 'ra_centre': 0.,
            'dec_delta': delta, 'dec_centre': 0.}
    cube1.info = info
    cube2 = copy.deepcopy(cube1)

    weight1 = algebra.ones_like(cube1)
    weight2 = algebra.ones_like(cube2)

    bin_left, bin_center, bin_right, counts_histo, binavg = \
                    calculate_xspec(cube1, cube2, weight1, weight2,
                                    window=True, truncate=False, nbins=40,
                                    unitless=unitless, logbins=True)

    if unitless:
        pwrspec_input = bin_center ** 3. / 2. / math.pi / math.pi
    else:
        pwrspec_input = np.ones_like(bin_center)

    pwrspec_input *= delta ** 3.

    for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, binavg,
                        pwrspec_input):
        print ("%10.15g " * 6) % specdata


if __name__ == '__main__':
    #test_with_simulation(unitless=True)
    test_with_random(unitless=True)
