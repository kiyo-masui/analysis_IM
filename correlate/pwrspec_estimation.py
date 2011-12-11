"""
    TODO:
        test that the array, weight shapes and axes are compatible
        calculate the full window function instead of just the diagonal
            i, j -> delta k -> nearest index k
        make plots of the real data's 3D power
        make uniform noise unit test
        make radius array unit test (anisotropic axes)
"""
import numpy as np
import math
from core import algebra
from utils import data_paths
from simulations import ps_estimation
from simulations import corr21cm
import multiprocessing
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
                    window=True, unitless=True, bins=None,
                    truncate=False, nbins=40, logbins=True):

    print "finding the signal power spectrum"
    pwrspec3d_signal = crossps(cube1, cube2, weight1, weight2,
                               window=window)
    radius_arr = radius_array(pwrspec3d_signal)

    if unitless:
        print "making the power spectrum unitless"
        pwrspec3d_signal = make_unitless(pwrspec3d_signal,
                                         radius_arr=radius_arr)

    print "binning the power spectrum"
    if bins is None:
        bins = suggest_bins(pwrspec3d_signal, truncate=truncate, logbins=logbins,
                            nbins=nbins, radius_arr=radius_arr)

    counts_histo, binavg = binpwrspec(pwrspec3d_signal, bins,
                                      radius_arr=radius_arr)

    bin_left, bin_center, bin_right = bin_edges(bins, log=logbins)

    return bin_left, bin_center, bin_right, counts_histo, binavg


def calculate_xspec_file(cube1_file, cube2_file, bins,
                    weight1_file=None, weight2_file=None,
                    window=True, unitless=True):
    """TODO: merge this with the code above"""

    cube1 = algebra.make_vect(algebra.load(cube1_file))
    cube2 = algebra.make_vect(algebra.load(cube2_file))

    if weight1_file is None:
        weight1 = algebra.ones_like(cube1)
    else:
        weight1 = algebra.make_vect(algebra.load(weight1_file))

    if weight2_file is None:
        weight2 = algebra.ones_like(cube2)
    else:
        weight2 = algebra.make_vect(algebra.load(weight2_file))

    return calculate_xspec(cube1, cube2, weight1, weight2, bins=bins,
                           window=window, unitless=unitless)


def wrap_xspec(param):
    """helper for multiprocessing.map; should toast"""
    (cube1_file, cube2_file, \
     weight1_file, weight2_file, \
     bins, window, unitless) = param

    return calculate_xspec_file(cube1_file, cube2_file, bins,
                                weight1_file=weight1_file,
                                weight2_file=weight2_file,
                                window=window, unitless=unitless)


def test_with_agg_simulation(unitless=True, parallel=True):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()
    filename = datapath_db.fetch('simideal_15hr_physical', intend_read=True)
    zfilename = datapath_db.fetch('simideal_15hr', intend_read=True,
                                 pick='1')

    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    # give no weights; just use Blackman window
    runlist = [(filename[1][index], filename[1][index], None, None, bins,
                True, unitless) for index in filename[0]]

    if parallel:
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-4)
        results = pool.map(wrap_xspec, runlist)
    else:
        for runitem in runlist:
            wrap_xspec(runitem)

    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    bin_left, bin_center, bin_right = bin_edges(bins, log=True)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    agg_array = np.zeros((bin_center.size, len(filename[0])))
    counter = 0
    for spec_output in results:
        (bin_left, bin_center, bin_right, counts_histo, binavg) = spec_output
        agg_array[:,counter] = binavg
        counter += 1
        #for specdata in zip(bin_left, bin_center,
        #                    bin_right, counts_histo, binavg,
        #                    pwrspec_input):
        #    print ("%10.15g " * 6) % specdata

    meanbin = np.mean(agg_array, axis=1)
    stdbin = np.std(agg_array, axis=1)

    for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, meanbin, stdbin,
                        pwrspec_input):
        print ("%10.15g " * 7) % specdata


def test_with_simulation(unitless=True):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()
    filename = datapath_db.fetch('simideal_15hr_physical', intend_read=True,
                                 pick='1')
    zfilename = datapath_db.fetch('simideal_15hr', intend_read=True,
                                 pick='1')
    print filename

    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    bin_left, bin_center, bin_right, counts_histo, binavg = \
                    calculate_xspec_file(filename, filename, bins,
                                    window=True, unitless=unitless)

    #volume = 1.
    #for axis_name in cube1.axes:
    #    axis_vector = cube1.get_axis(axis_name)
    #    volume *= abs(axis_vector[1]-axis_vector[0])

    #k_vec = np.logspace(math.log10(1.e-2),
    #                    math.log10(5.),
    #                    num=55, endpoint=True)
    zspace_cube = algebra.make_vect(algebra.load(zfilename))
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
    #test_with_agg_simulation(unitless=True)
    test_with_simulation(unitless=True)
    #test_with_random(unitless=True)
