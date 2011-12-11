"""
    TODO:
        test that the array, weight shapes and axes are compatible
        FIX THE VOLUME NORMALIZATION!
        calculate the full window function instead of just the diagonal
            i, j -> delta k -> nearest index k
"""
import numpy as np
import math
from core import algebra
from utils import radialprofile
from utils import data_paths
from simulations import ps_estimation
from simulations import corr21cm


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

    return np.sum(scale_array**2., axis=-1)**0.5


def crossps(arr1, arr2, weight1, weight2, window=True, unitless=True):
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
    unitless: boolean
        Multiply by e.g. k^3 / 2 pi^2 in 3D
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
        delta_axis = abs(axis_vector[1]-axis_vector[0])
        width[axis_index] = delta_axis

        k_axis = np.fft.fftshift(np.fft.fftfreq(n_axis, d=delta_axis))
        k_axis *= 2. * math.pi
        delta_k_axis = abs(k_axis[1]-k_axis[0])

        k_name = k_axes[axis_index]
        info[k_name + "_delta"] = delta_k_axis
        info[k_name + "_centre"] = 0.
        #print k_axis
        #print k_name, n_axis, delta_axis

    xspec_arr.info = info
    #print xspec_arr.get_axis("k_dec")

    if unitless:
        radius = radius_array(xspec_arr)
        # multiply by surface area in ndim sphere / 2pi^ndim
        factor = 2. * math.pi ** (ndim / 2.) / math.gamma(ndim / 2.)
        factor /= (2. * math.pi) ** ndim
        xspec_arr *= radius ** ndim * factor

    xspec_arr *= width.prod()

    return xspec_arr


def binpwrspec(input_array, truncate=True, nbins=40, logbins=True):
    """Bin the points in an array by radius
    """
    radius = radius_array(input_array)
    sort_ind = np.argsort(radius.flat)
    radius_sorted = radius.flat[sort_ind]
    arr_sorted = input_array.flat[sort_ind]

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
                           num=(nbins+1), endpoint=True)

        bin_left, bin_right = bins[:-1], bins[1:]
        bin_center = 10**(0.5*(np.log10(bin_left) +
                               np.log10(bin_right)))
    else:
        bins = np.linspace(radius_sorted[1], max_r,
                           num=(nbins+1), endpoint=True)

        bin_left, bin_right = bins[:-1], bins[1:]
        bin_center = 0.5*(np.log10(bin_left) +
                          np.log10(bin_right))

    counts_histo, bin_edges = np.histogram(radius_sorted, bins)
    binsum_histo, bin_edges = np.histogram(radius_sorted, bins,
                                           weights=arr_sorted)

    binavg = binsum_histo / counts_histo.astype(float)

    return bin_left, bin_center, bin_right, counts_histo, binavg


def test_with_simulation(unitless=True):
    datapath_db = data_paths.DataPath()
    filename = datapath_db.fetch('sim_15hr_physical', intend_read=True,
                                 pick='0')
    zfilename = datapath_db.fetch('sim_15hr', intend_read=True,
                                 pick='0')
    print filename
    cube1 = algebra.make_vect(algebra.load(filename))
    cube2 = algebra.make_vect(algebra.load(filename))
    zspace_cube = algebra.make_vect(algebra.load(zfilename))

    weight1 = algebra.ones_like(cube1)
    weight2 = algebra.ones_like(cube2)

    print "finding the signal power spectrum"
    pwrspec3d_signal = crossps(cube1, cube2, weight1, weight2,
                               window=True, unitless=unitless)
    print "binning the power spectrum"
    bin_left, bin_center, bin_right, counts_histo, binavg = \
                        binpwrspec(pwrspec3d_signal, truncate=True)

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
        pwrspec_input *= bin_center**3./2./math.pi/math.pi

    for bl, bc, br, ct, pk, th in zip(bin_left, bin_center,
                                  bin_right, counts_histo, binavg,
                                  pwrspec_input):
        print bl, bc, br, ct, pk, th


def test_with_random(unitless=True):
    import copy

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


    print "finding the signal power spectrum"
    pwrspec3d_signal = crossps(cube1, cube2, weight1, weight2,
                               window=True, unitless=unitless)
    print "binning the power spectrum"
    bin_left, bin_center, bin_right, counts_histo, binavg = \
                        binpwrspec(pwrspec3d_signal, truncate=True)

    if unitless:
        pwrspec_input = bin_center**3./2./math.pi/math.pi*delta**3.
    else:
        pwrspec_input = np.ones_like(bin_center)*delta**3.

    for bl, bc, br, ct, pk, th in zip(bin_left, bin_center,
                                  bin_right, counts_histo, binavg,
                                  pwrspec_input):
        print bl, bc, br, ct, pk, th


if __name__ == '__main__':
    test_with_simulation(unitless=False)
