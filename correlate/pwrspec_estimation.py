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
        print factor

    print width
    xspec_arr *= width.prod()

    return xspec_arr


def binpwrspec(input_array, truncate=True, nbins=40):
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
    bins = np.logspace(math.log10(radius_sorted[1]),
                       math.log10(max_r),
                       num=(nbins+1), endpoint=True)

    bin_left, bins_right = bins[:-1], bins[1:]
    # TODO: this is not correct in log-bins, but does it matter to anyone?
    bin_center = 0.5*(bin_left + bin_right)

    # calculate the binning useing a histogram
    counts_histo, bin_edges = np.histogram(radius_sorted, bins)
    binsum_histo, bin_edges = np.histogram(radius_sorted, bins,
                                           weights=arr_sorted)

    binavg = binsum_histo / counts_histo.astype(float)

    for bl, bc, br, ct, pk in zip(bin_left, bin_center,
                                  bin_right, counts, binavg):
        print bl, bc, br, ct, pk

def binpwrspec_alt(input_array, truncate=True, nbins=40):
    """
    The alternate implementation with digitize, cumsum, etc. has problems for
    bins with zero points, other reliability issues.
    """
    radius = radius_array(input_array)

    sort_ind = np.argsort(radius.flat)
    radius_sorted = radius.flat[sort_ind]
    arr_sorted = input_array.flat[sort_ind]

    axis_range = np.zeros((input_array.ndim))
    for axis_index in range(input_array.ndim):
        axis_name = input_array.axes[axis_index]
        axis = input_array.get_axis(axis_name)
        axis_range[axis_index] = axis.max()

    if truncate:
        max_r = axis_range.min()
        maxind = np.searchsorted(radius_sorted, max_r, side='right')
        radius_sorted = radius_sorted[:maxind]
        arr_sorted = arr_sorted[:maxind]

    # TODO: should the enpoint be here?
    # ignore the bin at k=0
    bins = np.logspace(math.log10(radius_sorted[1]),
                       math.log10(radius_sorted[-1]),
                       num=(nbins+1), endpoint=True)

    bin_left, bins_right = bins[:-1], bins[1:]
    # TODO: this is not correct in log-bins, but does it matter to anyone?
    bin_center = 0.5*(bin_left + bin_right)

    # calculate the binning another way
    bin_radii_ind = np.digitize(radius_sorted, bins)
    delta_radius = bin_radii_ind[1:] - bin_radii_ind[:-1]
    boundaries = np.where(delta_radius)[0]
    #boundaries = np.insert(boundaries, 0, -1)
    counts = boundaries[1:] - boundaries[:-1]

    cumulative_signal = np.cumsum(arr_sorted, dtype=input_array.dtype)
    cumulative_signal = np.insert(cumulative_signal, 0, 0.)
    binsum = cumulative_signal[boundaries[1:]+1] - \
             cumulative_signal[boundaries[:-1]+1]

    binavg = binsum / counts.astype(float)

    for bl, bc, br, ct, pk in zip(bin_left, bin_center,
                                  bin_right, counts, binavg):
        print bl, bc, br, ct, pk


if __name__ == '__main__':
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
    pwrspec3d_signal = crossps(cube1, cube2, weight1, weight2, window=True)
    print "binning the power spectrum"
    binpwrspec(pwrspec3d_signal, truncate=True)
    print "fetching the P(k) of the input sims"
    k_vec = np.logspace(math.log10(1.e-2),
                        math.log10(5.),
                        num=55, endpoint=True)

    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(k_vec)
    pwrspec_input *= k_vec**3./2./math.pi/math.pi
    #for k, p in zip(k_vec, pwrspec_input):
    #    print "%10.15g %10.15g" % (k, p)

