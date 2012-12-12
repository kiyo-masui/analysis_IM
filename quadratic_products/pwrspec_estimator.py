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
import scipy.special
import math
from core import algebra
import copy
import gc
from utils import fftutil
from utils import binning
import h5py


def cross_power_est(arr1, arr2, weight1, weight2,
                    window="blackman", nonorm=False):
    """Calculate the cross-power spectrum of a two nD fields.

    The arrays must be identical and have the same length (physically
    and in pixel number) along each axis.

    inputs are clobbered to save memory
    """

    # make the axes
    ndim = arr1.ndim
    k_axes = tuple(["k_" + axis_name for axis_name in arr1.axes])
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

    if window:
        window_function = fftutil.window_nd(arr1.shape, name=window)
        weight1 *= window_function
        weight2 *= window_function
        del window_function

    arr1 *= weight1
    arr2 *= weight2

    # correct for the weighting
    fisher_diagonal = np.sum(weight1 * weight2)

    fft_arr1 = np.fft.fftshift(np.fft.fftn(arr1))

    fft_arr2 = np.fft.fftshift(np.fft.fftn(arr2))

    fft_arr1 *= fft_arr2.conj()
    xspec = fft_arr1.real
    del fft_arr1, fft_arr2
    gc.collect()

    xspec /= fisher_diagonal

    # make the axes
    xspec_arr = algebra.make_vect(xspec, axis_names=k_axes)

    xspec_arr.info = info
    #print xspec_arr.get_axis("k_dec")

    if not nonorm:
        xspec_arr *= width.prod()

    return xspec_arr


def cross_power_est_highmem(arr1, arr2, weight1, weight2,
                    window="blackman", nonorm=False):
    """Calculate the cross-power spectrum of a two nD fields.

    The arrays must be identical and have the same length (physically
    and in pixel number) along each axis.

    Same goal as above without the emphasis on saving memory.
    This is the "tried and true" legacy function.
    """
    if window:
        window_function = fftutil.window_nd(arr1.shape, name=window)
        weight1 *= window_function
        weight2 *= window_function

    warr1 = arr1 * weight1
    warr2 = arr2 * weight2
    ndim = arr1.ndim

    fft_arr1 = np.fft.fftshift(np.fft.fftn(warr1))
    fft_arr2 = np.fft.fftshift(np.fft.fftn(warr2))
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

    if not nonorm:
        xspec_arr *= width.prod()

    return xspec_arr


def make_unitless(xspec_arr, radius_arr=None, ndim=None):
    """multiply by  surface area in ndim sphere / 2pi^ndim * |k|^D
    (e.g. k^3 / 2 pi^2 in 3D)
    """
    if radius_arr is None:
        radius_arr = binning.radius_array(xspec_arr)

    if ndim is None:
        ndim = xspec_arr.ndim

    factor = 2. * math.pi ** (ndim / 2.) / scipy.special.gamma(ndim / 2.)
    factor /= (2. * math.pi) ** ndim
    return xspec_arr * radius_arr ** ndim * factor


def convert_2d_to_1d_pwrspec(pwr_2d, counts_2d, bin_kx, bin_ky, bin_1d,
                             weights_2d=None, nullval=np.nan,
                             null_zero_counts=True):
    """take a 2D power spectrum and the counts matrix (number of modex per k
    cell) and return the binned 1D power spectrum
    pwr_2d is the 2D power
    counts_2d is the counts matrix
    bin_kx is the x-axis
    bin_ky is the x-axis
    bin_1d is the k vector over which to return the result
    weights_2d is an optional weight matrix in 2d; otherwise use counts

    null_zero_counts sets elements of the weight where there are zero counts to
    zero.
    """
    # find |k| across the array
    index_array = np.indices(pwr_2d.shape)
    scale_array = np.zeros(index_array.shape)
    scale_array[0, ...] = bin_kx[index_array[0, ...]]
    scale_array[1, ...] = bin_ky[index_array[1, ...]]
    scale_array = np.rollaxis(scale_array, 0, scale_array.ndim)
    radius_array = np.sum(scale_array ** 2., axis=-1) ** 0.5

    radius_flat = radius_array.flatten()
    pwr_2d_flat = pwr_2d.flatten()
    counts_2d_flat = counts_2d.flatten()
    if weights_2d is not None:
        weights_2d_flat = weights_2d.flatten()
    else:
        weights_2d_flat = counts_2d_flat.astype(float)

    weight_pwr_prod = weights_2d_flat * pwr_2d_flat
    weight_pwr_prod[np.isnan(weight_pwr_prod)] = 0.
    weight_pwr_prod[np.isinf(weight_pwr_prod)] = 0.
    weights_2d_flat[np.isnan(weight_pwr_prod)] = 0.
    weights_2d_flat[np.isinf(weight_pwr_prod)] = 0.

    if null_zero_counts:
        weight_pwr_prod[counts_2d_flat == 0] = 0.
        weights_2d_flat[counts_2d_flat == 0] = 0.

    counts_histo = np.histogram(radius_flat, bin_1d,
                                weights=counts_2d_flat)[0]

    weights_histo = np.histogram(radius_flat, bin_1d,
                                weights=weights_2d_flat)[0]

    binsum_histo = np.histogram(radius_flat, bin_1d,
                                weights=weight_pwr_prod)[0]

    # explicitly handle cases where the counts are zero
    #binavg = np.zeros_like(binsum_histo)
    #binavg[weights_histo > 0.] = binsum_histo[weights_histo > 0.] / \
    #                             weights_histo[weights_histo > 0.]
    #binavg[weights_histo <= 0.] = nullval
    old_settings = np.seterr(invalid="ignore")
    binavg = binsum_histo / weights_histo
    binavg[np.isnan(binavg)] = nullval
    # note that if the weights are 1/sigma^2, then the variance of the weighted
    # sum is just 1/sum of weights; so return Gaussian errors based on that
    gaussian_errors = np.sqrt(1./weights_histo)
    gaussian_errors[np.isnan(binavg)] = nullval

    np.seterr(**old_settings)

    return counts_histo, gaussian_errors, binavg


def calculate_xspec(cube1, cube2, weight1, weight2,
                    window="blackman", unitless=True, bins=None,
                    truncate=False, nbins=40, logbins=True, return_3d=False):

    print "finding the signal power spectrum"
    pwrspec3d_signal = cross_power_est(cube1, cube2, weight1, weight2,
                                       window=window)

    radius_arr = binning.radius_array(pwrspec3d_signal)

    if bins is None:
        bins = binning.suggest_bins(pwrspec3d_signal,
                                          truncate=truncate,
                                          logbins=logbins,
                                          nbins=nbins,
                                          radius_arr=radius_arr)

    if unitless:
        print "making the power spectrum unitless"
        pwrspec3d_signal = make_unitless(pwrspec3d_signal,
                                         radius_arr=radius_arr)

    print "calculating the 1D histogram"
    counts_histo, binavg = binning.bin_an_array(pwrspec3d_signal, bins,
                                                radius_arr=radius_arr)

    del radius_arr
    gc.collect()

    print "calculating the 2D histogram"
    # find the k_perp by not including k_nu in the distance
    radius_arr_perp = binning.radius_array(pwrspec3d_signal,
                                           zero_axes=[0])

    # find the k_perp by not including k_RA,Dec in the distance
    radius_arr_parallel = binning.radius_array(pwrspec3d_signal,
                                               zero_axes=[1, 2])

    # TODO: do better independent binning; for now:
    bins_x = copy.deepcopy(bins)
    bins_y = copy.deepcopy(bins)
    counts_histo_2d, binavg_2d = binning.bin_an_array_2d(pwrspec3d_signal,
                                                         radius_arr_perp,
                                                         radius_arr_parallel,
                                                         bins_x, bins_y)

    del radius_arr_perp
    del radius_arr_parallel
    gc.collect()

    bin_left_x, bin_center_x, bin_right_x = binning.bin_edges(bins_x,
                                                              log=logbins)

    bin_left_y, bin_center_y, bin_right_y = binning.bin_edges(bins_y,
                                                              log=logbins)

    bin_left, bin_center, bin_right = binning.bin_edges(bins, log=logbins)

    pwrspec2d_product = {}
    pwrspec2d_product['bin_x_left'] = bin_left_x
    pwrspec2d_product['bin_x_center'] = bin_center_x
    pwrspec2d_product['bin_x_right'] = bin_right_x
    pwrspec2d_product['bin_y_left'] = bin_left_y
    pwrspec2d_product['bin_y_center'] = bin_center_y
    pwrspec2d_product['bin_y_right'] = bin_right_y
    pwrspec2d_product['counts_histo'] = counts_histo_2d
    pwrspec2d_product['binavg'] = binavg_2d

    pwrspec1d_product = {}
    pwrspec1d_product['bin_left'] = bin_left
    pwrspec1d_product['bin_center'] = bin_center
    pwrspec1d_product['bin_right'] = bin_right
    pwrspec1d_product['counts_histo'] = counts_histo
    pwrspec1d_product['binavg'] = binavg

    if not return_3d:
        return pwrspec2d_product, pwrspec1d_product
    else:
        return pwrspec3d_signal, pwrspec2d_product, pwrspec1d_product


def calculate_xspec_file(cube1_file, cube2_file, bins,
                    weight1_file=None, weight2_file=None,
                    truncate=False, window="blackman",
                    return_3d=False, unitless=True):

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

    print cube1.shape, cube2.shape, weight1.shape, weight2.shape
    return calculate_xspec(cube1, cube2, weight1, weight2, bins=bins,
                           window=window, unitless=unitless,
                           truncate=truncate, return_3d=return_3d)


def load_transferfunc(beamtransfer_file, modetransfer_file, treatment_list,
                      beam_treatment="0modes"):
    r"""Given two hd5 files containing a beam and mode loss transfer function,
    load them and make a composite transfer function"""
    # TODO: check that both files have the same treatment cases

    if modetransfer_file is not None:
        print "Applying 2d transfer from " + modetransfer_file
        modetransfer_2d = h5py.File(modetransfer_file, "r")

    if beamtransfer_file is not None:
        print "Applying 2d transfer from " + beamtransfer_file
        beamtransfer_2d = h5py.File(beamtransfer_file, "r")

    # make a list of treatments or cross-check if given one
    if modetransfer_file is not None:
        if treatment_list is None:
            treatment_list = modetransfer_2d.keys()
        else:
            assert modetransfer_2d.keys() == treatment_list, \
                    "mode transfer treatments do not match data"
    else:
        if treatment_list is None:
            treatment_list = [beam_treatment]

    # given both
    if (beamtransfer_file is not None) and (modetransfer_file is not None):
        print "using the product of beam and mode transfer functions"
        transfer_dict = {}

        for treatment in modetransfer_2d:
            transfer_dict[treatment] = modetransfer_2d[treatment].value
            transfer_dict[treatment] *= \
                                    beamtransfer_2d[beam_treatment].value

    # given mode only
    if (beamtransfer_file is None) and (modetransfer_file is not None):
        print "using just the mode transfer function"
        transfer_dict = {}

        for treatment in modetransfer_2d:
            transfer_dict[treatment] = modetransfer_2d[treatment].value

    # given beam only
    if (beamtransfer_file is not None) and (modetransfer_file is None):
        print "using just the beam transfer function"
        transfer_dict = {}
        for treatment in treatment_list:
            transfer_dict[treatment] = beamtransfer_2d[beam_treatment].value

    # no transfer function
    if (beamtransfer_file is None) and (modetransfer_file is None):
        print "not using transfer function"
        transfer_dict = None

    return transfer_dict


def test_with_random(unitless=True):
    """Test the power spectral estimator using a random noise cube"""

    delta = 1.33333
    cube1 = algebra.make_vect(np.random.normal(0, 1,
                                size=(257, 124, 68)))

    info = {'axes': ["freq", "ra", "dec"], 'type': 'vect',
            'freq_delta': delta / 3.78, 'freq_centre': 0.,
            'ra_delta': delta / 1.63, 'ra_centre': 0.,
            'dec_delta': delta, 'dec_centre': 0.}
    cube1.info = info
    cube2 = copy.deepcopy(cube1)

    weight1 = algebra.ones_like(cube1)
    weight2 = algebra.ones_like(cube2)

    bin_left, bin_center, bin_right, counts_histo, binavg = \
                    calculate_xspec(cube1, cube2, weight1, weight2,
                                    window="blackman",
                                    truncate=False,
                                    nbins=40,
                                    unitless=unitless,
                                    logbins=True)

    if unitless:
        pwrspec_input = bin_center ** 3. / 2. / math.pi / math.pi
    else:
        pwrspec_input = np.ones_like(bin_center)

    volume = 1.
    for axis_name in cube1.axes:
        axis_vector = cube1.get_axis(axis_name)
        volume *= abs(axis_vector[1] - axis_vector[0])

    pwrspec_input *= volume

    for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, binavg,
                        pwrspec_input):
        print ("%10.15g " * 6) % specdata


if __name__ == '__main__':
    test_with_random(unitless=False)
