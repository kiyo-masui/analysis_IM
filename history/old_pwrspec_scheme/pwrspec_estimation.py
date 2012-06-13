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
from utils import fftutil
from utils import binning


def cross_power_est(arr1, arr2, weight1, weight2, window="blackman"):
    """Calculate the radially average cross-power spectrum of a two nD fields.

    The arrays must be identical and have the same length (physically
    and in pixel number) along each axis.

    Parameters
    ----------
    arr1, arr2: np.ndarray
        The cubes to calculate the cross-power spectrum of. If arr2 is
        None then return the standard (auto-)power spectrum.
    window: boolean
        Apply an additional named window
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


def calculate_xspec(cube1, cube2, weight1, weight2,
                    window="blackman", unitless=True, bins=None,
                    truncate=False, nbins=40, logbins=True, return_3d=False):

    print "finding the signal power spectrum"
    pwrspec3d_signal = cross_power_est(cube1, cube2, weight1, weight2,
                               window=window)

    radius_arr = binning.radius_array(pwrspec3d_signal)

    # find the k_perp by not including k_nu in the distance
    radius_arr_perp = binning.radius_array(pwrspec3d_signal,
                                           zero_axes=[0])

    # find the k_perp by not including k_RA,Dec in the distance
    radius_arr_parallel = binning.radius_array(pwrspec3d_signal,
                                               zero_axes=[1, 2])

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

    print "calculating the 2D histogram"
    # TODO: do better independent binning; for now:
    bins_x = copy.deepcopy(bins)
    bins_y = copy.deepcopy(bins)
    counts_histo_2d, binavg_2d = binning.bin_an_array_2d(pwrspec3d_signal,
                                                         radius_arr_perp,
                                                         radius_arr_parallel,
                                                         bins_x, bins_y)

    print "calculating the 1D histogram"
    counts_histo, binavg = binning.bin_an_array(pwrspec3d_signal, bins,
                                                radius_arr=radius_arr)

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


def agg_stat_1d_pwrspec(pwr_1d, apply_1d_transfer=None):
    pwrshp_1d = pwr_1d[0]['binavg'].shape
    n_runs = len(pwr_1d)

    pwrmat_1d = np.zeros((n_runs, pwrshp_1d[0]))
    for index in range(n_runs):
        if apply_1d_transfer is not None:
            print apply_1d_transfer
            pwrmat_1d[index, :] = pwr_1d[index]['binavg'] / apply_1d_transfer
        else:
            pwrmat_1d[index, :] = pwr_1d[index]['binavg']

    mean_1d = np.mean(pwrmat_1d, axis=0)
    std_1d = np.std(pwrmat_1d, axis=0, ddof=1)
    corrmat_1d = np.corrcoef(np.transpose(pwrmat_1d))
    covmat_1d = np.cov(np.transpose(pwrmat_1d))

    return (mean_1d, std_1d, corrmat_1d, covmat_1d)


def summarize_1d_agg_pwrspec(pwr_1d, filename, corr_file=None,
                             apply_1d_transfer=None):
    r"""Summarize the 1D power spectrum from a list of one-dimensional power
    spectrum outputs.
    """
    (mean_1d, std_1d, corrmat_1d, covmat_1d) = agg_stat_1d_pwrspec(pwr_1d,
                                      apply_1d_transfer=apply_1d_transfer)

    # assume that they all have the same binning
    bin_left = pwr_1d[0]['bin_left']
    bin_center = pwr_1d[0]['bin_center']
    bin_right = pwr_1d[0]['bin_right']
    counts_histo = pwr_1d[0]['counts_histo']

    outfile = open(filename, "w")
    for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, mean_1d, std_1d):
        outfile.write(("%10.15g " * 6 + "\n") % specdata)
    outfile.close()

    bin_left_lt = np.log10(bin_left)
    bin_center_lt = np.log10(bin_center)
    bin_right_lt = np.log10(bin_right)

    if corr_file is not None:
        outfile = open(corr_file, "w")
        for xind in range(len(bin_center)):
            for yind in range(len(bin_center)):
                outstr = ("%10.15g " * 7 + "\n") % \
                        (bin_left_lt[xind], bin_center_lt[xind], bin_right_lt[xind], \
                         bin_left_lt[yind], bin_center_lt[yind], bin_right_lt[yind], \
                         corrmat_1d[xind, yind])
                outfile.write(outstr)

        outfile.close()

    retval = {}
    retval["mean_1d"] = mean_1d
    retval["std_1d"] = std_1d
    retval["covmat_1d"] = covmat_1d
    retval["bin_left"] = bin_left
    retval["bin_center"] = bin_center
    retval["bin_right"] = bin_right
    retval["counts_histo"] = counts_histo

    return retval


def summarize_1d_pwrspec(pwr_1d, filename):
    r"""Write out a 1d power spectrum
    """
    bin_left = pwr_1d['bin_left']
    bin_center = pwr_1d['bin_center']
    bin_right = pwr_1d['bin_right']

    outfile = open(filename, "w")
    for specdata in zip(bin_left, bin_center,
                        bin_right, pwr_1d['binavg']):
        outfile.write(("%10.15g " * 4 + "\n") % specdata)
    outfile.close()


def agg_stat_2d_pwrspec(pwr_2d, dataname='binavg'):
    pwrshp_2d = pwr_2d[0][dataname].shape
    n_runs = len(pwr_2d)

    pwrmat_2d = np.zeros((n_runs, pwrshp_2d[0], pwrshp_2d[1]))
    for index in range(n_runs):
        pwrmat_2d[index, :, :] = pwr_2d[index][dataname]

    mean_2d = np.mean(pwrmat_2d, axis=0)
    std_2d = np.std(pwrmat_2d, axis=0, ddof=1)

    return (mean_2d, std_2d)


def summarize_2d_pwrspec(pwr_2d, filename, dataname='binavg', resetnan=0.):
    r"""Write out a single pwrspec
    """
    bin_x_left = np.log10(pwr_2d['bin_x_left'])
    bin_x_center = np.log10(pwr_2d['bin_x_center'])
    bin_x_right = np.log10(pwr_2d['bin_x_right'])
    bin_y_left = np.log10(pwr_2d['bin_y_left'])
    bin_y_center = np.log10(pwr_2d['bin_y_center'])
    bin_y_right = np.log10(pwr_2d['bin_y_right'])

    outfile = open(filename, "w")
    reset_2d = copy.deepcopy(pwr_2d[dataname])
    reset_2d[np.isnan(reset_2d)] = resetnan
    for xind in range(len(bin_x_center)):
        for yind in range(len(bin_y_center)):
            outstr = ("%10.15g " * 7 + "\n") % \
                    (bin_x_left[xind], bin_x_center[xind], bin_x_right[xind], \
                     bin_y_left[yind], bin_y_center[yind], bin_y_right[yind], \
                     reset_2d[xind, yind])
            outfile.write(outstr)

    outfile.close()


def summarize_2d_agg_pwrspec(pwr_2d, filename, dataname='binavg', resetnan=0.):
    r"""Combine a list of 2D power runs and write out
    """
    (mean_2d, std_2d) = agg_stat_2d_pwrspec(pwr_2d, dataname=dataname)

    bin_x_left = np.log10(pwr_2d[0]['bin_x_left'])
    bin_x_center = np.log10(pwr_2d[0]['bin_x_center'])
    bin_x_right = np.log10(pwr_2d[0]['bin_x_right'])
    bin_y_left = np.log10(pwr_2d[0]['bin_y_left'])
    bin_y_center = np.log10(pwr_2d[0]['bin_y_center'])
    bin_y_right = np.log10(pwr_2d[0]['bin_y_right'])

    plot_mean_2d = copy.deepcopy(mean_2d)
    plot_std_2d = copy.deepcopy(std_2d)
    plot_mean_2d[np.isnan(mean_2d)] = resetnan
    plot_std_2d[np.isnan(std_2d)] = resetnan
    outfile = open(filename, "w")
    for xind in range(len(bin_x_center)):
        for yind in range(len(bin_y_center)):
            outstr = ("%10.15g " * 8 + "\n") % \
                    (bin_x_left[xind], bin_x_center[xind], bin_x_right[xind], \
                     bin_y_left[yind], bin_y_center[yind], bin_y_right[yind], \
                     plot_mean_2d[xind, yind], plot_std_2d[xind, yind])
            outfile.write(outstr)

    outfile.close()

    return mean_2d, std_2d


def summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d,
                      tag, outdir="./plot_data"):
    r"""Plot the 1D and 2D power spectra from a run
    """
    fileout = outdir + "/" + tag + "_from2d.dat"
    summarize_1d_pwrspec(pwr_1d_from_2d, fileout)

    fileout = outdir + "/" + tag + ".dat"
    summarize_1d_pwrspec(pwr_1d, fileout)

    fileout = outdir + "/" + tag + "_2d.dat"
    summarize_2d_pwrspec(pwr_2d, fileout, dataname="binavg")

    fileout = outdir + "/" + tag + "_2d_counts.dat"
    summarize_2d_pwrspec(pwr_2d, fileout, dataname="counts_histo")


def summarize_agg_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d,
                          tag, outdir="./plot_data", apply_1d_transfer=None):
    r"""call various power spectral aggregation functions to make 2D, 1D, etc.
    P(k)s to plot.
    """
    fileout = outdir + "/" + tag + "_avg_from2d.dat"
    corr_fileout = outdir + "/" + tag + "_corr_from2d.dat"
    agg_1d_pwrspec_f2d = summarize_1d_agg_pwrspec(pwr_1d_from_2d,
                                fileout,
                                corr_file=corr_fileout,
                                apply_1d_transfer=apply_1d_transfer)

    fileout = outdir + "/" + tag + "_avg.dat"
    corr_fileout = outdir + "/" + tag + "_corr.dat"
    agg_1d_pwrspec = summarize_1d_agg_pwrspec(pwr_1d, fileout,
                                             corr_file=corr_fileout,
                                apply_1d_transfer=apply_1d_transfer)

    fileout = outdir + "/" + tag + "_avg_2d.dat"
    summarize_2d_agg_pwrspec(pwr_2d, fileout, dataname="binavg")

    fileout = outdir + "/" + tag + "_avg_2d_counts.dat"
    summarize_2d_agg_pwrspec(pwr_2d, fileout, dataname="counts_histo")

    return agg_1d_pwrspec_f2d


def convert_2d_to_1d_driver(pwr_2d, counts_2d, bin_kx, bin_ky, bin_1d,
                            transfer=None):
    """take a 2D power spectrum and the counts matrix (number of modex per k
    cell) and return the binned 1D power spectrum
    pwr_2d is the 2D power
    counts_2d is the counts matrix
    bin_kx is the x-axis
    bin_ky is the x-axis
    bin_1d is the k vector over which to return the result
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

    if transfer is not None:
        orig_counts = copy.deepcopy(counts_2d_flat)
        trans_flat = transfer.flatten()
        pwr_2d_flat /= trans_flat
        counts_2d_flat *= trans_flat * trans_flat
        counts_2d_flat[np.isnan(counts_2d_flat)] = 0
        counts_2d_flat[np.isinf(counts_2d_flat)] = 0
        counts_2d_flat[orig_counts == 0] = 0

    count_pwr_prod = counts_2d_flat * pwr_2d_flat
    count_pwr_prod[np.isnan(count_pwr_prod)] = 0
    count_pwr_prod[np.isinf(count_pwr_prod)] = 0
    count_pwr_prod[counts_2d_flat == 0] = 0

    counts_histo = np.histogram(radius_flat, bin_1d,
                                weights=counts_2d_flat)[0]
    binsum_histo = np.histogram(radius_flat, bin_1d,
                                weights=count_pwr_prod)[0]

    binavg = binsum_histo / counts_histo.astype(float)

    return counts_histo, binavg


def convert_2d_to_1d(pwrspec2d_product, logbins=True,
                     bins=None, transfer=None):
    """if bins is not given, just use the x axis"""
    nxbins = len(pwrspec2d_product['bin_x_center'])
    bins_kx = np.zeros(nxbins + 1)
    bins_kx[0: -1] = pwrspec2d_product['bin_x_left']
    bins_kx[-1] = pwrspec2d_product['bin_x_right'][-1]

    nybins = len(pwrspec2d_product['bin_y_center'])
    bins_ky = np.zeros(nybins + 1)
    bins_ky[0: -1] = pwrspec2d_product['bin_y_left']
    bins_ky[-1] = pwrspec2d_product['bin_y_right'][-1]

    if bins is None:
        bins = bins_kx

    entry = {}
    (entry['counts_histo'], entry['binavg']) = \
                        convert_2d_to_1d_driver(pwrspec2d_product['binavg'],
                        pwrspec2d_product['counts_histo'],
                        pwrspec2d_product['bin_x_center'],
                        pwrspec2d_product['bin_y_center'],
                        bins, transfer=transfer)

    bin_left, bin_center, bin_right = binning.bin_edges(bins, log=logbins)
    entry['bin_left'] = bin_left
    entry['bin_center'] = bin_center
    entry['bin_right'] = bin_right

    return entry

if __name__ == '__main__':
    test_with_random(unitless=False)
