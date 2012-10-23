"""various functions to calculate and bin correlation functions"""
import numpy as np
import scipy as sp
from numpy import linalg
import time
import copy


def corr_est(map1, map2, noise1, noise2, freq1, freq2,
             lags=(), speedup=False, verbose=False):
    r"""Calculate the cross correlation function of the maps.

    The cross correlation function is a function of f1, f2 and angular lag.
    The angular lag bins are passed, all pairs of frequencies are
    calculated.

    Parameters
    ----------
    lags: array like
        Angular lags bins (upper side bin edges).
    speedup: boolean
        Speeds up the correlation. This works fine, yes? Should be the
        normal way if so.

    Returns
    -------
    corr: array
        The correlation between 2 maps.
    counts: array
        The weighting of the correlation based on the maps' weights.

    """
    map1_ra = map1.get_axis('ra')
    map2_ra = map2.get_axis('ra')
    map1_dec = map1.get_axis('dec')
    map2_dec = map2.get_axis('dec')

    input_map1 = map1[freq1, :, :]
    input_map2 = map2[freq2, :, :]
    input_noise1 = noise1[freq1, :, :]
    input_noise2 = noise2[freq2, :, :]

    # Noise weight
    input_map1 *= input_noise1
    input_map2 *= input_noise2

    nlags = len(lags)
    nfreq = len(freq1)
    corr = sp.zeros((nfreq, nfreq, nlags), dtype=float)
    counts = sp.zeros(corr.shape, dtype=float)
    # Noting that if DEC != 0, then a degree of RA is less than a degree.
    ra_fact = sp.cos(sp.pi * map1.info['dec_centre'] / 180.0)

    # Calculate the pairwise lags.
    dra = (map1_ra[:, None] - map2_ra[None, :]) * ra_fact
    ddec = map1_dec[:, None] - map2_dec[None, :]
    lag = dra[:, None, :, None] ** 2 + ddec[None, :, None, :] ** 2
    lag = sp.sqrt(lag)
    # Bin this up.
    lag_inds = sp.digitize(lag.flatten(), lags)

    if speedup:
        print "Starting Correlation (sparse version)"

        (nr1, nd1) = (len(map1_ra), len(map1_dec))
        (nr2, nd2) = (len(map2_ra), len(map2_dec))
        (r1ind, d1ind) = (sp.arange(nr1), sp.arange(nd1))
        (r2ind, d2ind) = (sp.arange(nr2), sp.arange(nd2))
        ra1_pairind = r1ind.repeat(nr2 * nd1 * nd2)
        ra2_pairind = sp.tile(r2ind.repeat(nd2), (1, nr1 * nd1)).flatten()
        dec1_pairind = sp.tile(d1ind.repeat(nr2 * nd2), (1, nr1)).flatten()
        dec2_pairind = sp.tile(d2ind, (1, nr1 * nr2 * nd1)).flatten()

        # precalculate the pair indices for a given lag
        # could also imagine calculating the map slices here
        posmaskdict = {}
        for klag in range(nlags):
            mask = (lag_inds == klag)
            posmaskdict[repr(klag)] = (ra1_pairind[mask],
                                       ra2_pairind[mask],
                                       dec1_pairind[mask],
                                       dec2_pairind[mask])

        for if1 in range(len(freq1)):
            for jf2 in range(len(freq2)):
                start = time.time()

                data1 = input_map1[if1, :, :]
                data2 = input_map2[jf2, :, :]
                weights1 = input_noise1[if1, :, :]
                weights2 = input_noise2[jf2, :, :]

                for klag in range(nlags):
                    (r1m, r2m, d1m, d2m) = posmaskdict[repr(klag)]
                    dprod = data1[r1m, d1m] * data2[r2m, d2m]
                    wprod = weights1[r1m, d1m] * weights2[r2m, d2m]
                    corr[if1, jf2, klag] += sp.sum(dprod)
                    counts[if1, jf2, klag] += sp.sum(wprod)

                if verbose:
                    print if1, jf2, (time.time() - start)
                    print counts[if1, jf2, :]
    else:
        print "Starting Correlation (full version)"
        for if1 in range(len(freq1)):
            for jf2 in range(len(freq2)):
                start = time.time()
                # Calculate the pairwise products.
                data1 = input_map1[if1, :, :]
                data2 = input_map2[jf2, :, :]
                weights1 = input_noise1[if1, :, :]
                weights2 = input_noise2[jf2, :, :]
                dprod = data1[..., None, None] * data2[None, None, ...]
                wprod = weights1[..., None, None] * \
                        weights2[None, None, ...]
                for klag in range(nlags):
                    mask = (lag_inds == klag)
                    corr[if1, jf2, klag] += sp.sum(dprod.flatten()[mask])
                    counts[if1, jf2, klag] += sp.sum(wprod.flatten()[mask])

                if verbose:
                    print if1, jf2, (time.time() - start)
                    print counts[if1, jf2, :]

    mask = (counts < 1e-20)
    counts[mask] = 1
    corr /= counts
    corr[mask] = 0
    counts[mask] = 0

    return corr, counts


def get_corr_and_std_3d(corr_list):
    '''Return the average correlation and its standard deviation.

    Parameters
    ----------
    corr_list: list of 3D arrays
        This list contains a correlation for each `map_pair`.
        len must be > 0.

    Returns
    -------
    corr_avg: 3D array
        The average of the correlations in `corr_list`.
    corr_std: 3D array
        The standard deviation on each point `in corr_avg`.

    '''
    # Get average.
    corr_sum = copy.deepcopy(corr_list[0])
    for corr_index in range(1, len(corr_list)):
        corr_sum += corr_list[corr_index]

    corr_avg = corr_sum / float(len(corr_list))
    # Get std. Note it will be all 0 if only one pair was used.
    corr_x_size, corr_y_size, corr_z_size = corr_list[0].shape
    corr_std = sp.zeros((corr_x_size, corr_y_size, corr_z_size))

    for corr_x_index in range(0, corr_x_size):
        for corr_y_index in range(0, corr_y_size):
            for corr_z_index in range(0, corr_z_size):
                value_list = []
                for corr in corr_list:
                    value_list.append(corr[corr_x_index, corr_y_index,
                                           corr_z_index])

                corr_std[corr_x_index, corr_y_index, corr_z_index] = \
                         sp.std(value_list)

    return corr_avg, corr_std


def subtract_modes_corr(corr, n_modes):
    r"""Get the modes at subtract them directly from the correlation.

    Similar to get_freq_svd_modes, but the modes are subtracted and from the
    correlation and the cleaned correlation is returned.

    Parameters
    ----------
    corr: 3D array
        The correlation. Only the 1st lag is used.
    n_modes: int
        The number of modes wanted.

    Returns
    -------
    corr: 3D array
        The input `corr` with its modes directly subtracted.

    Notes
    -----
    This works only for the 1st lag.
    (Was that it? There was something holding this back.)

    """

    corr = copy.deepcopy(corr)

    singular_values, left_vectors, right_vectors = \
                get_freq_svd_modes(corr, n_modes)

    for mode_index in range(n_modes):
        corr -= singular_values[mode_index] * \
                left_vectors[mode_index][:, None, None] * \
                right_vectors[mode_index][None, :, None]

    return corr


def normalize_corr(corr):
    r"""Return the normalized correlation along the diagonal.

    Paramters
    ---------
    corr: 3D array
        The correlation.

    Returns
    -------
    corr_norm: 3D array
        The normalized `corr`. The normalization is such that the values
        along the diagonal (f=f_prime) are 1.

    """
    # Get dimensions.
    freqs, freqs_prime, lags = corr.shape
    corr_norm = sp.zeros((freqs, freqs_prime, lags))
    # Each angular lag is separate.
    for lag in range(0, lags):

        # At each pixel, get the geometric mean of the 2 freqs along the
        # diagonal [ie. at C[f, f, lag] and C[f_prime, f_prime, lag]). Divide
        # the pixel by that value to make it normalized.
        for freq in range(0, freqs):
            print lag,
            print freq

            for freq_prime in range(0, freqs_prime):
                value = corr[freq, freq, lag] * \
                        corr[freq_prime, freq_prime, lag]

                factor = sp.sqrt(value)

                corr_norm[freq, freq_prime, lag] = \
                    corr[freq, freq_prime, lag] / factor

    return corr_norm


def rebin_corr_freq_lag(corr, freq1, freq2=None, weights=None, nfbins=20,
                        return_fbins=False):
    r"""Collapses frequency pair correlation function to frequency lag.

    Basically this constructs the 2D correlation function.

    Parameters
    ----------
    corr: 3D array
        Covariance matrix which is a function of frequency and frequency prime
        and angular lag.
    freq1, freq2: tuple of floats
        The REAL frequencies. ie. 744000Hz, not 0, 1, 2...
        freq2 is used if using a map at a different redshift, but we haven't
        looked at this yet.
    weights: 3D array
        The weights of the correlation. It is found in pair.counts right now.
    nfbins: int
        How many lag bins out you go in frequency. A higher number means a
        more accurate result at high lag.
    return_fbins: bool
        If `True`, `fbins` is returned.

    Returns
    -------
    out_corr: 2D array
        `corr` from before but now only in terms of frequency lag and
        angular lag.
    out_weights: 2D array
        `weights` from before but now in 2D. The weights for `out_corr`
    fbins: 1D array
        The frequency lags in terms of Hz. Much like how `lags` in the rest of
        this module is angular lag in degrees.

    """

    if freq2 is None:
        freq2 = freq1
    # Default is equal weights.
    if weights is None:
        weights = sp.ones_like(corr)
    corr = corr * weights

    nf1 = corr.shape[0]
    nf2 = corr.shape[1]
    nlags = corr.shape[2]
    # Frequency bin size.
    delta_freq = min(abs(sp.diff(freq1)))
    # Frequency bin upper edges.
    fbins = (sp.arange(nfbins) + 0.5) * delta_freq
    # Allowcate memory for outputs.
    out_corr = sp.zeros((nfbins, nlags))
    out_weights = sp.zeros((nfbins, nlags))

    # Loop over all frequency pairs and bin by lag.
    for freq1_index in range(nf1):
        for freq2_index in range(nf2):
            f_lag = abs(freq1[freq1_index] - freq2[freq2_index])
            bin_ind = sp.digitize([f_lag], fbins)[0]
            if bin_ind < nfbins:
                out_corr[bin_ind, :] += corr[freq1_index, freq2_index, :]
                out_weights[bin_ind, :] += weights[freq1_index, freq2_index, :]

    # Normalize dealing with 0 weight points explicitly.
    bad_inds = out_weights < 1.0e-20
    out_weights[bad_inds] = 1.0
    out_corr /= out_weights
    out_weights[bad_inds] = 0.0
    out_corr[bad_inds] = 0.0

    if return_fbins:
        return out_corr, out_weights, fbins - delta_freq * 0.5
    else:
        return out_corr, out_weights


def collapse_correlation_1d(corr, f_lags, a_lags, weights=None):
    r"""Takes a 2D correlation function and collapses to a 1D correlation
    function.

    Parameters
    ----------
    corr: 2D array
        Covariance matrix in terms of frequency lag and angular lag.
        The first output from `rebin_corr_freq_lag` right now.
    f_lags: 1D array
        The frequency lags in terms of Hz.
        The third output from `rebin_corr_freq_lag` right now.
    a_lags: 1D array
        The angular lags in terms of degrees.
    weights: 2D array
        The weights of `corr`.
        The second output from `rebin_corr_freq_lag` right now.

    Returns
    -------
    out_corr: 1D array
        The 1D autocorrelation.
    out_weights:
        The weights for `out_corr`.
    x_axis: tuple of 3 1D arrays
        `x_axis[1]` is the x - values that correspond to `out_corr`.
        `x_axis[0]` and `x_axis[2]` are the left and rightmost points
         covered by each lag bin.

    Notes
    -----
    `a_lags` are not the same as the lags from the .ini file.
    The lags from the .ini file are the right side of each lag bin,
    but you want the centre of the bin when you plot.
    To get the right values, you must do: (ask Eric or Liviu)
        lags = sp.array(F.params['lags'])
        a_lags = copy.deepcopy(lags)
        a_lags[0] = 0
        a_lags[1:] -= sp.diff(lags)/2.0
    """

    if corr.ndim != 2:
        msg = "Must start with a 2D correlation function."
        raise ValueError(msg)

    if len(f_lags) != corr.shape[0] or len(a_lags) != corr.shape[1]:
        msg = ("corr.shape must be (len(f_lags), len(a_lags)).  Passed: "
               + repr(corr.shape) + " vs (" + repr(len(f_lags)) + ", "
               + repr(len(a_lags)) + ").")

        raise ValueError(msg)

    if weights is None:
        weights = sp.ones_like(corr)

    corr = corr * weights
    # Hard code conversion factors to MPc/h for now.
    a_fact = 34.0  # Mpc/h per degree at 800MHz.
    f_fact = 4.5   # Mpc/h per MHz at 800MHz.
    # Hard code lags in MPc/h.
    #nbins = 10
    nbins = 15
    lags = sp.empty(nbins)
    lags[0] = 2.0
    lags[1] = 4.0

    for bin_index in range(2, nbins):
        lags[bin_index] = 1.5 * lags[bin_index - 1]

    # Calculate the total 1D lags.
    separation = a_lags
    separation = (a_fact * separation[sp.newaxis, :]) ** 2
    separation = separation + (f_fact * f_lags[:, sp.newaxis] / 1.0e6) ** 2
    separation = sp.sqrt(separation)

    # Initialize memory for outputs.
    out_corr = sp.zeros(nbins)
    out_weights = sp.zeros(nbins)

    # Rebin.
    for lag_index in range(separation.shape[0]):
        bin_inds = sp.digitize(separation[lag_index, :], lags)
        for bin_index in range(nbins):
            out_corr[bin_index] += sp.sum(corr[lag_index,
                                               bin_inds == bin_index])
            out_weights[bin_index] += sp.sum(weights[lag_index,
                                                     bin_inds == bin_index])
    # Normalize.
    bad_inds = out_weights < 1.0e-20
    out_weights[bad_inds] = 1.0
    out_corr /= out_weights
    out_weights[bad_inds] = 0.0

    # Make real lags to be returned.
    x_left = sp.empty(nbins)
    x_left[0] = 0
    x_left[1:] = lags[:-1]
    x_right = lags
    x_centre = (x_right + x_left) / 2.0

    return out_corr, out_weights, (x_left, x_centre, x_right)
