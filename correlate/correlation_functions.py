"""various functions to calculate and bin correlation functions"""
import scipy as sp
import copy
from numpy import linalg

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


def get_freq_svd_modes(corr, n_modes):
    r"""Same as get freq eigenmodes, but treats left and right maps
    separatly with an SVD.

    Parameters
    ----------
    corr: 3D array
        The correlation. Only the 1st lag is used.
    n_modes: int
        The number of modes wanted.

    Returns
    -------
    s: 1D array
        The amplitude of the modes. length = `n_modes`
    left_vectors, right_vectors: 2D array
        The first `n_modes` from the svd.

    """
    u_matrix, singular_values, v_matrix = linalg.svd(corr[:, :, 0])
    v_matrix = v_matrix.T
    sorted_singular_values = list(singular_values)
    sorted_singular_values.sort()
    left_vectors = []
    right_vectors = []
    for mode_index in range(n_modes):
        ind, = sp.where(abs(singular_values) ==
                        sorted_singular_values[-mode_index - 1])
        if len(ind) > 1:
            raise NotImplementedError('2 eigenvalues bitwise equal.')

        left_vectors.append(u_matrix[:, ind[0]])
        right_vectors.append(v_matrix[:, ind[0]])

    return singular_values, left_vectors, right_vectors


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
