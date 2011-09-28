"""Program that calculates the correlation function across frequency slices.
"""

import copy
import multiprocessing
import time
import os
import cPickle

import scipy as sp
from numpy import linalg
import matplotlib.pyplot as plt

from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from correlate import map_pair

params_init = {
               # IO:
               'input_root': './',
               # The unique part of every fname
               'file_middles': ("testfile", ),
               'input_end_map': "_map.fits",
               'input_end_noise': "_noise.fits",
               'output_root': "./testoutput",
               # Options of saving.
               'save_maps': False,
               'save_noises': False,
               'save_modes': False,
               'pickle_slices': False,
               # What frequencies to correlate:
               'freq': (),
               # Angular lags at which to calculate the correlation.  Upper
               # edge bins in degrees.
               'lags': (0.1, 0.2),
               'convolve': False,
               'sub_weighted_mean': False,
               'modes': 10,
               # How long to run the program.
               'first_pass_only': False,
               'skip_fore_corr': False,
               # saving an svd file used when wanting to skip fore corr.
               'save_svd_info': False,
               # Saving and loading svd stuff occurs to the same file but
               # never at the same time.
               'svd_file': '',  # Must be a cPickle file
               'make_plots': False,
               'factorizable_noise': False,
               'no_weights': False
               }
prefix = 'fs_'


class NewSlices(object):
    """Pipeline module that scripts together the cross correlation of maps.

    Parameters
    ----------
    parameter_file_or_dict: file or dict
        Loads parameters for execution and stores in a dictionary.

    Attributes
    ----------
    params: dict
        A dictionary containing all the information from correlate_slices.ini
    fore_pairs: list of map_pair
        Keeps track of all pair information until just after 1st correlation.
    pairs: list of map_pair
        Keeps track of all pairs. The ordering of the pairs in this list is
        based on the double for loop that goes over file_middles.
    svd_info_list: list of svd_infos
        Contains the svd_info for each `map_pair`.
        svd_info has 3 elements - `vals`, `all_modes1`, and `all_modes2` (see
        map_pair documention for what they are).
    corr_final: 3D array
        The average correlation from all pairs.
    corr_std: 3D array
        The standard deviation of `corr_final`.

    """

    def __init__(self, parameter_file_or_dict=None):
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                 prefix=prefix)

    def execute(self):
        '''Clean the maps of foregrounds, save the results, and get the
        autocorrelation.'''

        params = self.params
        freq = sp.array(params['freq'], dtype=int)
        lags = sp.array(params['lags'])
        # Write parameter file.
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        # Get the map data from file as well as the noise inverse.
        if len(params['file_middles']) == 1:
            fmid_name = params['file_middles'][0]
            params['file_middles'] = (fmid_name, fmid_name)
        if len(params['file_middles']) >= 2:
            # Deal with multiple files.
            num_maps = len(params['file_middles'])
            maps = []
            noise_invs = []
            # Load all maps and noises once.
            for map_index in range(0, num_maps):
                map_file = (params['input_root'] +
                    params['file_middles'][map_index] + params['input_end_map'])
                print "Loading map %d of %d." % (map_index + 1, num_maps)
                map_in = algebra.make_vect(algebra.load(map_file))
                maps.append(map_in)
                if not params["no_weights"]:
                    noise_file = (params['input_root'] +
                                  params['file_middles'][map_index] +
                                  params['input_end_noise'])
                    print "Loading noise %d of %d." % (map_index + 1, num_maps)
                    noise_inv = algebra.make_mat(
                                    algebra.open_memmap(noise_file, mode='r'))
                    noise_inv = noise_inv.mat_diag()
                else:
                    noise_inv = algebra.ones_like(map_in)
                noise_invs.append(noise_inv)

            pairs = []
            # Make pairs with deepcopies to not make mutability mistakes.
            for map1_index in range(0, num_maps):
                for map2_index in range(0, num_maps):
                    if (map2_index > map1_index):
                        map1 = copy.deepcopy(maps[map1_index])
                        map2 = copy.deepcopy(maps[map2_index])
                        noise_inv1 = copy.deepcopy(noise_invs[map1_index])
                        noise_inv2 = copy.deepcopy(noise_invs[map2_index])
                        pair = map_pair.MapPair(map1, map2,
                            noise_inv1, noise_inv2, freq)
                        pair.lags = lags
                        pair.params = params
                        # Keep track of the names of maps in pairs so
                        # it knows what to save later.
                        pair.set_names(params['file_middles'][map1_index],
                                       params['file_middles'][map2_index])
                        pairs.append(pair)
            num_map_pairs = len(pairs)
            print "%d map pairs created from %d maps." % (len(pairs), num_maps)
        # Hold a reference in self.
        self.pairs = pairs

        # Get maps/ noise inv ready for running.
        if params["convolve"]:
            for pair in pairs:
                pair.degrade_resolution()
        if params['factorizable_noise']:
            for pair in pairs:
                pair.make_noise_factorizable()
        if params['sub_weighted_mean']:
            for pair in pairs:
                pair.subtract_weighted_mean()

        self.pairs = pairs
        # Since correlating takes so long, if you already have the svds
        # you can skip this first correlation [since that's all it's really
        # for and it is the same no matter how many modes you want].
        # Note: map_pairs will not have anything saved in 'fore_corr' if you
        # skip this correlation.
        if not params['skip_fore_corr']:
            # Correlate the maps with multiprocessing. Note that the
            # correlations are saved to file separately then loaded in
            # together because that's (one way) how multiprocessing works.
            fore_pairs = []
            processes_list = []
            for pair_index in range(0, num_map_pairs):
                # Calls 1 multiproc (which governs the correlating) for each
                # pair on a new CPU so you can have all pairs working at once.
                multi = multiprocessing.Process(target=multiproc,
                                                args=([pairs[pair_index],
                                                       params['output_root'],
                                                       pair_index, False]))
                processes_list.append(multi)
                multi.start()

            # Waits for all correlations to finish before continuing.
            while True in [multi.is_alive() for multi in processes_list]:
                print "processing"
                time.sleep(5)
            # just to be safe
            time.sleep(1)

            # Load the correlations and save them to each pair. The pairs that
            # got passed to multiproc are not the same ones as ones in
            # self.pairs, so this must be done to have actual values.
            print "Loading map pairs back into program."
            file_name = params['output_root']
            file_name += "map_pair_for_freq_slices_fore_corr_"
            for count in range(0, num_map_pairs):
                print "Loading correlation for pair %d" % (count)
                pickle_handle = open(file_name + str(count) + ".pkl", "r")
                correlate_results = cPickle.load(pickle_handle)
                pairs[count].fore_corr = correlate_results[0]
                pairs[count].fore_counts = correlate_results[1]
                fore_pairs.append(pairs[count])
                pickle_handle.close()
            self.fore_pairs = copy.deepcopy(fore_pairs)
            # With this, you do not need fore_pairs anymore.
            self.pairs = copy.deepcopy(fore_pairs)
            print "gung ho!"

            pairs = self.pairs

            # Get foregrounds.

            # svd_info_list keeps track of all of the modes of all maps in
            # all pairs. This means if you want to subract a different number
            # of modes for the same maps/noises/frequencies, you have the modes
            # already saved and do not need to run the first correlation again.
            svd_info_list = []
            for pair in pairs:
                vals, modes1, modes2 = get_freq_svd_modes(pair.fore_corr,
                                                   len(freq))
                pair.vals = vals
                # Save ALL of the modes for reference.
                pair.all_modes1 = modes1
                pair.all_modes2 = modes2
                svd_info = (vals, modes1, modes2)
                svd_info_list.append(svd_info)
                # Save only the modes you want to subtract.
                n_modes = params['modes']
                pair.modes1 = modes1[:n_modes]
                pair.modes2 = modes2[:n_modes]
            self.svd_info_list = svd_info_list
            self.pairs = pairs
            if params['save_svd_info']:
                save_svd_info(self.svd_info_list, params['svd_file'])
        else:
            # The first correlation and svd has been skipped.
            # This means you already have the modes so you can just load
            # them from file.
            self.svd_info_list = load_svd_info(params['svd_file'])
            # Set the svd info to the pairs.
            for i in range(0, len(pairs)):
                svd_info = self.svd_info_list[i]
                pairs[i].vals = svd_info[0]
                pairs[i].all_modes1 = svd_info[1]
                pairs[i].all_modes2 = svd_info[2]
                n_modes = params['modes']
                pairs[i].modes1 = svd_info[1][:n_modes]
                pairs[i].modes2 = svd_info[2][:n_modes]
            self.pairs = pairs

        # Subtract foregrounds.
        for pair_index in range(0, len(pairs)):
            pairs[pair_index].subtract_frequency_modes(pairs[pair_index].modes1,
                pairs[pair_index].modes2)

        # Save cleaned clean maps, cleaned noises, and modes.
        save_data(self, params['save_maps'], params['save_noises'],
            params['save_modes'])

        # Finish if this was just first pass.
        if params['first_pass_only']:
            self.pairs = pairs
            return

        # Correlate the cleaned maps.
        # Here we could calculate the power spectrum instead eventually.

        # Same as previous multiproc, but this time, you get the final
        # correlation for each pair.
        temp_pair_list = []
        processes_list = []
        for pair_index in range(0, num_map_pairs):
            multi = multiprocessing.Process(target=multiproc,
                                            args=([pairs[pair_index],
                                                   params['output_root'],
                                                   pair_index, True]))
            processes_list.append(multi)
            multi.start()

        while True in [multi.is_alive() for multi in processes_list]:
            print "processing"
            time.sleep(5)
        # just to be safe
        time.sleep(1)
        print "Loading map pairs back into program."
        file_name = params['output_root']
        file_name += "map_pair_for_freq_slices_corr_"
        for count in range(0, num_map_pairs):
            print "Loading correlation for pair %d" % (count)
            pickle_handle = open(file_name + str(count) + ".pkl", "r")
            correlate_results = cPickle.load(pickle_handle)
            pairs[count].corr = correlate_results[0]
            pairs[count].counts = correlate_results[1]
            temp_pair_list.append(pairs[count])
            pickle_handle.close()
        self.pairs = copy.deepcopy(temp_pair_list)
        print "gung ho!"

        # Get the average correlation and its standard deviation.
        corr_list = []
        for pair in self.pairs:
            corr_list.append(pair.corr)
        self.corr_final, self.corr_std = get_corr_and_std_3d(corr_list)

        # Used be useful a long time ago but now plotting has moved to it's
        # own folder. It really moved up in the world this summer.
        if params["make_plots"]:
            print "Plots not supported for multiple pairs."
 #           self.make_plots()

        if params['pickle_slices']:
            pickle_slices(self)

        return


def multiproc(pair, save_dir, pair_number, final):
    """Do the correlation for a map_pair `pair`.

    The correlation is then saved to a numbered file corresponding to
    `pair_number`. This is done since this function gets multiprocessed and
    must have a common information ground somewhere.

    Parameters
    ----------
    pair: map_pair
        A pair with 2 maps that will be correlated.
    save_dir: str
        The full pathname to where the correlations will be saved.
    pair_number: int
        A number from 0 to (# of pairs)-1. Must be unique for each pair.
    final: bool
        If `True`, then the fore correlation is being done. If `False`, then
        the final correlation is being done. Dont mix these up!

    """

    print "I am starting."
    # Having a function looks nicer than an if statement.
    control_correlation(pair, pair.lags, final)
    file_name = save_dir
    # Make sure folder is there.
    if not os.path.isdir(file_name):
        os.mkdir(file_name)
    # Save the correlation and its weights. Note that using an int for
    # each pair is enough because of the for loop that passes each pair into
    # multiproc.
    if final:
        file_name += "map_pair_for_freq_slices_corr_" + \
                        str(pair_number) + ".pkl"
        to_save = (pair.corr, pair.counts)
    else:
        file_name += "map_pair_for_freq_slices_fore_corr_" + \
                        str(pair_number) + ".pkl"
        to_save = (pair.fore_corr, pair.fore_counts)
    pickle_handle = open(file_name, "w")
    print "Writing to: ",
    print file_name
    cPickle.dump(to_save, pickle_handle)
    pickle_handle.close()
    return


def control_correlation(pair, lags, final):
    """Call the 1st or 2nd correlation on the map_pair.

    Saves the correlation and its weight in `pair`.
    For multiprocessing since it cannot deal with return values.

    Parameters
    ----------
    pair: map_pair
        A pair with 2 maps that will be correlated.
    lags: tuple of ints
        The angular lags.
    final: bool
        If True, then the fore correlation is being done. If False, then the
        final correlation is being done. Dont mix these up!

    """

    if final:
        pair.corr, pair.counts = pair.correlate(lags, speedup=True)
    else:
        pair.fore_corr, pair.fore_counts = pair.correlate(lags, speedup=True)


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
    """Same as get freq eigenmodes, but treats left and right maps
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
    """Get the modes at subtract them directly from the correlation.

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


def plot_svd(vals):
    """Deprecated.

    Plots the svd values and prints out some statistics."""

    n_vals = len(vals)
    plt.semilogy(abs(sp.sort(-vals / n_vals)), marker='o',
                 linestyle='None')
    print 'Mean noise: ', sp.sum(vals) / n_vals
    print 'Largest eigenvalues/vals: ',
    print sp.sort(vals / n_vals)[-10:]


def normalize_corr(corr):
    """Return the normalized correlation along the diagonal.

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
                value = corr[freq, freq, lag] * corr[freq_prime, freq_prime, lag]
                factor = sp.sqrt(value)
                corr_norm[freq, freq_prime, lag] = corr[freq, freq_prime, lag] / factor
    return corr_norm


def rebin_corr_freq_lag(corr, freq1, freq2=None, weights=None, nfbins=20,
                        return_fbins=False):
    """Collapses frequency pair correlation function to frequency lag.

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
    """Takes a 2D correlation function and collapses to a 1D correlation
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


def save_data(slice_obj, save_maps=False, save_noises=False, save_modes=False):
    ''' Save the all of the clean data.

    Saves the cleaned data and modes to the output directory specified
    in the parameters of the configuration file.

    Parameters
    ----------

    slice_obj: NewSlices
        `slice_obj` contains ALL of the data.
    save_maps: bool
        Save the cleaned clean maps to output directory if `True`.
    save_noises: bool
        Save the cleaned noise invs to output directory if `True`.
    save_modes: bool
        Save what was subtracted from the maps to output directory if `True`.

    '''
    # Get the output directory and filenames.
    out_root = slice_obj.params['output_root']
    # Make sure folder is there.
    if not os.path.isdir(out_root):
        os.mkdir(out_root)
    for pair in slice_obj.pairs:
        map1_savename = out_root + pair.map1_name + \
            "_cleaned_clean_map_I_with_" + pair.map2_code + ".npy"
        map2_savename = out_root + pair.map2_name + \
            "_cleaned_clean_map_I_with_" + pair.map1_code + ".npy"
        noise_inv1_savename = out_root + pair.map1_name + \
            "_cleaned_noise_inv_I_with_" + pair.map2_code + ".npy"
        noise_inv2_savename = out_root + pair.map2_name + \
            "_cleaned_noise_inv_I_with_" + pair.map1_code + ".npy"
        left_modes_savename = out_root + pair.map1_name + \
            "_modes_clean_map_I_with_" + pair.map2_code + ".npy"
        right_modes_savename = out_root + pair.map2_name + \
            "_modes_clean_map_I_with_" + pair.map1_code + ".npy"
        # Save data.
        if save_maps:
            algebra.save(map1_savename, pair.map1)
            algebra.save(map2_savename, pair.map2)
        if save_noises:
            algebra.save(noise_inv1_savename, pair.noise_inv1)
            algebra.save(noise_inv2_savename, pair.noise_inv2)
        if save_modes:
            algebra.save(left_modes_savename, pair.left_modes)
            algebra.save(right_modes_savename, pair.right_modes)


def save_svd_info(svd_info_list, svd_file):
    """cPickle the `svd_info_list` to file with name `svd_file`.

    Parameters
    ----------
    svd_info_list: list of svd_infos
        See `NewSlices` documentation for what this is.
    svd_file: str
        The full pathname of where to save `svd_info_list`.

    """
    # Check saving folder exists.
    pickle_out_root = os.path.dirname(svd_file)
    if not os.path.isdir(pickle_out_root):
        os.mkdir(pickle_out_root)
    # Save.
    pickle_handle = open(svd_file, 'w')
    cPickle.dump(svd_info_list, pickle_handle)
    pickle_handle.close()


def load_svd_info(svd_file):
    """Return `svd_info_list` saved in file with name `svd_file`.

    Parameters
    ----------
    svd_file: str
        The full pathname of where to load from.

    Returns
    -------
    svd_info_list: list of svd_infos
        See `NewSlices` documentation for what this is.

    """
    pickle_handle = open(svd_file, 'r')
    svd_info_list = cPickle.load(pickle_handle)
    pickle_handle.close()
    return svd_info_list


def pickle_slices(slice_obj):
    """Save the NewSlices object with ALL the info to file.

    Pickle slice_obj to the output directory from the ini file.

    Parameters
    ----------
    slice_obj: NewSlices
        The object to get pickled. The output directory is saved in itself.

    Notes
    -----
    If you pickle slice_obj from outide of this class, to use it,
    you can import: 'from correlate import freq_slices as fs'
    But since the pickle has been done from inside this class, you must:
    'from correlate.freq_slices import *'

    """
    # Check folder exists.
    out_root = slice_obj.params['output_root']
    if not os.path.isdir(out_root):
        os.mkdir(out_root)
    # Save.
    pickle_file = out_root + 'New_Slices_object.pkl'
    pickle_handle = open(pickle_file, 'w')
    cPickle.dump(slice_obj, pickle_handle)
    pickle_handle.close()


# For running this module from the command line
if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        NewSlices(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else:
        NewSlices().execute()
