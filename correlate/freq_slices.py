"""Clean the radio cubes
"""

import copy
import multiprocessing
import time
import os
import cPickle
import scipy as sp
from numpy import linalg
from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from correlate import map_pair
from utils import io_wrap
from correlate import correlation_functions as cf

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
               'freq_list': (),
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
               'factorizable_noise': False,
               'no_weights': False
               }
prefix = 'fs_'


class NewSlices(object):
    r"""Pipeline module that scripts together the cross correlation of maps.

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

        # main derived quantities:
        self.pairs = None
        self.corr_std = None
        self.fore_pairs = None
        self.svd_info_list = None
        self.corr_final = None

    def execute(self):
        '''Clean the maps of foregrounds, save the results, and get the
        autocorrelation.'''

        params = self.params
        freq_list = sp.array(params['freq_list'], dtype=int)
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
                            params['file_middles'][map_index] +
                            params['input_end_map'])

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
                                                noise_inv1, noise_inv2,
                                                freq_list)

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
                                                          len(freq_list))
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
                io_wrap.save_pickle(self.svd_info_list, params['svd_file'])
        else:
            # The first correlation and svd has been skipped.
            # This means you already have the modes so you can just load
            # them from file.
            self.svd_info_list = io_wrap.load_pickle(params['svd_file'])
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
        self.save_data(self, params['save_maps'], params['save_noises'],
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

        self.corr_final, self.corr_std = cf.get_corr_and_std_3d(corr_list)

        if params['pickle_slices']:
            pickle_slices(self)

        return

    def save_data(self, save_maps=False, save_noises=False, save_modes=False):
        ''' Save the all of the clean data.

        Saves the cleaned data and modes to the output directory specified
        in the parameters of the configuration file.

        Parameters
        ----------
        save_maps: bool
            Save the cleaned clean maps to output directory if `True`.
        save_noises: bool
            Save the cleaned noise invs to output directory if `True`.
        save_modes: bool
            Save what was subtracted from the maps to output directory if `True`.

        '''
        # Get the output directory and filenames.
        out_root = self.params['output_root']
        # Make sure folder is there.
        if not os.path.isdir(out_root):
            os.mkdir(out_root)

        for pair in self.pairs:
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

            if save_maps:
                algebra.save(map1_savename, pair.map1)
                algebra.save(map2_savename, pair.map2)
            if save_noises:
                algebra.save(noise_inv1_savename, pair.noise_inv1)
                algebra.save(noise_inv2_savename, pair.noise_inv2)
            if save_modes:
                algebra.save(left_modes_savename, pair.left_modes)
                algebra.save(right_modes_savename, pair.right_modes)


def multiproc(pair, save_dir, pair_number, final):
    r"""Do the correlation for a map_pair `pair`.

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
    r"""Call the 1st or 2nd correlation on the map_pair.

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
        final correlation is being done. Don't mix these up!

    """
    if final:
        pair.corr, pair.counts = pair.correlate(lags, speedup=True)
    else:
        pair.fore_corr, pair.fore_counts = pair.correlate(lags, speedup=True)


# TODO: phase this out because of annoying import dependencies
def pickle_slices(slice_obj):
    r"""Save the NewSlices object with ALL the info to file.

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
    else:
        print 'Need one argument: parameter file name.'
