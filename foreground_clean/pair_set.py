r"""perform operations on sets of map pairs

Some improvements to consider:
    tag the outputs with the same 15hr_session etc as the data; right now these
    are fully specified by the directory
    pro to change: what if files move to a differently-named dir?
    con: what if tags grow sec_A_15hr_test1_4waysplit_etc.etc_etc._etc._with_B
    currently prefer more verbose directory names registered in the file db

    have the cleaning act on the weights -- this is probably only relevant when
    the full N^-1 rather than its diagonal is considered
"""
import os
import sys
import copy
import time
import scipy as sp
import numpy as np
from utils import file_tools as ft
from utils import data_paths as dp
from quadratic_products import corr_estimation as ce
from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from foreground_clean import map_pair
from multiprocessing import Process, current_process
from utils import batch_handler
# TODO: make map cleaning multiprocess; could also use previous cleaning, e.g.
# 5 modes to clean 10 modes = modes 5 to 10 (but no need to do this)
# TODO: move all magic strings to __init__ or params
# TODO: replace print with logging


params_init = {
               'SVD_root': None,
               'output_root': "local_test",
               'map1': 'GBT_15hr_map',
               'map2': 'GBT_15hr_map',
               'noise_inv1': 'GBT_15hr_map',
               'noise_inv2': 'GBT_15hr_map',
               'simfile': None,
               'sim_multiplier': 1.,
               'subtract_inputmap_from_sim': False,
               'subtract_sim_from_inputmap': False,
               'freq_list': (),
               'tack_on': None,
               'convolve': True,
               'factorizable_noise': True,
               'sub_weighted_mean': True,
               'regenerate_noise_inv': True,
               'modes': [10, 15],
               'no_weights': False
               }
prefix = 'fs_'


def wrap_find_weight(filename, regenerate=False):
    if regenerate:
        retval = find_weight(filename)
    else:
        retval = memoize_find_weight(filename)

    return batch_handler.repackage_kiyo(retval)


@batch_handler.memoize_persistent
def memoize_find_weight(filename):
    return find_weight(filename)


def find_weight(filename):
    r"""rather than read the full noise_inv and find its diagonal, cache the
    diagonal values.

    Note that the .info does not get shelved (class needs to be made
    serializeable). Return the info separately.
    """
    #print "loading noise: " + filename
    #noise_inv = algebra.make_mat(algebra.open_memmap(filename, mode='r'))
    #noise_inv_diag = noise_inv.mat_diag()
    # if optimal map:
    noise_inv_diag = algebra.make_vect(algebra.load(filename))

    return noise_inv_diag, noise_inv_diag.info


class PairSet():
    r"""Class to manage a set of map pairs
    """

    @batch_handler.log_timing
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        # recordkeeping
        self.pairs = {}
        self.pairs_parallel_track = {}
        self.pairlist = []
        self.datapath_db = dp.DataPath()

        self.params = params_dict
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

        self.freq_list = sp.array(self.params['freq_list'], dtype=int)
        self.output_root = self.datapath_db.fetch(self.params['output_root'],
                                                  tack_on=self.params['tack_on'])
        #self.output_root = self.params['output_root']
        if not os.path.isdir(self.output_root):
            os.mkdir(self.output_root)

        if self.params['SVD_root']:
            self.SVD_root = self.datapath_db.fetch(self.params['SVD_root'],
                                                   intend_write=True)
            print "WARNING: using %s to clean (intended?)" % self.SVD_root
        else:
            self.SVD_root = self.output_root

        # Write parameter file.
        parse_ini.write_params(self.params, self.output_root + 'params.ini',
                               prefix=prefix)

    @batch_handler.log_timing
    def execute(self, processes):
        r"""main call to execute the various steps in foreground removal"""
        self.load_pairs()
        self.preprocess_pairs()
        self.calculate_correlation()
        self.calculate_svd()

        mode_list_stop = self.params['modes']
        mode_list_start = copy.deepcopy(self.params['modes'])
        mode_list_start[1:] = mode_list_start[:-1]

        #self.uncleaned_pairs = copy.deepcopy(self.pairs)
        for (n_modes_start, n_modes_stop) in zip(mode_list_start,
                                             mode_list_stop):
            self.subtract_foregrounds(n_modes_start, n_modes_stop)
            self.save_data(n_modes_stop)

            # NOTE: if you use this you also need to copy the parallel pairs!
            #self.pairs = copy.deepcopy(self.uncleaned_pairs)

    @batch_handler.log_timing
    def load_pairs(self):
        r"""load the set of map/noise pairs specified by keys handed to the
        database. This sets up operations on the quadratic product
            Q = map1^T noise_inv1 B noise_inv2 map2
        """
        par = self.params
        (self.pairlist, pairdict) = dp.cross_maps(par['map1'], par['map2'],
                                             par['noise_inv1'],
                                             par['noise_inv2'],
                                             noise_inv_suffix=";noise_weight",
                                             verbose=False,
                                             db_to_use=self.datapath_db)

        for pairitem in self.pairlist:
            pdict = pairdict[pairitem]
            print "-" * 80
            dp.print_dictionary(pdict, sys.stdout,
                                key_list=['map1', 'noise_inv1',
                                          'map2', 'noise_inv2'])

            map1 = algebra.make_vect(algebra.load(pdict['map1']))
            map2 = algebra.make_vect(algebra.load(pdict['map2']))
            if par['simfile'] is not None:
                print "adding %s with multiplier %s" % (par['simfile'],
                                                        par['sim_multiplier'])

                sim = algebra.make_vect(algebra.load(par['simfile']))
                sim *= par['sim_multiplier']
                print sim.shape, map1.shape
            else:
                sim = algebra.zeros_like(map1)

            if not par['no_weights']:
                noise_inv1 = wrap_find_weight(pdict['noise_inv1'],
                                regenerate=par['regenerate_noise_inv'])

                noise_inv2 = wrap_find_weight(pdict['noise_inv2'],
                                regenerate=par['regenerate_noise_inv'])
            else:
                noise_inv1 = algebra.ones_like(map1)
                noise_inv2 = algebra.ones_like(map2)

            pair = map_pair.MapPair(map1 + sim, map2 + sim,
                                    noise_inv1, noise_inv2,
                                    self.freq_list)

            pair.set_names(pdict['tag1'], pdict['tag2'])

            pair.params = self.params
            self.pairs[pairitem] = pair

            if par['subtract_inputmap_from_sim'] or \
               par['subtract_sim_from_inputmap']:
                if par['subtract_inputmap_from_sim']:
                    pair_parallel_track = map_pair.MapPair(map1, map2,
                                                  noise_inv1, noise_inv2,
                                                  self.freq_list)

                if par['subtract_sim_from_inputmap']:
                    pair_parallel_track = map_pair.MapPair(sim, sim,
                                                  noise_inv1, noise_inv2,
                                                  self.freq_list)

                pair_parallel_track.set_names(pdict['tag1'], pdict['tag2'])
                pair_parallel_track.params = self.params
                self.pairs_parallel_track[pairitem] = pair_parallel_track


    @batch_handler.log_timing
    def preprocess_pairs(self):
        r"""perform several preparation tasks on the data
        1. convolve down to a common beam
        2. make the noise factorizable
        3. subtract the weighted mean from each freq. slice
        """
        if self.params["convolve"]:
            self.call_pairs("degrade_resolution")

        if self.params["factorizable_noise"]:
            self.call_pairs("make_noise_factorizable")

        if self.params["sub_weighted_mean"]:
            self.call_pairs("subtract_weighted_mean")

    @batch_handler.log_timing
    def calculate_correlation(self):
        r"""Note that multiprocessing's map() is more elegant than Process,
        but fails for handing in complex map_pair objects
        """
        process_list = []
        for pairitem in self.pairlist:
            filename = self.output_root
            filename += "foreground_corr_pair_%s.pkl" % pairitem
            multi = Process(target=wrap_corr, args=([self.pairs[pairitem],
                            filename]), name=pairitem)

            process_list.append(multi)

            multi.start()

        for process in process_list:
            process.join()

    @batch_handler.log_timing
    def calculate_svd(self):
        r"""calculate the SVD of all pairs"""
        for pairitem in self.pairlist:
            filename = self.output_root
            filename_corr = filename + "foreground_corr_pair_%s.pkl" % pairitem
            filename_svd = filename + "SVD_pair_%s.pkl" % pairitem
            print filename_corr
            if os.access(filename_corr, os.F_OK):
                print "SVD loading corr. functions: " + filename
                (freq_cov, counts) = ft.load_pickle(filename_corr)

                # (vals, modes1, modes2)
                svd_info = ce.get_freq_svd_modes(freq_cov, len(self.freq_list))
                ft.save_pickle(svd_info, filename_svd)
            else:
                print "ERROR: in SVD, correlation functions not loaded"
                sys.exit()

    @batch_handler.log_timing
    def subtract_foregrounds(self, n_modes_start, n_modes_stop):
        for pairitem in self.pairlist:
            filename_svd = "%s/SVD_pair_%s.pkl" % (self.SVD_root, pairitem)
            print "subtracting %d to %d modes from %s using %s" % (n_modes_start, \
                                                                n_modes_stop, \
                                                                pairitem, \
                                                                filename_svd)

            # svd_info: 0 is vals, 1 is modes1 (left), 2 is modes2 (right)
            svd_info = ft.load_pickle(filename_svd)

            self.pairs[pairitem].subtract_frequency_modes(
                                    svd_info[1][n_modes_start:n_modes_stop],
                                    svd_info[2][n_modes_start:n_modes_stop])

            if self.params['subtract_inputmap_from_sim'] or \
               self.params['subtract_sim_from_inputmap']:
                self.pairs_parallel_track[pairitem].subtract_frequency_modes(
                                        svd_info[1][n_modes_start:n_modes_stop],
                                        svd_info[2][n_modes_start:n_modes_stop])

    @batch_handler.log_timing
    def save_data(self, n_modes):
        prodmap_list = []
        weight_list = []

        n_modes = "%dmodes" % n_modes
        for pairitem in self.pairlist:
            pair = self.pairs[pairitem]
            (tag1, tag2) = (pair.map1_name, pair.map2_name)
            clnoise = "cleaned_noise_inv"
            map1_file = "%s/sec_%s_cleaned_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, tag2, n_modes)
            map2_file = "%s/sec_%s_cleaned_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, tag1, n_modes)
            noise_inv1_file = "%s/sec_%s_%s_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, clnoise, tag2, n_modes)
            noise_inv2_file = "%s/sec_%s_%s_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, clnoise, tag1, n_modes)
            modes1_file = "%s/sec_%s_modes_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, tag2, n_modes)
            modes2_file = "%s/sec_%s_modes_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, tag1, n_modes)

            if self.params['subtract_inputmap_from_sim'] or \
               self.params['subtract_sim_from_inputmap']:
                map1 = pair.map1 - self.pairs_parallel_track[pairitem].map1

                map2 = pair.map2 - self.pairs_parallel_track[pairitem].map2
            else:
                map1 = copy.deepcopy(pair.map1)
                map2 = copy.deepcopy(pair.map2)

            prodmap_list.append(map1 * pair.noise_inv1)
            prodmap_list.append(map2 * pair.noise_inv2)
            weight_list.append(pair.noise_inv1)
            weight_list.append(pair.noise_inv2)

            algebra.save(map1_file, map1)
            algebra.save(map2_file, map2)
            algebra.save(noise_inv1_file, pair.noise_inv1)
            algebra.save(noise_inv2_file, pair.noise_inv2)
            algebra.save(modes1_file, pair.left_modes)
            algebra.save(modes2_file, pair.right_modes)

        cumulative_product = algebra.zeros_like(prodmap_list[0])
        cumulative_weight = algebra.zeros_like(prodmap_list[0])
        for mapind in range(0, len(prodmap_list)):
            cumulative_product += prodmap_list[mapind]
            cumulative_weight += weight_list[mapind]

        algebra.compressed_array_summary(cumulative_weight, "weight map")
        algebra.compressed_array_summary(cumulative_product, "product map")

        cumulative_weight[cumulative_weight < 1.e-20] = 0.
        cumulative_product[cumulative_weight < 1.e-20] = 0.

        newmap = cumulative_product / cumulative_weight

        # if the new map is nan or inf, set it and the wieghts to zero
        nan_array = np.isnan(newmap)
        newmap[nan_array] = 0.
        cumulative_product[nan_array] = 0.
        cumulative_weight[nan_array] = 0.
        inf_array = np.isinf(newmap)
        newmap[inf_array] = 0.
        cumulative_product[inf_array] = 0.
        cumulative_weight[inf_array] = 0.
        algebra.compressed_array_summary(newmap, "new map")
        algebra.compressed_array_summary(cumulative_product, "final map * weight")
        algebra.compressed_array_summary(cumulative_weight, "final weight map")

        combined = "combined_clean"
        combined_map_file = "%s/%s_map_%s.npy" % \
                            (self.output_root, combined, n_modes)
        combined_weight_file = "%s/%s_weight_%s.npy" % \
                            (self.output_root, combined, n_modes)
        combined_product_file = "%s/%s_product_%s.npy" % \
                            (self.output_root, combined, n_modes)
        combined_ones_file = "%s/%s_ones_%s.npy" % \
                            (self.output_root, combined, n_modes)

        algebra.save(combined_map_file, newmap)
        algebra.save(combined_product_file, cumulative_product)
        algebra.save(combined_weight_file, cumulative_weight)
        algebra.save(combined_ones_file, algebra.ones_like(newmap))

    # Service functions ------------------------------------------------------
    def call_pairs(self, call):
        r"""call some operation on the map pairs"""
        for pairitem in self.pairlist:
            pair = self.pairs[pairitem]
            try:
                method_to_call = getattr(pair, call)
            except AttributeError:
                print "ERROR: %s missing call %s" % (pairitem, call)
                sys.exit()

            print "calling %s() on pair %s" % (call, pairitem)
            method_to_call()

        if self.params['subtract_inputmap_from_sim'] or \
           self.params['subtract_sim_from_inputmap']:
            for pairitem in self.pairlist:
                pair_parallel_track = self.pairs_parallel_track[pairitem]
                try:
                    method_to_call_parallel_track = getattr(pair_parallel_track, call)
                except AttributeError:
                    print "ERROR: %s missing call %s" % (pairitem, call)
                    sys.exit()

                print "calling %s() on pair %s" % (call, pairitem)
                method_to_call_parallel_track()


def wrap_corr(pair, filename):
    r"""Do the correlation for a map_pair `pair`.
    Correlations in the `pair` (map_pair type) are saved to `filename`
    """
    name = current_process().name
    (freq_cov, counts) = pair.freq_covariance()
    ft.save_pickle((freq_cov, counts), filename)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        PairSet(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'
