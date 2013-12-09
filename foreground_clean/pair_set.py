r"""perform operations on sets of map pairs
"""
import os
import sys
import copy
import time
import scipy as sp
import numpy as np
from utils import data_paths as dp
from foreground_clean import find_modes
from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from foreground_clean import map_pair
from multiprocessing import Process, current_process
from utils import batch_handler
import h5py

params_init = {
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
               'tack_on_input': None,
               'tack_on_output': None,
               'convolve': True,
               'factorizable_noise': True,
               'sub_weighted_mean': True,
               'svd_filename': None,
               'modes': [10, 15]
               }
prefix = 'fs_'


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
        self.tack_on_input = self.params['tack_on_input']
        self.output_root = self.datapath_db.fetch(self.params['output_root'],
                                            tack_on=self.params['tack_on_output'])

        #self.output_root = self.params['output_root']
        print "foreground cleaning writing to output root", self.output_root

        if not os.path.isdir(self.output_root):
            os.mkdir(self.output_root)

        if self.params['svd_filename'] is not None:
            self.svd_filename = self.params['svd_filename']
            print "WARNING: using %s to clean (intended?)" % self.svd_filename
        else:
            self.svd_filename = self.output_root + "/" + "SVD.hd5"

        # Write parameter file.
        parse_ini.write_params(self.params, self.output_root + 'params.ini',
                               prefix=prefix)

    @batch_handler.log_timing
    def execute(self, processes):
        r"""main call to execute the various steps in foreground removal"""
        self.load_pairs()
        self.preprocess_pairs()

        if self.params['svd_filename'] is not None:
            print "WARNING: skipping correlation/SVD and using existing"
        else:
            self.calculate_correlation()

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
                                             tack_on=self.tack_on_input,
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

            noise_inv1 = algebra.make_vect(algebra.load(pdict['noise_inv1']))
            noise_inv2 = algebra.make_vect(algebra.load(pdict['noise_inv2']))

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
        r"""find the covariance in frequency space, take SVD"""
        svd_data_out = h5py.File(self.svd_filename, "w")

        cov_out = svd_data_out.create_group("cov")
        counts_out = svd_data_out.create_group("cov_counts")
        svd_vals_out = svd_data_out.create_group("svd_vals")
        svd_modes1_out = svd_data_out.create_group("svd_modes1")
        svd_modes2_out = svd_data_out.create_group("svd_modes2")

        nfreq = len(self.freq_list)

        for pairitem in self.pairlist:
            (freq_cov, counts) = self.pairs[pairitem].freq_covariance()
            cov_out[pairitem] = freq_cov
            counts_out[pairitem] = counts

            # (vals, modes1, modes2)
            svd_info = find_modes.get_freq_svd_modes(freq_cov, nfreq)
            svd_vals_out[pairitem] = svd_info[0]
            svd_modes1_out[pairitem] = svd_info[1]
            svd_modes2_out[pairitem] = svd_info[2]

        svd_data_out.close()

    @batch_handler.log_timing
    def subtract_foregrounds(self, n_modes_start, n_modes_stop):
        r"""take the SVD modes from above and clean each LOS with them"""
        svd_data = h5py.File(self.svd_filename, "r")
        svd_modes1 = svd_data["svd_modes1"]
        svd_modes2 = svd_data["svd_modes2"]

        for pairitem in self.pairlist:
            print "subtracting %d to %d modes from %s using %s" % \
                    (n_modes_start, n_modes_stop, pairitem, self.svd_filename)

            svd_modes1_pair = svd_modes1[pairitem].value
            svd_modes2_pair = svd_modes2[pairitem].value
            svd_modes1_use = svd_modes1_pair[n_modes_start: n_modes_stop]
            svd_modes2_use = svd_modes2_pair[n_modes_start: n_modes_stop]

            self.pairs[pairitem].subtract_frequency_modes(svd_modes1_use,
                                                          svd_modes2_use)

            if self.params['subtract_inputmap_from_sim'] or \
               self.params['subtract_sim_from_inputmap']:
                self.pairs_parallel_track[pairitem].subtract_frequency_modes(
                                        svd_modes1_use, svd_modes2_use)

        svd_data.close()

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


if __name__ == '__main__':
    if len(sys.argv) == 2:
        PairSet(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'
