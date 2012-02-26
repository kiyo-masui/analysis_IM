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
from utils import file_tools as ft
from utils import data_paths as dp
from correlate import corr_estimation as ce
from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from correlate import map_pair
from multiprocessing import Process, current_process
# TODO: make map cleaning multiprocess; could also use previous cleaning, e.g.
# 5 modes to clean 10 modes = modes 5 to 10 (but no need to do this)
# TODO: move all magic strings to __init__ or params
# TODO: replace print with logging

params_init = {
               'output_root': "local_test",
               # Options of saving.
               # What frequencies to correlate:
               'map1': 'GBT_15hr_map',
               'map2': 'GBT_15hr_map',
               'noise_inv1': 'GBT_15hr_map',
               'noise_inv2': 'GBT_15hr_map',
               'simfile': '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr/sim_beam_0.npy',
               'freq_list': (),
               # Angular lags at which to calculate the correlation.  Upper
               # edge bins in degrees.
               'lags': (0.1, 0.2),
               'convolve': True,
               'factorizable_noise': True,
               'sub_weighted_mean': True,
               'modes': [10, 15],
               'no_weights': False
               }
prefix = 'fs_'


class PairSet(ft.ClassPersistence):
    r"""Class to manage a set of map pairs
    """

    def __init__(self, parameter_file_or_dict=None):
        # recordkeeping
        self.pairs = {}
        self.pairs_nosim = {}
        self.pairlist = []
        self.noisefiledict = {}
        self.datapath_db = dp.DataPath()


        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=prefix)

        self.freq_list = sp.array(self.params['freq_list'], dtype=int)
        self.lags = sp.array(self.params['lags'])
        self.output_root = self.datapath_db.fetch(self.params['output_root'],
                                                  intend_write=True)

        # Write parameter file.
        kiyopy.utils.mkparents(self.output_root)
        parse_ini.write_params(self.params, self.output_root + 'params.ini',
                               prefix=prefix)

    def execute(self):
        r"""main call to execute the various steps in foreground removal"""
        self.calculate_corr_svd()
        self.clean_maps(freestanding=False) # data are already loaded

    def calculate_corr_svd(self):
        r""" "macro" which finds the correlation functions for the pairs of
        given maps, then the SVD.
        """
        self.load_pairs(regenerate=True)
        self.preprocess_pairs()
        # If the correlation is already calculated, this can be commented out
        # and the SVD will load the previously-calculated correlations
        self.calculate_correlation()
        self.calculate_svd()

    def clean_maps(self, freestanding=True):
        r""" "macro" which open the data, does pre-processing, loads the
        pre-calculated SVD modes, subtracts them, and saves the data.
        This is the second stage. If `freestanding` assume that the data have
        not been loaded into pairs and pre-processed.
        """
        if freestanding:
            self.load_pairs(regenerate=False)  # use previous diag(N^-1)
            self.preprocess_pairs()

        self.uncleaned_pairs = copy.deepcopy(self.pairs)
        for n_modes in self.params['modes']:
            # clean self.pairs and save its components
            self.subtract_foregrounds(n_modes)
            self.save_data(n_modes)
            # reset the all the pair data
            self.pairs = copy.deepcopy(self.uncleaned_pairs)

    def load_pairs(self, regenerate=True):
        r"""load the set of map/noise pairs specified by keys handed to the
        database. This sets up operations on the quadratic product
            Q = map1^T noise_inv1 B noise_inv2 map2
        """
        par = self.params
        (self.pairlist, pairdict) = dp.cross_maps(par['map1'], par['map2'],
                                             par['noise_inv1'],
                                             par['noise_inv2'],
                                             verbose=False)

        for pairitem in self.pairlist:
            pdict = pairdict[pairitem]
            print "-" * 80
            dp.print_dictionary(pdict, sys.stdout,
                                key_list=['map1', 'noise_inv1',
                                          'map2', 'noise_inv2'])

            map1 = algebra.make_vect(algebra.load(pdict['map1']))
            map2 = algebra.make_vect(algebra.load(pdict['map2']))
            sim = algebra.make_vect(algebra.load(par['simfile']))

            if not par['no_weights']:
                noise_inv1 = self.process_noise_inv(pdict['noise_inv1'],
                                                    regenerate=regenerate)

                noise_inv2 = self.process_noise_inv(pdict['noise_inv2'],
                                                    regenerate=regenerate)
            else:
                noise_inv1 = algebra.ones_like(map1)
                noise_inv2 = algebra.ones_like(map2)

            pair = map_pair.MapPair(map1 + sim, map2 + sim,
                                    noise_inv1, noise_inv2,
                                    self.freq_list)

            pair_nosim = map_pair.MapPair(map1, map2,
                                    noise_inv1, noise_inv2,
                                    self.freq_list)

            pair.set_names(pdict['tag1'], pdict['tag2'])
            pair_nosim.set_names(pdict['tag1'], pdict['tag2'])

            pair.lags = self.lags
            pair.params = self.params
            pair_nosim.lags = self.lags
            pair_nosim.params = self.params
            self.pairs[pairitem] = pair
            self.pairs_nosim[pairitem] = pair_nosim

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

    def calculate_correlation(self):
        r"""Note that multiprocessing's map() is more elegant than Process,
        but fails for handing in complex map_pair objects
        """
        process_list = []
        for pairitem in self.pairlist:
            filename = self.output_root
            filename += "foreground_corr_pair_%s.pkl" % pairitem
            #wrap_corr(self.pairs[pairitem], filename)
            multi = Process(target=wrap_corr, args=([self.pairs[pairitem],
                            filename]), name=pairitem)

            process_list.append(multi)

            multi.start()

        for process in process_list:
            process.join()

    def calculate_svd(self):
        r"""calculate the SVD of all pairs"""
        for pairitem in self.pairlist:
            filename = self.output_root
            filename_corr = filename + "foreground_corr_pair_%s.pkl" % pairitem
            filename_svd = filename + "SVD_pair_%s.pkl" % pairitem
            print filename_corr
            if os.access(filename_corr, os.F_OK):
                print "SVD loading corr. functions: " + filename
                (corr, counts) = ft.load_pickle(filename_corr)

                # (vals, modes1, modes2)
                svd_info = ce.get_freq_svd_modes(corr, len(self.freq_list))
                ft.save_pickle(svd_info, filename_svd)
            else:
                print "ERROR: in SVD, correlation functions not loaded"
                sys.exit()

    def subtract_foregrounds(self, n_modes):
        r"""call subtract_frequency_modes on the maps with the modes as found
        in the svd, removing the first `n_modes`
        """
        for pairitem in self.pairlist:
            print "subtracting %d modes from %s" % (n_modes, pairitem)
            filename_svd = "%s/SVD_pair_%s.pkl" % (self.output_root, pairitem)
            # svd_info: 0 is vals, 1 is modes1 (left), 2 is modes2 (right)
            svd_info = ft.load_pickle(filename_svd)

            self.pairs[pairitem].subtract_frequency_modes(
                                    svd_info[1][:n_modes],
                                    svd_info[2][:n_modes])

            self.pairs_nosim[pairitem].subtract_frequency_modes(
                                    svd_info[1][:n_modes],
                                    svd_info[2][:n_modes])

    def save_data(self, n_modes):
        ''' Save the all of the clean data.
        '''
        # Make sure folder is there.
        if not os.path.isdir(self.output_root):
            os.mkdir(self.output_root)

        n_modes = "%dmodes" % n_modes
        for pairitem in self.pairlist:
            pair = self.pairs[pairitem]
            pair_nosim = self.pairs_nosim[pairitem]
            (tag1, tag2) = (pair.map1_name, pair.map2_name)
            map1_file = "%s/sec_%s_cleaned_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, tag2, n_modes)
            map2_file = "%s/sec_%s_cleaned_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, tag1, n_modes)
            noise_inv1_file = "%s/sec_%s_cleaned_noise_inv_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, tag2, n_modes)
            noise_inv2_file = "%s/sec_%s_cleaned_noise_inv_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, tag1, n_modes)
            modes1_file = "%s/sec_%s_modes_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, tag2, n_modes)
            modes2_file = "%s/sec_%s_modes_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, tag1, n_modes)

            algebra.save(map1_file, pair.map1 - pair_nosim.map1)
            algebra.save(map2_file, pair.map2 - pair_nosim.map2)
            algebra.save(noise_inv1_file, pair.noise_inv1)
            algebra.save(noise_inv2_file, pair.noise_inv2)
            algebra.save(modes1_file, pair.left_modes)
            algebra.save(modes2_file, pair.right_modes)

    # Service functions ------------------------------------------------------
    def call_pairs(self, call):
        r"""call some operation on the map pairs"""
        for pairitem in self.pairlist:
            pair = self.pairs[pairitem]
            pair_nosim = self.pairs_nosim[pairitem]
            try:
                method_to_call = getattr(pair, call)
                method_to_call_nosim = getattr(pair_nosim, call)
            except AttributeError:
                print "ERROR: %s missing call %s" % (pairitem, call)
                sys.exit()

            print "calling %s() on pair %s" % (call, pairitem)
            method_to_call()
            method_to_call_nosim()

    # TODO: this could probably replaced with a memoize
    def process_noise_inv(self, filename, regenerate=True):
        r"""buffer reading the noise inverse files for speed and also
        save to a file in the intermediate output path.

        If the cached file exists as an intermediate product, load it else
        produce it.
        """
        if filename not in self.noisefiledict:
            basename = filename.split("/")[-1].split(".npy")[0]
            filename_diag = "%s/%s_diag.npy" % \
                           (self.output_root, basename)
            exists = os.access(filename_diag, os.F_OK)
            if exists and not regenerate:
                print "loading pre-diagonalized noise: " + filename_diag
                self.noisefiledict[filename] = algebra.make_vect(
                                                algebra.load(filename_diag))
            else:
                print "loading noise: " + filename
                noise_inv = algebra.make_mat(
                                    algebra.open_memmap(filename, mode='r'))
                self.noisefiledict[filename] = noise_inv.mat_diag()
                #self.noisefiledict[filename] = algebra.make_vect(
                #                               algebra.load(filename))
                algebra.save(filename_diag, self.noisefiledict[filename])

        return copy.deepcopy(self.noisefiledict[filename])


def wrap_corr(pair, filename):
    r"""Do the correlation for a map_pair `pair`.
    Correlations in the `pair` (map_pair type) are saved to `filename`
    """
    name = current_process().name
    print "starting at %s on pair %s => %s" % (name, time.asctime(), filename)
    (corr, counts) = pair.correlate(pair.lags, speedup=True)
    ft.save_pickle((corr, counts), filename)
    print "%s finished at %s" % (name, time.asctime())


if __name__ == '__main__':
    if len(sys.argv) == 2:
        PairSet(str(sys.argv[1])).execute()
        #PairSet(str(sys.argv[1])).clean_maps()
    else:
        print 'Need one argument: parameter file name.'
