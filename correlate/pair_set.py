import os
import sys
import copy
import time
import scipy as sp
from utils import file_tools as ft
from utils import data_paths as dp
from correlate import correlation_functions as cf
from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from correlate import map_pair
from multiprocessing import Process, current_process
# TODO: replace output_root with intermediate product directory
# TODO: move all magic strings to __init__ or params
# TODO: replace print with logging
# TODO: make function which subtracts 0-50
# TODO: have the cleaning act on the weights too?
# TODO: tag the outputs with the same 15hr_etc etc. as the data?

params_init = {
               'output_root': "./testoutput/",
               # Options of saving.
               # What frequencies to correlate:
               'map1': 'GBT_15hr_map',
               'map2': 'GBT_15hr_map',
               'noise_inv1': 'GBT_15hr_map',
               'noise_inv2': 'GBT_15hr_map',
               'freq_list': (),
               # Angular lags at which to calculate the correlation.  Upper
               # edge bins in degrees.
               'lags': (0.1, 0.2),
               'convolve': True,
               'factorizable_noise': True,
               'sub_weighted_mean': True,
               'modes': 10,
               'no_weights': False
               }
prefix = 'fs_'


class PairSet(ft.ClassPersistence):
    r"""Class to manage a set of map pairs
    """

    def __init__(self, parameter_file_or_dict=None):
        # recordkeeping
        self.pairs = {}
        self.pairlist = []
        self.noisefiledict = {}

        self.params = parse_ini.parse(parameter_file_or_dict, params_init,
                                      prefix=prefix)

        self.freq_list = sp.array(self.params['freq_list'], dtype=int)
        self.lags = sp.array(self.params['lags'])

        # Write parameter file.
        kiyopy.utils.mkparents(self.params['output_root'])
        parse_ini.write_params(self.params,
                               self.params['output_root'] + 'params.ini',
                               prefix=prefix)

    def execute(self):
        r"""main call to execute the various steps in foreground removal"""
        self.load_pairs(regenerate=True)
        self.preprocess_pairs()
        #self.calculate_correlation()
        self.calculate_svd()

    def clean_maps(self):
        self.load_pairs(regenerate=False)  # use previous diag(N^-1)
        self.preprocess_pairs()
        n_modes = self.params['modes']
        self.subtract_foregrounds(n_modes)
        self.save_data()

    def load_pairs(self, regenerate=True):
        r"""load the set of map/noise pairs specified by keys handed to the
        database. This sets up operations on the quadratic product
            Q = map1^T noise_inv1 B noise_inv2 map2
        """
        par = self.params
        (self.pairlist, pairdict) = dp.cross_maps(par['map1'], par['map2'],
                                             par['noise_inv1'], par['noise_inv2'],
                                             verbose=False)

        for pairitem in self.pairlist:
            pdict = pairdict[pairitem]
            print "-"*80
            dp.print_dictionary(pdict, sys.stdout,
                                key_list=['map1', 'noise_inv1', 'map2', 'noise_inv2'])

            map1 = algebra.make_vect(algebra.load(pdict['map1']))
            map2 = algebra.make_vect(algebra.load(pdict['map2']))

            if not par['no_weights']:
                noise_inv1 = self.process_noise_inv(pdict['noise_inv1'],
                                                    regenerate=regenerate)

                noise_inv2 = self.process_noise_inv(pdict['noise_inv2'],
                                                    regenerate=regenerate)
            else:
                noise_inv1 = algebra.ones_like(map1)
                noise_inv2 = algebra.ones_like(map2)

            pair = map_pair.MapPair(map1, map2,
                                    noise_inv1, noise_inv2,
                                    self.freq_list)

            pair.set_names(pdict['tag1'], pdict['tag2'])

            pair.lags = self.lags
            pair.params = self.params
            self.pairs[pairitem] = pair

    def preprocess_pairs(self):
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
            filename = self.params['output_root']
            filename += "map_pair_for_freq_slices_fore_corr_%s.pkl" % pairitem
            multi = Process(target=wrap_corr, args=([self.pairs[pairitem],
                            filename]), name=pairitem)

            process_list.append(multi)

            multi.start()

        for process in process_list:
            process.join()

    def calculate_svd(self):
        r"""calculate the SVD of all pairs"""
        for pairitem in self.pairlist:
            filename = self.params['output_root']
            filename_corr = filename + "foreground_corr_pair_%s.pkl" % pairitem
            filename_svd = filename + "SVD_pair_%s.pkl" % pairitem
            print filename_corr
            if os.access(filename_corr, os.F_OK):
                print "SVD loading corr. functions: " + filename
                (corr, counts) = ft.load_pickle(filename_corr)

                # (vals, modes1, modes2)
                svd_info = cf.get_freq_svd_modes(corr, len(self.freq_list))
                ft.save_pickle(svd_info, filename_svd)
            else:
                print "ERROR: in SVD, correlation functions not loaded"
                sys.exit()

    def subtract_foregrounds(self, n_modes):
        for pairitem in self.pairlist:
            print "subtracting %d modes from %s" % (n_modes, pairitem)
            filename = self.params['output_root']
            filename_svd = filename + "SVD_pair_%s.pkl" % pairitem
            svd_info = ft.load_pickle(filename_svd)

            vals = svd_info[0]
            (all_modes1, all_modes2) = (svd_info[1], svd_info[2])

            self.pairs[pairitem].subtract_frequency_modes(
                                    svd_info[1][:n_modes],
                                    svd_info[2][:n_modes])

    def save_data(self):
        ''' Save the all of the clean data.
        '''
        out_root = self.params['output_root']
        # Make sure folder is there.
        if not os.path.isdir(out_root):
            os.mkdir(out_root)

        for pairitem in self.pairlist:
            pair = self.pairs[pairitem]
            (tag1, tag2) = (pair.map1_name, pair.map2_name)
            map1_savename = "%s/sec_%s_cleaned_clean_map_I_with_%s.npy" % \
                            (out_root, tag1, tag2)
            map2_savename = "%s/sec_%s_cleaned_clean_map_I_with_%s.npy" % \
                            (out_root, tag2, tag1)
            noise_inv1_savename = "%s/sec_%s_cleaned_noise_inv_I_with_%s.npy" % \
                            (out_root, tag1, tag2)
            noise_inv2_savename = "%s/sec_%s_cleaned_noise_inv_I_with_%s.npy" % \
                            (out_root, tag2, tag1)
            modes1_savename = "%s/sec_%s_modes_clean_map_I_with_%s.npy" % \
                            (out_root, tag1, tag2)
            modes2_savename = "%s/sec_%s_modes_clean_map_I_with_%s.npy" % \
                            (out_root, tag2, tag1)

            algebra.save(map1_savename, pair.map1)
            algebra.save(map2_savename, pair.map2)
            algebra.save(noise_inv1_savename, pair.noise_inv1)
            algebra.save(noise_inv2_savename, pair.noise_inv2)
            algebra.save(modes1_savename, pair.left_modes)
            algebra.save(modes2_savename, pair.right_modes)

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
                           (self.params['output_root'], basename)
            exists = os.access(filename_diag, os.F_OK)
            if exists and not regenerate:
                print "loading diagonal noise: " + filename_diag
                self.noisefiledict[filename] = algebra.make_vect(
                                                algebra.load(filename_diag))
            else:
                print "loading noise: " + filename
                noise_inv = algebra.make_mat(
                                    algebra.open_memmap(filename, mode='r'))
                self.noisefiledict[filename] = noise_inv.mat_diag()
                #self.noisefiledict[filename] = sp.zeros((10,10))
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
    import sys
    if len(sys.argv) == 2:
        #PairSet(str(sys.argv[1])).execute()
        PairSet(str(sys.argv[1])).clean_maps()
    else:
        print 'Need one argument: parameter file name.'
