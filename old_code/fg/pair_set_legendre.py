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
from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from foreground_clean import map_pair
from multiprocessing import Process, current_process
from utils import batch_handler

import pair_set
import replace_svd_mode as rpmode

# TODO: make map cleaning multiprocess; could also use previous cleaning, e.g.
# 5 modes to clean 10 modes = modes 5 to 10 (but no need to do this)
# TODO: move all magic strings to __init__ or params
# TODO: replace print with logging
# Special parameters that rarely need to get used
# no_weights: do not weight before finding nu-nu' covariance (False by default)
# SVD_root: use the SVD from another cleaning run on this run (None by default)
# regenerate_noise_inv: do not use memoize to save previous diag(N^-1) calc


params_init = {
               'SVD_root': None,
               'SVD_file': None,
               'output_root': "local_test",
               'pairlist' : None,
               'pairdict' : None,

               'map1': 'GBT_15hr_map',
               'map2': 'GBT_15hr_map',
               'noise_inv1': 'GBT_15hr_map',
               'noise_inv2': 'GBT_15hr_map',
               'calc_diagnal' : False,
               'simnum' : None,
               'simfile1': None,
               'simfile2': None,
               'sim_multiplier': 1.,
               'subtract_inputmap_from_sim': False,
               'subtract_sim_from_inputmap': False,
               'subtract_realmap_from_sim': False,
               'realmap_dir' : '',
               'freq_list1': (),
               'freq_list2': (),
               'freq_n_all1': 256,
               'freq_n_all2': -1,
                # in deg: (unused)
               'tack_on': None,
               'convolve': True,
               'clip_weight_percent': None,
               'degrade_factor': 1.1,
               'factorizable_noise': True,
               'sub_weighted_mean': True,
               'regenerate_noise_inv': True,
               'modes': [10, 15],
               'good_modes': 0, # number of modes remained. 
                                # set it <=0, if no modes need to be replaced. 
               'no_weights': False,
               'save_section': True,
               }
prefix = 'fs_'



class PairSet_LegendreSVD(pair_set.PairSet):
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

        self.freq_list1 = sp.array(self.params['freq_list1'], dtype=int)
        if len(self.params['freq_list2']) == 0:
            self.freq_list2 = self.freq_list1
        else:
            self.freq_list2 = sp.array(self.params['freq_list2'], dtype=int)

    @batch_handler.log_timing
    def subtract_foregrounds(self, n_modes_start, n_modes_stop):
        for pairitem in self.pairlist:
            if self.params['SVD_file'] != None:
                filename_svd = "%s/%s" % (self.SVD_root, self.params['SVD_file'])
            else:
                filename_svd = "%s/SVD_pair_%s.pkl" % (self.SVD_root, pairitem)
            print "subtracting %d to %d modes from %s using %s" % (n_modes_start, \
                                                                n_modes_stop, \
                                                                pairitem, \
                                                                filename_svd)

            # svd_info: 0 is vals, 1 is modes1 (left), 2 is modes2 (right)
            svd_info = ft.load_pickle(filename_svd)

            # replace svd modes by Legendre polymodial
            if self.params['good_modes'] > 0:
                good_modes = self.params['good_modes']
                print "replace SVD mode from %d to the end" % good_modes
                freq_n_all1 = self.params['freq_n_all1']
                if freq_n_all2 == -1:
                    freq_n_all2 = freq_n_all1
                else:
                    freq_n_all2 = self.params['freq_n_all2']
                mode_n = len(svd_info[1])
                svd_info_all = [np.zeros(shape=(mode_n, freq_n_all1)), 
                                np.zeros(shape=(mode_n, freq_n_all2))]
                for i in range(mode_n):
                    np.put(svd_info_all[0][i], self.freq_list1, svd_info[1][i])
                    np.put(svd_info_all[1][i], self.freq_list2, svd_info[2][i])
                print 'replace mode from %d to the end' %good_modes

                if good_modes == 1:
                    svd_info = (svd_info[0],
                                rpmode.replace_modes(svd_info_all[0], 
                                                     good_modes,),
                                rpmode.replace_modes(svd_info_all[1], 
                                                     good_modes,))
                else:
                    svd_info = (svd_info[0],
                                rpmode.replace_modes(svd_info_all[0], 
                                                     good_modes, 
                                                     m=n_modes_stop,
                                                     weight=svd_info_all[0][0]),
                                rpmode.replace_modes(svd_info_all[1], 
                                                     good_modes, 
                                                     m=n_modes_stop,
                                                     weight=svd_info_all[1][0]))
                #svd_info = (svd_info[0],
                #    rpmode.replace_modes(svd_info_all[0], good_modes, m=n_modes_stop),
                #    rpmode.replace_modes(svd_info_all[1], good_modes, m=n_modes_stop))

                svd_info = (svd_info[0],
                    np.take(svd_info[1], self.freq_list1, axis=1),
                    np.take(svd_info[2], self.freq_list2, axis=1))

                ft.save_pickle(svd_info, 
                               filename_svd.replace('pair', 'pair_legendre'))

            self.pairs[pairitem].subtract_frequency_modes(
                                    svd_info[1][n_modes_start:n_modes_stop],
                                    svd_info[2][n_modes_start:n_modes_stop])

            if self.params['subtract_inputmap_from_sim'] or \
               self.params['subtract_sim_from_inputmap']:
                self.pairs_parallel_track[pairitem].subtract_frequency_modes(
                                        svd_info[1][n_modes_start:n_modes_stop],
                                        svd_info[2][n_modes_start:n_modes_stop])



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
