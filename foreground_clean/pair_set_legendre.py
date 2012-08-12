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
               'simfile': None,
               'sim_multiplier': 1.,
               'subtract_inputmap_from_sim': False,
               'subtract_sim_from_inputmap': False,
               'freq_list': (),
               'freq_n_all': 256,
                # in deg: (unused)
               'tack_on': None,
               'convolve': True,
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


#def wrap_find_weight(filename, regenerate=False):
#    if regenerate:
#        retval = find_weight(filename)
#    else:
#        retval = memoize_find_weight(filename)
#
#    return batch_handler.repackage_kiyo(retval)
#
#
#@batch_handler.memoize_persistent
#def memoize_find_weight(filename):
#    return find_weight(filename)


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

        self.freq_list = sp.array(self.params['freq_list'], dtype=int)
        #self.output_root = self.datapath_db.fetch(self.params['output_root'],
        #                                          tack_on=self.params['tack_on'])
        self.output_root = self.params['output_root']
        if not os.path.isdir(self.output_root):
            os.makedirs(self.output_root)
            #os.mkdir(self.output_root)

        if self.params['SVD_root']:
            if os.path.exists(self.params['SVD_root']):
                self.SVD_root = self.params['SVD_root']
            else:
                self.SVD_root = self.datapath_db.fetch(self.params['SVD_root'],
                                                   intend_write=True)
            print "WARNING: using %s to clean (intended?)" % self.SVD_root
        else:
            self.SVD_root = self.output_root

        # Write parameter file.
        parse_ini.write_params(self.params, self.output_root + 'params.ini',
                               prefix=prefix)


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
                freq_n_all = self.params['freq_n_all']
                mode_n = len(svd_info[1])
                svd_info_all = np.zeros(shape=(2, mode_n, freq_n_all))
                for i in range(mode_n):
                    np.put(svd_info_all[0][i], self.freq_list, svd_info[1][i])
                    np.put(svd_info_all[1][i], self.freq_list, svd_info[2][i])
                print 'replace mode from %d to the end' %good_modes

                svd_info = (svd_info[0],
                    rpmode.replace_modes(svd_info_all[0], good_modes, m=n_modes_stop,
                                         weight=svd_info_all[0][0]),
                    rpmode.replace_modes(svd_info_all[1], good_modes, m=n_modes_stop,
                                         weight=svd_info_all[1][0]))
                #svd_info = (svd_info[0],
                #    rpmode.replace_modes(svd_info_all[0], good_modes, m=n_modes_stop),
                #    rpmode.replace_modes(svd_info_all[1], good_modes, m=n_modes_stop))
                
                svd_info = (svd_info[0],
                    np.take(svd_info[1], self.freq_list, axis=1),
                    np.take(svd_info[2], self.freq_list, axis=1))

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
