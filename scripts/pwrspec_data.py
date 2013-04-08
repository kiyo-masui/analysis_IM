r"""Code to run large batches of quadratic estimators on combinations of
data/sims
"""
import numpy as np
import math
from core import algebra
from utils import data_paths
from simulations import corr21cm
import multiprocessing
from map import physical_gridding
from utils import binning
from correlate import map_pair as mp
from correlate import pwrspec_estimation as pe
import sys
from utils import batch_handler

def gather_batch_data_run():
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "reading from to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    for mode_num in range(0,55,5):
        map1_key = "GBT_15hr_map_cleaned_%dmode" % mode_num
        map2_key = "GBT_15hr_map_cleaned_%dmode" % mode_num
        noise1_key = "GBT_15hr_map_cleaned_%dmode" % mode_num
        noise2_key = "GBT_15hr_map_cleaned_%dmode" % mode_num

        (pairlist, pairdict) = \
                data_paths.cross_maps(map1_key, map2_key,
                              noise1_key, noise2_key,
                              map_suffix=";map",
                              noise_inv_suffix=";noise_inv",
                              cross_sym="_x_",
                              pair_former="GBTauto_cross_pairs",
                              ignore=['param'],
                              tag1prefix=map1_key + "_",
                              tag2prefix=map2_key + "_",
                              verbose=False)

        pwr_1d = []
        pwr_2d = []
        for item in pairdict.keys():
            pairrun = pairdict[item]
            print pairrun['tag1']
            print pairrun['map1']
            print pairrun['noise_inv1']
            print pairdict[item].keys()

            pwr2d_run, pwr1d_run = caller.execute(pairrun['map1'],
                                                  pairrun['map2'],
                                                  pairrun['noise_inv1'],
                                                  pairrun['noise_inv2'])
            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)

        pe.summarize_1d_agg_pwrspec(pwr_1d,
                                    "pwr_data_%dmodes_1d.dat" % mode_num)
        pe.summarize_2d_agg_pwrspec(pwr_2d,
                                    "pwr_data_%dmodes_2d.dat" % mode_num)


if __name__ == '__main__':
    gather_batch_data_run()

