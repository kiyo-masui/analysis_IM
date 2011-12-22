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

# TODO: do smarter loops using fetch lists
def gather_batch_sim_run():
    datapath_db = data_paths.DataPath()
    n_runs = 100

    outpath = datapath_db.fetch("quadratic_batch_simulations")
    print "reading from to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    pwr_1d = []
    pwr_2d = []
    for index in range(n_runs):
        sim_key = "simideal_15hr_beam"
        map1_key = "db:%s:%s" % (sim_key, index)
        map2_key = "db:%s:%s" % (sim_key, index)
        noiseinv1_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"
        noiseinv2_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"
        # could also use ./observed_window.npy as the weight

        pwrspec2d_product, pwrspec1d_product = caller.execute(map1_key,
                                                              map2_key,
                                                              noiseinv1_key,
                                                              noiseinv2_key)

        pwr_2d.append(pwrspec2d_product)
        pwr_1d.append(pwrspec1d_product)

    pe.summarize_1d_agg_pwrspec(pwr_1d, "pwr_map_pair_realmask_beam_avg.dat")
    pe.summarize_2d_agg_pwrspec(pwr_2d, "pwr_map_pair_realmask_beam_avg_2d.dat")

if __name__ == '__main__':
    gather_batch_sim_run()

