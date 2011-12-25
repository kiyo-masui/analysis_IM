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
def gather_batch_sim_run(sim_key, tag, subtract_mean=False, degrade_resolution=False,
                  unitless=True, return_3d=False,
                  truncate=False, window=None, n_modes=None,
                  refinement=2, pad=5, order=2, transfer=None):
    datapath_db = data_paths.DataPath()
    n_runs = 100

    outpath = datapath_db.fetch("quadratic_batch_simulations")
    #outpath = "/mnt/raid-project/gmrt/eswitzer/quadratic_products/simulations_v1"
    print "reading from to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in range(n_runs):
        map1_key = "db:%s:%s" % (sim_key, index)
        map2_key = "db:%s:%s" % (sim_key, index)
        noiseinv1_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"
        noiseinv2_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"
        # could also use ./observed_window.npy as the weight

        pwrspec2d_product, pwrspec1d_product = caller.execute(map1_key, map2_key,
                                               noiseinv1_key, noiseinv2_key,
                                               subtract_mean=subtract_mean,
                                               degrade_resolution=degrade_resolution,
                                               unitless=unitless,
                                               return_3d=return_3d,
                                               truncate=truncate,
                                               window=window, n_modes=n_modes,
                                               refinement=refinement,
                                               pad=pad, order=order)

        pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec2d_product, transfer=transfer))

        pwr_2d.append(pwrspec2d_product)
        pwr_1d.append(pwrspec1d_product)

    fileout = tag + "_avg_from2d.dat"
    corr_fileout = tag + "_corr_from2d.dat"
    pe.summarize_1d_agg_pwrspec(pwr_1d_from_2d, fileout,
                                corr_file=corr_fileout)

    fileout = tag + "_avg.dat"
    corr_fileout = tag + "_corr.dat"
    pe.summarize_1d_agg_pwrspec(pwr_1d, fileout,
                                corr_file=corr_fileout)

    fileout = tag + "_avg_2d.dat"
    pe.summarize_2d_agg_pwrspec(pwr_2d, fileout)

    return (pwr_1d, pwr_1d_from_2d, pwr_2d)


def gather_batch_data_run(tag, subtract_mean=False, degrade_resolution=False,
                  unitless=True, return_3d=False,
                  truncate=False, window=None, n_modes=None,
                  refinement=2, pad=5, order=2, transfer=None):
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
        pwr_1d_from_2d = []
        for item in pairdict.keys():
            pairrun = pairdict[item]
            print pairrun['tag1']
            print pairrun['map1']
            print pairrun['noise_inv1']
            print pairdict[item].keys()

            pwr2d_run, pwr1d_run = caller.execute(pairrun['map1'],
                                                  pairrun['map2'],
                                                  pairrun['noise_inv1'],
                                                  pairrun['noise_inv2'],
                                                  subtract_mean=subtract_mean,
                                                  degrade_resolution=degrade_resolution,
                                                  unitless=unitless,
                                                  return_3d=return_3d,
                                                  truncate=truncate,
                                                  window=window, n_modes=n_modes,
                                                  refinement=refinement,
                                                  pad=pad, order=order)

            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                  transfer=transfer))
            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)

        fileout = tag + "_avg_data_%dmodes_from2d.dat" % mode_num
        corr_fileout = tag + "_corr_data_%dmodes_from2d.dat" % mode_num
        pe.summarize_1d_agg_pwrspec(pwr_1d_from_2d, fileout,
                                corr_file=corr_fileout)

        fileout = tag + "_avg_data_%dmodes.dat" % mode_num
        corr_fileout = tag + "_corr_data_%dmodes.dat" % mode_num
        pe.summarize_1d_agg_pwrspec(pwr_1d, fileout,
                                corr_file=corr_fileout)

        fileout = tag + "_avg_data_%dmodes_2d.dat" % mode_num
        pe.summarize_2d_agg_pwrspec(pwr_2d, fileout)


def calculate_transfer_function(pwr_stack1, pwr_stack2, tag):
    #n_runs = len(pwr_stack1)
    #if n_runs != len(pwr_stack2):
    #    print "These runs are incompatible (different number)."
    #    sys.exit()
    n_runs = 100

    trans_stack = []
    for index in range(n_runs):
        entry = {}
        entry['bin_x_left'] = pwr_stack1[2][index]['bin_x_left']
        entry['bin_x_center'] = pwr_stack1[2][index]['bin_x_center']
        entry['bin_x_right'] = pwr_stack1[2][index]['bin_x_right']
        entry['bin_y_left'] = pwr_stack1[2][index]['bin_y_left']
        entry['bin_y_center'] = pwr_stack1[2][index]['bin_y_center']
        entry['bin_y_right'] = pwr_stack1[2][index]['bin_y_right']
        entry['counts_histo'] = pwr_stack1[2][index]['counts_histo']

        binavg1 = pwr_stack1[2][index]['binavg']
        binavg2 = pwr_stack2[2][index]['binavg']
        entry['binavg'] = binavg1/binavg2
        trans_stack.append(entry)

    fileout = tag + ".dat"
    trans_mean, trans_std = pe.summarize_2d_agg_pwrspec(trans_stack, fileout)

    return trans_mean

if __name__ == '__main__':
    simideal_15hr = gather_batch_sim_run("simideal_15hr", "simideal_15hr")
    sim_15hr = gather_batch_sim_run("sim_15hr", "sim_15hr")
    sim_15hr_beam = gather_batch_sim_run("sim_15hr_beam", "sim_15hr_beam")
    sim_15hr_beam_meanconv = gather_batch_sim_run("sim_15hr_beam", "sim_15hr_beam_meanconv", degrade_resolution=True, subtract_mean=True)
    trans_beam_meanconv = calculate_transfer_function(sim_15hr_beam_meanconv, sim_15hr, "trans_beam_meanconv")

    trans_beam = calculate_transfer_function(sim_15hr_beam, sim_15hr, "trans_beam")

    sim_15hr_beam_meanconv_corrected = gather_batch_sim_run("sim_15hr_beam",
                                                            "sim_15hr_beam_meanconv_corrected", degrade_resolution=True, subtract_mean=True,
                                                            transfer=trans_beam_meanconv)

    gather_batch_data_run("GBT15hr_wbeam")
    gather_batch_data_run("GBT15hr", transfer=trans_beam_meanconv)
    #gather_batch_data_run("GBT15hr", transfer=trans_beam)
