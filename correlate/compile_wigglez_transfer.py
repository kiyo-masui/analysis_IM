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
from correlate import transfer_function as tf
import sys
from utils import batch_handler
import copy


def gather_batch_gbtxwigglez_trans_run(tag,
                                       sim_key,
                                       sim_wigglez,
                                       base_sim_GBT,
                                       gbt_map_key,
                                       wigglez_selection_key,
                                       subtract_mean=False,
                                       degrade_resolution=False,
                                       unitless=True, return_3d=False,
                                       truncate=False, window=None, n_modes=None,
                                       refinement=2, pad=5, order=2,
                                       outdir="./plot_data"):
    r"""
    relevant parameters:
        tag -- name for output files
        sim_key -- cleaned GBT map simulations
        sim_wigglez -- wiggleZ sim
        base_sim_GBT -- GBT sim to compare to in the transfer, e.g. no beam
        wigglez_selection_key -- point to the selection function
    """
    datapath_db = data_paths.DataPath()
    outpath = datapath_db.fetch("quadratic_batch_data")

    print "reading from:" + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    map1_key = "db:%s:0" % base_sim_GBT
    map2_key = "db:%s:0" % sim_wigglez
    noiseinv1_key = "db:%s_0mode_weight" % (gbt_map_key)
    noiseinv2_key = "db:%s" % wigglez_selection_key

    pwr2d_run, pwr1d_run = caller.execute(map1_key, map2_key,
                   noiseinv1_key, noiseinv2_key,
                   subtract_mean=subtract_mean,
                   degrade_resolution=degrade_resolution,
                   unitless=unitless,
                   return_3d=return_3d,
                   truncate=truncate,
                   window=window, n_modes=n_modes,
                   refinement=refinement,
                   pad=pad, order=order)

    pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run)

    pwr_1d_zero = pwr1d_run['binavg']
    pwr_2d_zero = pwr2d_run['binavg']
    pwr_1d_from_2d_zero = pwr1d_from_2d['binavg']

    transfer_functions = {}
    for mode_num in range(0, 55, 5):
        print mode_num
        map1_key = "db:%s_%dmode_map" % (sim_key, mode_num)
        map2_key = "db:%s:0" % sim_wigglez
        noiseinv1_key = "db:%s_%dmode_weight" % (gbt_map_key, mode_num)
        noiseinv2_key = "db:%s" % wigglez_selection_key

        pwr2d_run, pwr1d_run = caller.execute(map1_key, map2_key,
                       noiseinv1_key, noiseinv2_key,
                       subtract_mean=subtract_mean,
                       degrade_resolution=degrade_resolution,
                       unitless=unitless,
                       return_3d=return_3d,
                       truncate=truncate,
                       window=window, n_modes=n_modes,
                       refinement=refinement,
                       pad=pad, order=order)

        pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run)

        pwr_1d = pwr1d_run['binavg']
        pwr_1d_from_2d = pwr1d_from_2d['binavg']
        trans1d_mode = pwr_1d/pwr_1d_zero
        trans1d_from2d_mode = pwr_1d_from_2d/pwr_1d_from_2d_zero
        trans2d_mode = pwr2d_run['binavg']/pwr_2d_zero

        transfer_functions[mode_num] = (trans1d_mode,
                                        trans1d_from2d_mode,
                                        trans2d_mode)

        # assume that they all have the same binning
        bin_left = pwr1d_run['bin_left']
        bin_center = pwr1d_run['bin_center']
        bin_right = pwr1d_run['bin_right']
        counts_histo = pwr1d_run['counts_histo']

        filename = "./plot_data/%s_%dmodes_trans.dat" % (tag, mode_num)
        outfile = open(filename, "w")
        for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, pwr_1d, trans1d_mode,
                        trans1d_from2d_mode):
            outfile.write(("%10.15g " * 7 + "\n") % specdata)
        outfile.close()

    return transfer_functions


def gather_batch_genericsim_run(simleft_key, simright_key,
                         weightleft_key, weightright_key,
                         tag, subtract_mean=False,
                         degrade_resolution=False,
                         unitless=True, return_3d=False,
                         truncate=False, window=None, n_modes=None,
                         refinement=2, pad=5, order=2, transfer=None,
                         res_path=None):

    if res_path is None:
        datapath_db = data_paths.DataPath()
        batchquad_path = datapath_db.fetch("quadratic_batch_simulations")

    print "reading from: " + batchquad_path

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, batchquad_path,
                                        verbose=True)

    (simlist, simlistfiles) = datapath_db.fetch(simleft_key, silent=True,
                                                intend_read=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in simlist:
        map1_key = "db:%s:%s" % (simleft_key, index)
        map2_key = "db:%s:%s" % (simright_key, index)
        noiseinv1_key = "db:%s" % weightleft_key
        noiseinv2_key = "db:%s" % weightright_key

        pwr2d_run, pwr1d_run = caller.execute(map1_key, map2_key,
                                              noiseinv1_key, noiseinv2_key,
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

    pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, tag,
                         outdir="./plot_data")

    return (pwr_1d, pwr_1d_from_2d, pwr_2d)



def find_crossbeam_trans():
    sim_15hr = gather_batch_genericsim_run("sim_15hr", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr")

    sim_15hr_beam = gather_batch_genericsim_run("sim_15hr_beam", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr_beam")

    sim_15hr_beam_mean = gather_batch_genericsim_run("sim_15hr_beam_meansub",
                         "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo",
                         "sim_wigglezxGBT15hr_beam_meansub")

    # TODO: should not use the noconv weight here... need to separate these
    # better
    sim_15hr_beam_meanconv = gather_batch_genericsim_run("sim_15hr_beam_meansubconv",
                         "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo",
                         "sim_wigglezxGBT15hr_beam_meansubconv")

    trans_beam = tf.calculate_2d_transfer_function(sim_15hr_beam[2],
                                                sim_15hr[2],
                                                "crosstrans_beam")

    trans_beam_mean = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_mean[2],
                                                sim_15hr[2],
                                                "crosstrans_beam_mean")

    trans_beam_meanconv = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_meanconv[2],
                                                sim_15hr[2],
                                                "crosstrans_beam_meanconv")

    # TODO: check that it makes sense by reconstructing the simulated power without
    # these effects:
    #sim_15hr_corr = gather_batch_sim_run("sim_15hr_beam",
    #                                     "sim_15hr_beam_meanconv_corrected",
    #                                     degrade_resolution=True,
    #                                     subtract_mean=True,
    #                                     transfer=trans_beam_meanconv)

    return (trans_beam, trans_beam_mean, trans_beam_meanconv)


def wigglez_x_GBT():
    gather_batch_genericsim_run("simvel_15hr", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "simvel_wigglezxGBT15hr")

    gather_batch_genericsim_run("simvel_15hr_beam", "simvel_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "simvel_wigglezxGBT15hr_beam")

    gather_batch_genericsim_run("sim_15hr", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr")

    gather_batch_genericsim_run("sim_15hr_beam", "sim_15hr_delta",
                         "GBT_15hr_map_combined_cleaned_noconv_0mode_weight",
                         "WiggleZ_15hr_montecarlo", "sim_wigglezxGBT15hr_beam")


if __name__ == '__main__':
    find_crossbeam_trans()
    wigglez_x_GBT()
    sys.exit()


