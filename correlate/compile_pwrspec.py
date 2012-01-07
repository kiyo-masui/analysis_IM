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
from correlate import compile_wigglez_transfer as cwt
from correlate import compile_transfer as ct
import sys
from utils import batch_handler
import copy


def gather_batch_data_run(tag, subtract_mean=False, degrade_resolution=False,
                  unitless=True, return_3d=False,
                  truncate=False, window=None, n_modes=None,
                  refinement=2, pad=5, order=2, beam_transfer=None,
                  outdir="./plot_data", mode_transfer_1d=None,
                  mode_transfer_2d=None, alt=""):
    datapath_db = data_paths.DataPath()
    outpath = datapath_db.fetch("quadratic_batch_data")

    print "reading from:" + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    for mode_num in range(0,55,5):
        maptag = "GBT_15hr_map"
        map1_key = "%s_cleaned_%s%dmode" % (maptag, alt, mode_num)
        map2_key = "%s_cleaned_%s%dmode" % (maptag, alt, mode_num)
        noise1_key = "%s_cleaned_%s%dmode" % (maptag, alt, mode_num)
        noise2_key = "%s_cleaned_%s%dmode" % (maptag, alt, mode_num)

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

        transfer_2d = None

        if (mode_transfer_2d is not None) and (beam_transfer is None):
            transfer_2d = mode_transfer_2d[mode_num][1]

        if (mode_transfer_2d is None) and (beam_transfer is not None):
            transfer_2d = beam_transfer

        if (mode_transfer_2d is not None) and (beam_transfer is not None):
            transfer_2d = mode_transfer_2d[mode_num][1] * beam_transfer

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
                                  transfer=transfer_2d))

            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)

        mtag = tag + "_%dmodes" % mode_num
        if mode_transfer_1d is not None:
            transfunc = mode_transfer_1d[mode_num][0]
        else:
            transfunc = None

        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, mtag,
                          outdir="./plot_data",
                          apply_1d_transfer=transfunc)


def gather_batch_gbtxwigglez_data_run(tag, gbt_map_key, wigglez_map_key,
                                      wigglez_mock_key, wigglez_selection_key,
                                      subtract_mean=False,
                                      degrade_resolution=False,
                                      unitless=True, return_3d=False,
                                      truncate=False, window=None, n_modes=None,
                                      refinement=2, pad=5, order=2, beam_transfer=None,
                                      outdir="./plot_data", mode_transfer_1d=None,
                                      mode_transfer_2d=None):
    datapath_db = data_paths.DataPath()
    outpath = datapath_db.fetch("quadratic_batch_data")

    print "reading from:" + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    transfer_2d = None

    if (mode_transfer_2d is not None) and (beam_transfer is None):
        transfer_2d = mode_transfer_2d[mode_num][1]

    if (mode_transfer_2d is None) and (beam_transfer is not None):
        transfer_2d = beam_transfer

    if (mode_transfer_2d is not None) and (beam_transfer is not None):
        transfer_2d = mode_transfer_2d[mode_num][1] * beam_transfer

    for mode_num in range(0, 55, 5):
        print mode_num
        pwr_1d = []
        pwr_2d = []
        pwr_1d_from_2d = []
        for index in range(100):
            map1_key = "db:%s_%dmode_map" % (gbt_map_key, mode_num)
            map2_key = "db:%s:%d" % (wigglez_mock_key, index)
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

            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                  transfer=transfer_2d))

            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)

        mtag = tag + "_%dmodes_mock" % mode_num
        if mode_transfer_1d is not None:
            transfunc = mode_transfer_1d[mode_num][0]
        else:
            transfunc = None

        mean1dmock, std1dmock = pe.summarize_pwrspec(pwr_1d,
                          pwr_1d_from_2d, pwr_2d, mtag,
                          outdir="./plot_data",
                          apply_1d_transfer=transfunc)

        # now recover the xspec with the real data
        map1_key = "db:%s_%dmode_map" % (gbt_map_key, mode_num)
        map2_key = "db:%s" % wigglez_map_key
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

        pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run, transfer=transfer_2d)

        pwr_1d = pwr1d_run['binavg']
        pwr_1d_from_2d = pwr1d_from_2d['binavg']
        if mode_transfer_1d is not None:
            pwr_1d /= mode_transfer_1d[mode_num][0]
            pwr_1d_from_2d /= mode_transfer_1d[mode_num][0]

        # assume that they all have the same binning
        bin_left = pwr1d_run['bin_left']
        bin_center = pwr1d_run['bin_center']
        bin_right = pwr1d_run['bin_right']
        counts_histo = pwr1d_run['counts_histo']

        filename = "./plot_data/%s_%dmodes.dat" % (tag, mode_num)
        outfile = open(filename, "w")
        for specdata in zip(bin_left, bin_center,
                            bin_right, counts_histo, pwr_1d, pwr_1d_from_2d,
                            mean1dmock, std1dmock):
            outfile.write(("%10.15g " * 8 + "\n") % specdata)
        outfile.close()


def process_wigglez():
    (crosstrans_beam, crosstrans_beam_mean, crosstrans_beam_meanconv) = cwt.find_crossbeam_trans()

    #wigglez_transfer_functions_noconv = cwt.gather_batch_gbtxwigglez_trans_run(
    #                                "sim_15hr_combined_cleaned_noconv",
    #                                "sim_15hr",
    #                                "GBT_15hr_map_combined_cleaned_noconv",
    #                                "WiggleZ_15hr_montecarlo")

    #gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_beammodecomp_noconv",
    #                           "GBT_15hr_map_combined_cleaned_noconv",
    #                           "WiggleZ_15hr_delta_binned_data",
    #                           "WiggleZ_15hr_delta_mock",
    #                           "WiggleZ_15hr_montecarlo",
    #                           mode_transfer_1d = wigglez_transfer_functions_noconv)

    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_noconv",
                               "GBT_15hr_map_combined_cleaned_noconv",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")

    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_beamcomp_noconv",
                               "GBT_15hr_map_combined_cleaned_noconv",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo",
                               beam_transfer=crosstrans_beam)

    sys.exit()
    # TODO: this should not be noconv
    wigglez_transfer_functions = cwt.gather_batch_gbtxwigglez_trans_run(
                                    "sim_15hr_combined_cleaned_noconv",
                                    "sim_15hr",
                                    "GBT_15hr_map_combined_cleaned",
                                    "WiggleZ_15hr_montecarlo")

    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT_beammodecomp",
                               "GBT_15hr_map_combined_cleaned",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo",
                               mode_transfer_1d = wigglez_transfer_functions)

    gather_batch_gbtxwigglez_data_run("wigglez_x_GBT",
                               "GBT_15hr_map_combined_cleaned",
                               "WiggleZ_15hr_delta_binned_data",
                               "WiggleZ_15hr_delta_mock",
                               "WiggleZ_15hr_montecarlo")

def process_autopower():
    # all of the beam, meansub and conv transfer functions
    (trans_beam, trans_beam_mean, trans_beam_meanconv) = ct.find_beam_trans()

    # the data with only mean subtraction
    transfer_functions_noconv = ct.gather_batch_datasim_run("GBT15hr_sim_noconv",
                                                     alt="noconv_")

    gather_batch_data_run("GBT15hr_noconv", alt="noconv_",
                          subtract_mean=True)

    gather_batch_data_run("GBT15hr_noconv_beamcomp", alt="noconv_",
                          beam_transfer=trans_beam_mean,
                          subtract_mean=True)

    gather_batch_data_run("GBT15hr_noconv_modecomp", alt="noconv_",
                          subtract_mean=True,
                          mode_transfer_1d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconv_beammodecomp", alt="noconv_",
                          subtract_mean=True,
                          beam_transfer=trans_beam_mean,
                          mode_transfer_1d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconv_2dmodecomp", alt="noconv_",
                          subtract_mean=True,
                          mode_transfer_2d=transfer_functions_noconv)

    gather_batch_data_run("GBT15hr_noconv_beam2dmodecomp", alt="noconv_",
                          subtract_mean=True,
                          beam_transfer=trans_beam_mean,
                          mode_transfer_2d=transfer_functions_noconv)

    #gather_batch_data_run("GBT15hr_nomeanconv", alt="nomeanconv_",
    #                      subtract_mean=True)

    # data with mean subtraction and degrading
    transfer_functions = ct.find_modeloss_transfer()

    gather_batch_data_run("GBT15hr")

    gather_batch_data_run("GBT15hr_beamcomp", beam_transfer=trans_beam_meanconv)

    gather_batch_data_run("GBT15hr_modecomp", mode_transfer_1d=transfer_functions)

    gather_batch_data_run("GBT15hr_beammodecomp", beam_transfer=trans_beam_meanconv,
                          mode_transfer_1d=transfer_functions)

    gather_batch_data_run("GBT15hr_2dmodecomp", mode_transfer_2d=transfer_functions)

    gather_batch_data_run("GBT15hr_beam2dmodecomp", beam_transfer=trans_beam_meanconv,
                          mode_transfer_2d=transfer_functions)

    #gather_batch_data_run("GBT15hr", transfer=trans_beam)


if __name__ == '__main__':
    #process_wigglez()
    process_autopower()
