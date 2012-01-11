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
from correlate import compile_simpwrspec as csp
from utils import batch_handler
import copy


def gather_batch_datasim_run(tag, subtract_mean=False,
                             degrade_resolution=False, unitless=True,
                             return_3d=False, truncate=False, window=None,
                             n_modes=None, refinement=2, pad=5, order=2,
                             outdir="./plot_data", alt=""):
    datapath_db = data_paths.DataPath()
    outpath = datapath_db.fetch("quadratic_batch_simulations")
    print "reading from: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    transfer_functions = {}
    for mode_num in range(0,55,5):
        mapsim = "sim_15hr"
        map1_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)
        map2_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)
        noise1_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)
        noise2_key = "%s_cleaned_%s%dmode" % (mapsim, alt, mode_num)

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
            print pairrun['tag1'], pairrun['map1'], pairrun['noise_inv1']

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

            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run))
            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)

        if (mode_num == 0):
            pwr_1d_from2d_zero_mode = copy.deepcopy(pwr_1d_from_2d)
            pwr_1d_zero_mode = copy.deepcopy(pwr_1d)
            pwr_2d_zero_mode = copy.deepcopy(pwr_2d)

        ttag = tag + "_%dmodes_2dtrans" % mode_num
        trans2d_mode = tf.calculate_2d_transfer_function(
                                            pwr_2d, pwr_2d_zero_mode,
                                            ttag)

        ttag = tag + "_%dmodes_1dtrans" % mode_num
        trans1d_mode = tf.calculate_1d_transfer_function(
                                            pwr_1d, pwr_1d_zero_mode,
                                            ttag)

        ttag = tag + "_%dmodes_1dtrans_from2d" % mode_num
        trans1d_mode = tf.calculate_1d_transfer_function(
                                            pwr_1d_from_2d, pwr_1d_from2d_zero_mode,
                                            ttag)

        transfer_functions[mode_num] = (trans1d_mode, trans2d_mode)

        mtag = tag + "_%dmodes" % mode_num

        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, mtag,
                          outdir="./plot_data")

    return transfer_functions


def gather_batch_one_sided_trans_run(tag,
                                     modeloss_simkey,
                                     base_simkey,
                                     modeloss_weight_root,
                                     subtract_mean=False,
                                     degrade_resolution=False,
                                     unitless=True, return_3d=False,
                                     truncate=False, window=None, n_modes=None,
                                     refinement=2, pad=5, order=2,
                                     outdir="./plot_data"):
    r"""
    relevant parameters:
        tag -- name for output files
        ...
    """
    datapath_db = data_paths.DataPath()
    outpath = datapath_db.fetch("quadratic_batch_data")

    print "reading from:" + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    map1_key = "db:%s:0" % base_simkey
    map2_key = "db:%s:0" % base_simkey
    noiseinv1_key = "db:%s_0mode_weight" % (modeloss_weight_root)
    noiseinv2_key = "db:%s_0mode_ones" % (modeloss_weight_root)

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
        map1_key = "db:%s_%dmode_map" % (modeloss_simkey, mode_num)
        map2_key = "db:%s:0" % base_simkey
        noiseinv1_key = "db:%s_%dmode_weight" % (modeloss_weight_root, mode_num)
        noiseinv2_key = "db:%s_%dmode_ones" % (modeloss_weight_root, mode_num)

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

def find_one_sided_trans():
    gather_batch_one_sided_trans_run("GBTxGBT15hr",
                                     "sim_15hr_combined_cleaned", "sim_15hr",
                                     "GBT_15hr_map_combined_cleaned")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_noconv",
                                     "sim_15hr_combined_cleaned_noconv", "sim_15hr",
                                     "GBT_15hr_map_combined_cleaned_noconv")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_beam",
                                     "sim_15hr_combined_cleaned", "sim_15hr_beam",
                                     "GBT_15hr_map_combined_cleaned")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_noconv_beam",
                                     "sim_15hr_combined_cleaned_noconv", "sim_15hr_beam",
                                     "GBT_15hr_map_combined_cleaned_noconv")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_beam_meansubconv",
                                     "sim_15hr_combined_cleaned", "sim_15hr_beam_meansubconv",
                                     "GBT_15hr_map_combined_cleaned")

    gather_batch_one_sided_trans_run("GBTxGBT15hr_noconv_beam_meansub",
                                     "sim_15hr_combined_cleaned_noconv", "sim_15hr_beam_meansub",
                                     "GBT_15hr_map_combined_cleaned_noconv")

def find_beam_trans():
    # derive the beam, preprocessing transfer function:
    sim_15hr = csp.gather_batch_sim_run("sim_15hr", "sim_15hr")

    sim_15hr_beam = csp.gather_batch_sim_run("sim_15hr_beam", "sim_15hr_beam")

    sim_15hr_beam_mean = csp.gather_batch_sim_run("sim_15hr_beam",
                                              "sim_15hr_beam_mean",
                                              degrade_resolution=False,
                                              subtract_mean=True)

    sim_15hr_beam_meanconv = csp.gather_batch_sim_run("sim_15hr_beam",
                                                  "sim_15hr_beam_meanconv",
                                                  degrade_resolution=True,
                                                  subtract_mean=True)

    trans_beam = tf.calculate_2d_transfer_function(sim_15hr_beam[2],
                                                sim_15hr[2],
                                                "trans_beam")

    trans_beam_mean = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_mean[2],
                                                sim_15hr[2],
                                                "trans_beam_mean")

    trans_beam_meanconv = tf.calculate_2d_transfer_function(
                                                sim_15hr_beam_meanconv[2],
                                                sim_15hr[2],
                                                "trans_beam_meanconv")

    # check that it makes sense by reconstructing the simulated power without
    # these effects:
    sim_15hr_corr = csp.gather_batch_sim_run("sim_15hr_beam",
                                         "sim_15hr_beam_meanconv_corrected",
                                         degrade_resolution=True,
                                         subtract_mean=True,
                                         transfer=trans_beam_meanconv)

    return (trans_beam, trans_beam_mean, trans_beam_meanconv)


#def find_modeloss_transfer(alt=""):
#    # now gather the mode loss transfer functions
#    return gather_batch_datasim_run("GBT15hr_sim", alt=alt)


if __name__ == '__main__':
    find_one_sided_trans()
    #find_crossbeam_trans()
    #check_pipeline()
    #find_beam_trans()
    #find_modeloss_transfer()

