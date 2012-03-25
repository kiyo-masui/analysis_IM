from utils import data_paths
from correlate import pwrspec_estimation as pe
from correlate import transfer_function as tf
from utils import batch_handler
import copy
from correlate import batch_quadratic as bq


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

        pe.summarize_agg_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, mtag,
                          outdir="./plot_data")

    return transfer_functions
