r"""Code to run large batches of quadratic estimators on combinations of
data/sims
"""
import numpy as np
import math
from core import algebra
from utils import data_paths
import multiprocessing
from map import physical_gridding
from utils import binning
from correlate import map_pair as mp
from correlate import pwrspec_estimation as pe
import sys
from utils import batch_handler


def call_xspec_run(map1_key, map2_key,
                   noiseinv1_key, noiseinv2_key,
                   subtract_mean=False, degrade_resolution=False,
                   unitless=True, return_3d=False,
                   truncate=False, window=None, n_modes=None,
                   refinement=2, pad=5, order=2):
    r"""a free-standing function which calls the xspec analysis
    batteries included
    """

    # define the bad frequency list
    cutlist = [6, 7, 8, 15, 16, 18, 19, 20, 21, 22, 37, 103, 104, 105, 106,
               107, 108, 130, 131, 132, 133, 134, 237, 244, 254, 255]

    # visual inspection for wacky data #1
    augmented = [177, 194, 195, 196, 197, 198, 201, 204, 209, 213, 229]
    for entry in augmented:
        cutlist.append(entry)

    # visual inspection for wacky data #2
    augmented = [80, 171, 175, 179, 182, 183, 187, 212, 218, 219]
    for entry in augmented:
        cutlist.append(entry)

    # visual inspection of weights
    badweights = [133, 189, 192, 193, 194, 195, 196, 197, 198, 208, 209, 213,
                  233]
    for entry in badweights:
        cutlist.append(entry)

    freq = range(256)
    for entry in cutlist:
        try:
            freq.remove(entry)
        except ValueError:
            print "can't cut %d" % entry

    print "%d of %d freq slices removed" % (len(cutlist), 256)
    print freq

    # define the 1D spectral bins
    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    # initialize and calculate the xspec
    simpair = mp.MapPair(map1_key, map2_key,
                         noiseinv1_key, noiseinv2_key,
                         freq)

    if subtract_mean:
        simpair.subtract_weighted_mean()

    if degrade_resolution:
        simpair.degrade_resolution()

    if n_modes is not None:
        print "mode subtraction not implemented yet"

    retval = simpair.pwrspec_summary(window=window, unitless=unitless,
                                     bins=bins, truncate=truncate,
                                     refinement=refinement, pad=pad,
                                     order=order, return_3d=return_3d)
    return retval


def call_phys_space_run(cube1_file, cube2_file,
                       unitless=True, return_3d=False,
                       truncate=False, window="blackman"):
    """Directly call the power spectral estimation on some physical vol"""
    # define the 1D spectral bins
    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    retval = pe.calculate_xspec_file(cube1_file, cube2_file, bins,
                    weight1_file=None, weight2_file=None,
                    truncate=truncate, window=window,
                    return_3d=return_3d, unitless=unitless)

    return retval


def batch_physical_sim_run(sim_key, unitless=True, return_3d=False,
                           truncate=False, window="blackman"):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_simulations")
    print "writing to: " + outpath
    fileset = datapath_db.fetch(sim_key, intend_read=True)

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for index in fileset[0]:
        map1_file = fileset[1][index]
        map2_file = fileset[1][index]

        caller.execute(map1_file, map2_file,
                       unitless=unitless,
                       return_3d=return_3d,
                       truncate=truncate,
                       window=window)

    caller.multiprocess_stack()


def batch_sim_run(sim_key, subtract_mean=False, degrade_resolution=False,
                  unitless=True, return_3d=False,
                  truncate=False, window=None, n_modes=None,
                  refinement=2, pad=5, order=2):
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_simulations")
    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)
    for index in range(100):
        map1_key = "db:%s:%s" % (sim_key, index)
        map2_key = "db:%s:%s" % (sim_key, index)
        noiseinv1_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"
        noiseinv2_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"
        # could also use ./observed_window.npy as the weight

        caller.execute(map1_key, map2_key,
                       noiseinv1_key, noiseinv2_key,
                       subtract_mean=subtract_mean,
                       degrade_resolution=degrade_resolution,
                       unitless=unitless,
                       return_3d=return_3d,
                       truncate=truncate,
                       window=window, n_modes=n_modes,
                       refinement=refinement,
                       pad=pad, order=order)

    caller.multiprocess_stack()


def batch_wigglez_automock_run(mock_key, sel_key, subtract_mean=False,
                               degrade_resolution=False,
                               unitless=True, return_3d=False,
                               truncate=False, window=None, n_modes=None,
                               refinement=2, pad=5, order=2):

    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath
    fileset = datapath_db.fetch(mock_key, silent=True, intend_read=True)

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for index in fileset[0]:
        map1_key = "db:%s:%s" % (mock_key, index)
        map2_key = "db:%s:%s" % (mock_key, index)
        noiseinv1_key = "db:%s" % sel_key
        noiseinv2_key = "db:%s" % sel_key

        caller.execute(map1_key, map2_key,
                       noiseinv1_key, noiseinv2_key,
                       subtract_mean=subtract_mean,
                       degrade_resolution=degrade_resolution,
                       unitless=unitless,
                       return_3d=return_3d,
                       truncate=truncate,
                       window=window, n_modes=n_modes,
                       refinement=refinement,
                       pad=pad, order=order)

    caller.multiprocess_stack()


def batch_data_run(subtract_mean=False, degrade_resolution=False,
                   unitless=True, return_3d=False,
                   truncate=False, window=None, n_modes=None,
                   refinement=2, pad=5, order=2, sim=False):
    datapath_db = data_paths.DataPath()

    if sim:
        mapsim = "sim_15hr"
        outpath = datapath_db.fetch("quadratic_batch_simulations")
    else:
        mapsim = "GBT_15hr_map"
        outpath = datapath_db.fetch("quadratic_batch_data")

    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for mode_num in range(0, 55, 5):
        map1_key = "%s_cleaned_%dmode" % (mapsim, mode_num)
        map2_key = "%s_cleaned_%dmode" % (mapsim, mode_num)
        noise1_key = "%s_cleaned_%dmode" % (mapsim, mode_num)
        noise2_key = "%s_cleaned_%dmode" % (mapsim, mode_num)

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
                              verbose=True)

        print pairlist, pairdict

        for item in pairdict.keys():
            pairrun = pairdict[item]
            print pairrun['map1']
            print pairrun['noise_inv1']
            print pairrun['map2']
            print pairrun['noise_inv2']
            print pairdict[item].keys()

            caller.execute(pairrun['map1'], pairrun['map2'],
                           pairrun['noise_inv1'], pairrun['noise_inv2'],
                           subtract_mean=subtract_mean,
                           degrade_resolution=degrade_resolution,
                           unitless=unitless,
                           return_3d=return_3d,
                           truncate=truncate,
                           window=window, n_modes=n_modes,
                           refinement=refinement,
                           pad=pad, order=order)

    caller.multiprocess_stack()


def batch_GBTxwigglez_data_run(subtract_mean=False, degrade_resolution=False,
                               unitless=True, return_3d=False,
                               truncate=False, window=None, n_modes=None,
                               refinement=2, pad=5, order=2):
    datapath_db = data_paths.DataPath()

    mapkey = "GBT_15hr_map"
    outpath = datapath_db.fetch("quadratic_batch_data")

    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for mode_num in range(0, 55, 5):
        map1_key = "%s_cleaned_%dmode" % (mapkey, mode_num)
        map2_key = "%s_cleaned_%dmode" % (mapkey, mode_num)
        noise1_key = "%s_cleaned_%dmode" % (mapkey, mode_num)
        noise2_key = "%s_cleaned_%dmode" % (mapkey, mode_num)

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
                              verbose=True)

        print pairlist, pairdict

        for item in pairdict.keys():
            pairrun = pairdict[item]
            print pairrun['map1']
            print pairrun['noise_inv1']
            print pairdict[item].keys()

            caller.execute(pairrun['map1'], "db:WiggleZ_15hr_binned_data",
                           pairrun['noise_inv1'], "db:WiggleZ_15hr_montecarlo",
                           subtract_mean=subtract_mean,
                           degrade_resolution=degrade_resolution,
                           unitless=unitless,
                           return_3d=return_3d,
                           truncate=truncate,
                           window=window, n_modes=n_modes,
                           refinement=refinement,
                           pad=pad, order=order)

    caller.multiprocess_stack()


if __name__ == '__main__':

    #batch_physical_sim_run("simideal_15hr_physical")
    #batch_physical_sim_run("sim_15hr_physical")

    # real data
    #batch_data_run()
    # modeloss sim
    #batch_data_run(sim=True)
    # ideal simulations without beam
    batch_sim_run("simideal_15hr")
    sys.exit()
    # vv + evo sims without beam
    batch_sim_run("sim_15hr")
    # vv + evo sims with beam
    batch_sim_run("sim_15hr_beam")
    # vv + evo sims with treatments
    batch_sim_run("sim_15hr_beam", degrade_resolution=True, subtract_mean=True)

    batch_wigglez_automock_run("WiggleZ_15hr_mock", "WiggleZ_15hr_montecarlo")
    batch_wigglez_automock_run("WiggleZ_15hr_complete_mock", "WiggleZ_15hr_complete_separable_selection")
    batch_GBTxwigglez_data_run()
