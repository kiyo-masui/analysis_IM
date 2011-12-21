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

def wrap_xspec(param):
    """helper for multiprocessing.map; should toast"""
    (cube1_file, cube2_file, \
     weight1_file, weight2_file, \
     bins, window, unitless) = param
    print cube1_file

    return pe.calculate_xspec_file(cube1_file, cube2_file, bins,
                                weight1_file=weight1_file,
                                weight2_file=weight2_file,
                                window=window, unitless=unitless)


def test_with_agg_simulation(unitless=True, parallel=True):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()
    #filename = datapath_db.fetch('simideal_15hr_physical', intend_read=True)
    #zfilename = datapath_db.fetch('simideal_15hr', intend_read=True,
    #                             pick='1')
    filename = datapath_db.fetch('sim_15hr_physical', intend_read=True)
    zfilename = datapath_db.fetch('sim_15hr', intend_read=True,
                                 pick='1')

    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    # give no weights; just use window
    runlist = [(filename[1][index], filename[1][index], None, None, bins,
                'blackman', unitless) for index in filename[0]]

    if parallel:
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-4)
        results = pool.map(wrap_xspec, runlist)
    else:
        for runitem in runlist:
            wrap_xspec(runitem)

    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    bin_left, bin_center, bin_right = binning.bin_edges(bins, log=True)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    agg_array = np.zeros((bin_center.size, len(filename[0])))
    counter = 0
    for spec_output in results:
        (bin_left, bin_center, bin_right, counts_histo, binavg) = spec_output
        agg_array[:,counter] = binavg
        counter += 1
        #for specdata in zip(bin_left, bin_center,
        #                    bin_right, counts_histo, binavg,
        #                    pwrspec_input):
        #    print ("%10.15g " * 6) % specdata

    meanbin = np.mean(agg_array, axis=1)
    stdbin = np.std(agg_array, axis=1)

    for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, meanbin, stdbin,
                        pwrspec_input):
        print ("%10.15g " * 7) % specdata


def test_with_simulation(unitless=True):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()
    pfilename = datapath_db.fetch('simideal_15hr_physical', intend_read=True,
                                 pick='1')
    zfilename = datapath_db.fetch('simideal_15hr', intend_read=True,
                                 pick='1')
    print pfilename
    ofilename = "./physical_cube.npy"
    obfilename = "./physical_cube_beam.npy"
    pwindowfile = "./physical_window.npy"
    owindowfile = "./observed_window_physreg.npy"

    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    bin_left, bin_center, bin_right, pcounts_histo, pbinavg = \
                    pe.calculate_xspec_file(pfilename, pfilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)

    bin_left, bin_center, bin_right, ocounts_histo, obinavg = \
                    pe.calculate_xspec_file(ofilename, ofilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)


    bin_left, bin_center, bin_right, obcounts_histo, obbinavg = \
                    pe.calculate_xspec_file(obfilename, obfilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)


    #k_vec = np.logspace(math.log10(1.e-2),
    #                    math.log10(5.),
    #                    num=55, endpoint=True)
    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    for specdata in zip(bin_left, bin_center,
                        bin_right, pcounts_histo, obinavg, obbinavg, pbinavg,
                        pwrspec_input):
        print ("%10.15g " * 8) % specdata


def test_with_map_pair(sim_key, sim_index, unitless=True):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()

    cutlist = [6, 7, 8, 15, 16, 18, 19, 20, 21, 22, 37, 103, 104, 105, 106,
               107, 108, 130, 131, 132, 133, 134, 237, 244, 254, 255]

    augmented = [177, 194, 195, 196, 197, 198, 201, 204, 209, 213, 229]
    for entry in augmented:
        cutlist.append(entry)

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

    print "%d of %d freq slices removed" % (len(cutlist), len(freq))
    print freq

    # could also use ./observed_window.npy as the weight
    map1_key = "db:%s:%s" % (sim_key, sim_index)
    simpair = mp.MapPair(map1_key,
                         map1_key,
                         "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv",
                         "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv",
                         freq)

    #simpair.subtract_weighted_mean()
    simpair.degrade_resolution()

    nbins=40
    bins = np.logspace(math.log10(0.00702349679605685),
                       math.log10(2.81187396154818),
                       num=(nbins + 1), endpoint=True)

    retval = simpair.pwrspec_1D(bins=bins, unitless=unitless)
    #bin_left, bin_center, bin_right, counts_histo, binavg = \
    #                simpair.pwrspec_1D(bins=bins, unitless=unitless)

    #for specdata in zip(bin_left, bin_center,
    #                    bin_right, counts_histo, binavg):
    #    print ("%10.15g " * 5) % specdata

def test_batch_sim_run():
    caller = batch_handler.MemoizeBatch("correlate.pwrspec_simulations.test_with_map_pair", "./testdata/",
                                        generate=True, verbose=True)
    for index in range(10):
        caller.execute("simideal_15hr_beam", index, unitless=True)

    caller.multiprocess_stack()

if __name__ == '__main__':
    #test_with_agg_simulation(unitless=True)
    #test_with_map_pair(unitless=True)
    #test_with_simulation(unitless=True)
    test_batch_sim_run()
