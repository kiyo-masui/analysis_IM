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
        (bin_left, bin_center, bin_right, counts_histo, binavg, pwrspec_3d) = spec_output
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

    pwrspec_2d_p, pwrspec_2d_p, pwrspec_1d_p = \
                    pe.calculate_xspec_file(pfilename, pfilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)

    pwrspec_2d_o, pwrspec_2d_o, pwrspec_1d_o = \
                    pe.calculate_xspec_file(ofilename, ofilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)

    pwrspec_2d_ob, pwrspec_2d_ob, pwrspec_1d_ob = \
                    pe.calculate_xspec_file(obfilename, obfilename, bins,
                                    weight1_file=owindowfile,
                                    weight2_file=owindowfile,
                                    window=None, unitless=unitless)

    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    # TODO: update this all to the new format
    #for specdata in zip(bin_left, bin_center,
    #                    bin_right, pcounts_histo, obinavg, obbinavg, pbinavg,
    #                    pwrspec_input):
    #    print ("%10.15g " * 8) % specdata


if __name__ == '__main__':
    #test_with_agg_simulation(unitless=True)
    #test_with_map_pair(unitless=True)
    #test_with_simulation(unitless=True)
