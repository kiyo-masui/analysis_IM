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


def gather_batch_sim_run(sim_key, tag, subtract_mean=False,
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

    (simlist, simlistfiles) = datapath_db.fetch(sim_key, silent=True,
                                                intend_read=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in simlist:
        map1_key = "db:%s:%s" % (sim_key, index)
        map2_key = "db:%s:%s" % (sim_key, index)
        noiseinv1_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"
        noiseinv2_key = "db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"

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


def gather_physical_sim_run(sim_key, tag, unitless=True, return_3d=False,
                            truncate=False, window="blackman",
                            outdir="./plot_data", transfer=None,
                            redshift_filekey='simideal_15hr'):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()

    batchquad_path = datapath_db.fetch("quadratic_batch_simulations")
    print "reading from: " + batchquad_path
    fileset = datapath_db.fetch(sim_key, silent=True, intend_read=True)

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    caller = batch_handler.MemoizeBatch(funcname, batchquad_path,
                                        generate=False, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in fileset[0]:
        map1_file = fileset[1][index]
        map2_file = fileset[1][index]

        pwr2d_run, pwr1d_run = caller.execute(map1_file, map2_file,
                                               unitless=unitless,
                                               return_3d=return_3d,
                                               truncate=truncate,
                                               window=window)

        pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run, transfer=transfer))
        pwr_2d.append(pwr2d_run)
        pwr_1d.append(pwr1d_run)

    pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, tag,
                         outdir="./plot_data")

    # now report the theory input curve over the same bins
    bin_left = pwr_1d[0]['bin_left']
    bin_center = pwr_1d[0]['bin_center']
    bin_right = pwr_1d[0]['bin_right']

    zfilename = datapath_db.fetch(redshift_filekey, intend_read=True,
                                 pick='1')

    fileout = outdir + "/" + "theory_power_spectrum.dat"

    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    outfile = open(fileout, "w")
    for specdata in zip(bin_left, bin_center, bin_right, pwrspec_input):
        outfile.write(("%10.15g " * 4 + "\n") % specdata)

    outfile.close()

    return (pwr_1d, pwr_1d_from_2d, pwr_2d)


def plot_all_simpwr(fileroot, basics_only=False):
    r"""currently fileroot is sim, simvel or simideal"""

    gather_physical_sim_run("%s_physical" % fileroot,
                            "%s_physical" % fileroot)

    gather_batch_sim_run(fileroot,
                         fileroot)

    if not basics_only:
        gather_batch_sim_run("%s_beam" % fileroot,
                             "%s_beam" % fileroot)

        gather_batch_sim_run("%s_beam" % fileroot,
                             "%s_beam_meanconv" % fileroot,
                             degrade_resolution=True, subtract_mean=True)

        #gather_batch_sim_run("%s_beam" % fileroot,
        #                     "%s_beam_mean" % fileroot,
        #                     degrade_resolution=False, subtract_mean=True)

        #gather_batch_sim_run("%s_beam" % fileroot,
        #                     "%s_beam_conv" % fileroot,
        #                     degrade_resolution=True, subtract_mean=False)


if __name__ == '__main__':
    plot_all_simpwr("sim_15hr")
    plot_all_simpwr("simvel_15hr")
    plot_all_simpwr("simideal_15hr", basics_only=True)
    sys.exit()


