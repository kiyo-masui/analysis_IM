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
import copy


# TODO: do smarter loops using fetch lists
def gather_batch_sim_run(sim_key, tag, subtract_mean=False,
                         degrade_resolution=False,
                         unitless=True, return_3d=False,
                         truncate=False, window=None, n_modes=None,
                         refinement=2, pad=5, order=2, transfer=None,
                         res_path=None):

    if res_path is None:
        datapath_db = data_paths.DataPath()
        outpath = datapath_db.fetch("quadratic_batch_simulations")

    print "reading from to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    n_runs = 100
    for index in range(n_runs):
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

    outpath = datapath_db.fetch("quadratic_batch_simulations")
    print "writing to: " + outpath
    fileset = datapath_db.fetch(sim_key, intend_read=True)

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=False, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    print fileset[0]
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


def calculate_2d_transfer_function(pwr_stack1, pwr_stack2, tag, outdir="./plot_data"):
    n_runs = len(pwr_stack1)
    if n_runs != len(pwr_stack2):
        print "These runs are incompatible (different number)."
        sys.exit()

    (mean1_2d, std1_2d) = pe.agg_stat_2d_pwrspec(pwr_stack1)
    (mean2_2d, std2_2d) = pe.agg_stat_2d_pwrspec(pwr_stack2)

    trans_stack = []
    for index in range(n_runs):
        entry = {}
        entry['bin_x_left'] = pwr_stack1[index]['bin_x_left']
        entry['bin_x_center'] = pwr_stack1[index]['bin_x_center']
        entry['bin_x_right'] = pwr_stack1[index]['bin_x_right']
        entry['bin_y_left'] = pwr_stack1[index]['bin_y_left']
        entry['bin_y_center'] = pwr_stack1[index]['bin_y_center']
        entry['bin_y_right'] = pwr_stack1[index]['bin_y_right']
        entry['counts_histo'] = pwr_stack1[index]['counts_histo']

        binavg1 = pwr_stack1[index]['binavg']
        binavg2 = pwr_stack2[index]['binavg']
        entry['binavg'] = binavg1/binavg2
        trans_stack.append(entry)

    fileout = outdir + "/" + tag + ".dat"
    trans_mean, trans_std = pe.summarize_2d_agg_pwrspec(trans_stack, fileout)

    # now derive the alternative transfer function by taking the ratio of
    # averages
    trans_alt = mean1_2d/mean2_2d

    bin_x_left = np.log10(entry['bin_x_left'])
    bin_y_left = np.log10(entry['bin_y_left'])
    bin_x_center = np.log10(entry['bin_x_center'])
    bin_y_center = np.log10(entry['bin_y_center'])
    bin_x_right = np.log10(entry['bin_x_right'])
    bin_y_right = np.log10(entry['bin_y_right'])
    fileout = outdir + "/" + tag + "_altavg.dat"
    outfile = open(fileout, "w")
    for xind in range(len(bin_x_center)):
        for yind in range(len(bin_y_center)):
            outstr = ("%10.15g " * 7 + "\n") % \
                    (bin_x_left[xind], bin_x_center[xind], bin_x_right[xind], \
                     bin_y_left[yind], bin_y_center[yind], bin_y_right[yind], \
                     trans_alt[xind, yind])
            outfile.write(outstr)

    outfile.close()

    return trans_mean


def calculate_1d_transfer_function(pwr_stack1, pwr_stack2, tag, outdir="./plot_data"):
    n_runs = len(pwr_stack1)
    if n_runs != len(pwr_stack2):
        print "These runs are incompatible (different number)."
        sys.exit()

    (mean1_1d, std1_1d, corrmat1_1d) = pe.agg_stat_1d_pwrspec(pwr_stack1)
    (mean2_1d, std2_1d, corrmat2_1d) = pe.agg_stat_1d_pwrspec(pwr_stack2)

    trans_stack = []
    for index in range(n_runs):
        entry = {}
        entry['bin_left'] = pwr_stack1[index]['bin_left']
        entry['bin_center'] = pwr_stack1[index]['bin_center']
        entry['bin_right'] = pwr_stack1[index]['bin_right']
        entry['counts_histo'] = pwr_stack1[index]['counts_histo']

        binavg1 = pwr_stack1[index]['binavg']
        binavg2 = pwr_stack2[index]['binavg']
        entry['binavg'] = binavg1/binavg2
        trans_stack.append(entry)

    fileout = outdir + "/" + tag + ".dat"
    trans_mean, trans_std = pe.summarize_1d_agg_pwrspec(trans_stack, fileout)

    # now derive the alternative transfer function by taking the ratio of
    # averages
    trans_alt = mean1_1d/mean2_1d

    bin_left = entry['bin_left']
    bin_center = entry['bin_center']
    bin_right = entry['bin_right']
    fileout = outdir + "/" + tag + "_altavg.dat"
    outfile = open(fileout, "w")
    for ind in range(len(bin_center)):
        outstr = ("%10.15g " * 4 + "\n") % \
                (bin_left[ind], bin_center[ind], bin_right[ind], \
                 trans_alt[ind])
        outfile.write(outstr)

    outfile.close()

    return trans_alt

if __name__ == '__main__':
    sim_15hr = gather_batch_sim_run("sim_15hr", "sim_15hr")
    sim_15hr_beam = gather_batch_sim_run("sim_15hr_beam", "sim_15hr_beam")
    sim_15hr_beam_meanconv = gather_batch_sim_run("sim_15hr_beam", "sim_15hr_beam_meanconv", degrade_resolution=True, subtract_mean=True)
    simideal_15hr = gather_batch_sim_run("simideal_15hr", "simideal_15hr")
    sim_phys = gather_physical_sim_run("sim_15hr_physical", "sim_15hr_physical")
    simideal_phys = gather_physical_sim_run("simideal_15hr_physical", "simideal_15hr_physical")
    sim_15hr_beam_meanconv_corrected = gather_batch_sim_run("sim_15hr_beam",
                                                            "sim_15hr_beam_meanconv_corrected", degrade_resolution=True, subtract_mean=True,
                                                            transfer=trans_beam_meanconv)

    trans_beam_meanconv = calculate_2d_transfer_function(sim_15hr_beam_meanconv[2], sim_15hr[2], "trans_beam_meanconv")

    trans_beam = calculate_2d_transfer_function(sim_15hr_beam[2], sim_15hr[2], "trans_beam")

