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
            pwr_1d_zero_mode = copy.deepcopy(pwr_1d)
            pwr_2d_zero_mode = copy.deepcopy(pwr_2d)

        ttag = tag + "_%dmodes_2dtrans" % mode_num
        trans2d_mode = calculate_2d_transfer_function(
                                            pwr_2d, pwr_2d_zero_mode,
                                            ttag)

        ttag = tag + "_%dmodes_1dtrans" % mode_num
        trans1d_mode = calculate_1d_transfer_function(
                                            pwr_1d, pwr_1d_zero_mode,
                                            ttag)

        transfer_functions[mode_num] = (trans1d_mode, trans2d_mode)

        mtag = tag + "_%dmodes" % mode_num

        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, mtag,
                          outdir="./plot_data")

    return transfer_functions


def calculate_2d_transfer_function(pwr_stack1, pwr_stack2, tag, outdir="./plot_data"):
    n_runs = len(pwr_stack1)
    if n_runs != len(pwr_stack2):
        print "These runs are incompatible (different number)."
        sys.exit()

    (mean1_2d, std1_2d) = pe.agg_stat_2d_pwrspec(pwr_stack1)
    (mean2_2d, std2_2d) = pe.agg_stat_2d_pwrspec(pwr_stack2)

    trans_stack = []
    for index in range(n_runs):
        entry = copy.deepcopy(pwr_stack1[index])  # TODO: copy needed?

        binavg1 = pwr_stack1[index]['binavg']
        binavg2 = pwr_stack2[index]['binavg']
        entry['binavg'] = binavg1/binavg2
        trans_stack.append(entry)

    fileout = outdir + "/" + tag + ".dat"
    trans_mean, trans_std = pe.summarize_2d_agg_pwrspec(trans_stack, fileout)

    # now derive the alternative transfer function by taking the ratio of
    # averages
    trans_alt = mean1_2d/mean2_2d
    trans_alt[np.isnan(trans_alt)] = 0.

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


def find_beam_trans(meanconv=True):
    # derive the beam, preprocessing transfer function:
    sim_15hr = gather_batch_sim_run("sim_15hr", "sim_15hr")

    sim_15hr_beam = gather_batch_sim_run("sim_15hr_beam", "sim_15hr_beam")

    sim_15hr_beam_mean = gather_batch_sim_run("sim_15hr_beam",
                                              "sim_15hr_beam_mean",
                                              degrade_resolution=False,
                                              subtract_mean=True)

    sim_15hr_beam_meanconv = gather_batch_sim_run("sim_15hr_beam",
                                                  "sim_15hr_beam_meanconv",
                                                  degrade_resolution=True,
                                                  subtract_mean=True)

    trans_beam = calculate_2d_transfer_function(sim_15hr_beam[2],
                                                sim_15hr[2],
                                                "trans_beam")

    trans_beam_mean = calculate_2d_transfer_function(
                                                sim_15hr_beam_mean[2],
                                                sim_15hr[2],
                                                "trans_beam_mean")

    trans_beam_meanconv = calculate_2d_transfer_function(
                                                sim_15hr_beam_meanconv[2],
                                                sim_15hr[2],
                                                "trans_beam_meanconv")

    # check that it makes sense by reconstructing the simulated power without
    # these effects:
    sim_15hr_corr = gather_batch_sim_run("sim_15hr_beam",
                                         "sim_15hr_beam_meanconv_corrected",
                                         degrade_resolution=True,
                                         subtract_mean=True,
                                         transfer=trans_beam_meanconv)

    return (trans_beam, trans_beam_mean, trans_beam_meanconv)


def check_pipeline():
    # make plots of the ideal cases in the physical volume, with
    # simplified P(k), etc., as an overall check.
    simideal_15hr = gather_batch_sim_run("simideal_15hr", "simideal_15hr")
    simideal_phys = gather_physical_sim_run("simideal_15hr_physical", "simideal_15hr_physical")
    sim_phys = gather_physical_sim_run("sim_15hr_physical", "sim_15hr_physical")


def find_modeloss_transfer(alt=""):
    # now gather the mode loss transfer functions
    return gather_batch_datasim_run("GBT15hr_sim", alt=alt)


if __name__ == '__main__':
    check_pipeline()
    find_beam_trans()
    find_modeloss_transfer()
