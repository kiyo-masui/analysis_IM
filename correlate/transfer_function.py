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
    trans_mean, trans_std, trans_cov = pe.summarize_1d_agg_pwrspec(trans_stack, fileout)

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
