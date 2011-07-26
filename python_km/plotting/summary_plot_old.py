"""Make summary plots of the binned correlation functions for large sets of
random catalogs """
import os
import sys
import string
import shelve
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from correlate import correlation_plots as cp
from correlate import freq_slices as fs
from core import bootstrap
import multiprocessing
# note: keep these around until the proposal is done (ERS will delete later)

def make_corr_old(filename, verbose=False, identifier=None):
    """wrap the plot correlation class which reads correlation object shelve
    files"""
    if identifier:
        print "binning the correlation function in: " + filename + \
              ".shelve" + " with id=" + identifier
    else:
        print "binning the correlation function in: " + filename + \
              ".shelve"

    output = {}
    corr_shelve = shelve.open(filename + ".shelve")

    corr = corr_shelve["corr"]
    run_params = corr_shelve["params"]

    if verbose:
        for key in run_params:
            print key + ": " + repr(run_params[key])
        #np.set_printoptions(threshold=np.nan)
        #print corr_shelve["corr"]

    corr[np.isnan(corr)] = 0.
    corr[np.isinf(corr)] = 0.

    plobj = cp.CorrelationPlot(run_params["lags"], corr,
                                     run_params["freq"],
                                     corr_shelve["freq_axis"])

    lag_inds = list(range(len(run_params["lags"])))
    freq_diffs = sp.arange(0.1e6, 100e6, 200.0 / 256 * 1e6)
    nbins = 10

    (sep_lags, pdat, pdatb) = plobj.bin_correlation_sep(lag_inds, freq_diffs,
                                                       nbins, norms=False,
                                                       cross_power=True)

    output["run_params"] = run_params
    output["lags"] = run_params["lags"]
    output["corr"] = corr
    output["freq"] = run_params["freq"]
    output["freq_axis"] = corr_shelve["freq_axis"]
    output["sep_lags"] = sep_lags
    output["pdat"] = pdat
    output["pdatb"] = pdatb
    if identifier:
        output["identifier"] = identifier

    return output


def plot_corr_old(shelve_entry, filename, title, coloraxis=None):
    """make plots for a single correlation function (internal)"""

    print "plotting " + filename
    outlog = open(filename + ".txt", 'w')
    run_params = shelve_entry["run_params"]
    outlog.write("run parameters \n" + "-" * 80 + "\n")
    for key in run_params:
        outlog.write(key + ": " + repr(run_params[key]) + "\n")
    outlog.write("\n")

    outlog.write("binned correlation function \n" + "-" * 80 + "\n")
    for (lag, cdat) in zip(shelve_entry["sep_lags"], shelve_entry["pdat"]):
        outlog.write(repr(lag) + repr(cdat) + "\n")

    plobj = cp.CorrelationPlot(shelve_entry["lags"], shelve_entry["corr"],
                                     shelve_entry["freq"],
                                     shelve_entry["freq_axis"])

    plobj.execute_plot_collapsed(filename + ".png", shelve_entry["sep_lags"],
                                 shelve_entry["pdat"], shelve_entry["pdatb"],
                                 cross_power=True, title=title)

    plobj.plot_contour(filename + "_contour.png",
                        lag_inds=range(len(run_params["lags"])),
                        cross_power=True, title=title,
                        coloraxis=coloraxis)


def batch_correlations_statistics_old(filename, batch_param, randtoken="rand"):
    """bin a large batch of correlation functions"""
    master = shelve.open(filename)
    filelist = make_shelve_names(batch_param)

    # make a sublist of calculated correlations for just the random trials
    randlist = []
    for item in filelist:
        if string.find(item[0], randtoken) != -1:
            randlist.append(item)

    print "number of random catalogs to stack: " + repr(len(randlist))
    # TODO remove 10 magic number and have it figure this out instead
    rancats = np.zeros((len(randlist), 10))
    index = 0
    for (run_id, run_file) in randlist:
        print run_file
        shelve_entry = master[run_id]
        rancats[index, :] = shelve_entry["pdat"]
        #ranaxis = shelve_entry["sep_lags"]
        index += 1

    # TODO: have it be smarter about this rather than assume "signal" exists
    shelve_signal = master["signal"]

    ranstd = np.std(rancats, axis=0)
    ranmean = np.mean(rancats, axis=0)
    print "average binned correlation function and signal \n" + "-" * 80
    output_package = zip(shelve_signal["sep_lags"], ranmean,
                                         ranstd, shelve_signal["pdat"])
    for (lag, cdat, cdaterr, sig) in output_package:
        print lag, cdat, cdaterr, sig

    master.close()
    return output_package



