import shelve
import sys
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
from correlate import correlation_plots as cp
from core import handythread as ht
import multiprocessing
import subprocess

# run1:
# node: sunnyvale
# selection function: not separable, 1000 catalogs
# radio map: A cleaned with D
# threading turned off
# notes: job did not complete because of 48 hour limit
# errors: filename has _sep even though separability was not assumed
# TODO: fix filenames for this run or just throw it out
batch1_param = {
    "randnum": [0, 1, 2, 3, 4, 6, 7, 15, 16, 17, 18, 19, 20, 21, \
              42, 43, 44, 45, 59, 60, 61, 62, 63, 64, 97, 98, 99],
    "rootdir": "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/data_run1",
    "basename": "opt_x_radio_mapA",
    "randbasename": "opt_x_radio_mapArand",
    "basename_suffix": "_noconv_sep.shelve"
    }

# run2:
# node: sunnyvale
# selection function: no separable, 1000 catalogs
# radio map: A cleaned with D
# threading turned on
# notes: should be same as run1, but all 100 rand. catalogs finished
batch2_param = {
    "randnum": range(100),
    "rootdir": "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/data_run2",
    "basename": "opt_x_radio_mapA",
    "randbasename": "opt_x_radio_mapArand",
    "basename_suffix": "_noconv.shelve"
    }

# run3:
# node: sunnyvale
# selection function: separable, 1000 catalogs
# radio map: A cleaned with D
# threading turned on
batch3_param = {
    "randnum": range(100),
    "rootdir": "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/data_run3",
    "basename": "opt_x_radio_mapA",
    "randbasename": "opt_x_radio_mapArand",
    "basename_suffix": "_noconv_sep.shelve"
    }

def make_corr(filename, printcorr=False, index=None):
    """wrap the plot correlation class which reads correlation object shelve
    files"""
    output = {}
    corr_shelve = shelve.open(filename)

    corr = corr_shelve["corr"]
    run_params = corr_shelve["params"]

    print "run parameters \n" + "-" * 80
    for key in run_params:
        print key + ": " + repr(run_params[key])

    if printcorr:
        np.set_printoptions(threshold=np.nan)
        print corr_shelve["corr"]

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
    output["index"] = index

    return output

def plot_corr(shelve_entry, filename, title, coloraxis=None):
    run_params = shelve_entry["run_params"]
    print "run parameters \n" + "-" * 80
    for key in run_params:
        print key + ": " + repr(run_params[key])
    print ""

    print "binned correlation function \n" + "-" * 80
    for (lag, cdat) in zip(shelve_entry["sep_lags"], shelve_entry["pdat"]):
        print lag, cdat
    print ""

    plobj = cp.CorrelationPlot(shelve_entry["lags"], shelve_entry["corr"],
                                     shelve_entry["freq"],
                                     shelve_entry["freq_axis"])

    plobj.execute_plot_collapsed(filename, shelve_entry["sep_lags"], 
                                 shelve_entry["pdat"], shelve_entry["pdatb"],
                                 cross_power=True, title=title)

    contourname = filename.split(".")
    contourname = contourname[0]
    contourname += "_contour.png"
    plobj.plot_contour(contourname,
                        lag_inds=range(len(run_params["lags"])),
                        cross_power=True, title=title,
                        coloraxis=coloraxis)

def make_shelve_names(batch_param, random=True):
    filelist = []
    batch_rand_indices = batch_param["randnum"]
    rootdir = batch_param["rootdir"]
    basename = batch_param["basename"]
    randbasename = batch_param["randbasename"]
    basename_suffix = batch_param["basename_suffix"]

    if random:
        for runindex in batch_rand_indices:
            runindex_string = "%03d" % runindex
            shelvename = rootdir+'/'+randbasename+runindex_string+basename_suffix
            filelist.append((runindex, shelvename))
        return filelist
    else:
        return rootdir+'/'+basename+basename_suffix

def tuple_list_to_dict(list_in):
    dict_out = {}
    for listitem in list_in:
        dict_out[repr(listitem[0])] = listitem[1]
    return dict_out

def compare_corr(batchlist1, batchlist2):
    randdict1 = tuple_list_to_dict(make_shelve_names(batchlist1))
    randdict2 = tuple_list_to_dict(make_shelve_names(batchlist2))

    randindices = set(batchlist1["randnum"]).union(set(batchlist2["randnum"]))

    for randindex in randindices:
        try:
            file1 = randdict1[repr(randindex)]
        except KeyError:
            #print "dataset 1 does not have index "+repr(randindex)
            file1 = None

        try:
            file2 = randdict2[repr(randindex)]
        except KeyError:
            #print "dataset 2 does not have index "+repr(randindex)
            file2 = None

        if file1 and file2:
            print file1, file2
            corr1_shelve = shelve.open(file1)
            corr2_shelve = shelve.open(file2)
            corr1 = corr1_shelve["corr"]
            corr2 = corr2_shelve["corr"]
            run_params1 = corr1_shelve["params"]
            run_params2 = corr2_shelve["params"]
            print "-"*80
            print np.max(corr1), np.max(corr2), np.max(corr1-corr2)
            print run_params1
            print run_params2


def wrap_make_corr(runitem):
    print runitem
    (run_num, run_file) = runitem
    return make_corr(run_file, index=run_num)


def process_batch_correlations(filename, batch_param):
    product = shelve.open(filename)

    signal = make_shelve_names(batch_param, random=False)
    product["signal"] = make_corr(signal)

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    randlist = make_shelve_names(batch_param)
    results = pool.map(wrap_make_corr, randlist)
    # TODO: repack as a dictionary, wasteful ... better to write to dict
    # directly  
    for item in results:
        product["rand"+repr(item["index"])] = item

    product.close()


def wrap_plot_corr(runitem):
    print runitem
    (shelve_entry, filename, title, coloraxis) = runitem
    plot_corr(shelve_entry, filename, title, coloraxis=coloraxis)


def plot_batch_correlations(filename, batch_param):
    product = shelve.open(filename)

    signal = make_shelve_names(batch_param, random=False)
    shelve_signal = product["signal"]
    signalcorr = shelve_signal["pdat"]
    coloraxis_a = np.linspace(-0.2, 0.2, 100, endpoint=True)
    plot_corr(shelve_signal, "plots/signal_xcorr.png", "Wigglez x GBT map A",
              coloraxis=coloraxis_a)

    randlist = make_shelve_names(batch_param)
    print "number of random catalogs to stack: "+repr(len(randlist))
    rancats = np.zeros((len(randlist), 10))
    index = 0
    for (run_num, run_file) in randlist:
        print run_file
        randbasename = batch_param["randbasename"]
        shelve_entry = product["rand"+repr(run_num)]
        #plot_corr(shelve_entry, "plots/"+randbasename+".png", randbasename, coloraxis=coloraxis_a)
        rancats[index, :] = shelve_entry["pdat"] 
        ranaxis = shelve_entry["sep_lags"]
        index += 1

    # pooled plotting 
    runlist = []
    for (run_num, run_file) in randlist:
        shelve_entry = product["rand"+repr(run_num)]
        runlist.append((shelve_entry,
                        "plots/"+randbasename+repr(run_num)+".png", 
                        randbasename+repr(run_num), 
                        coloraxis_a))

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(wrap_plot_corr, runlist)

    ranstd = np.std(rancats, axis=0)
    ranmean = np.mean(rancats, axis=0)

    print "average binned correlation function and signal \n" + "-" * 80
    for (lag, cdat, cdaterr, sig) in zip(shelve_signal["sep_lags"], ranmean, ranstd, signalcorr):
        print lag, cdat, cdaterr, sig

    product.close()


#process_batch_correlations("run1_correlations.shelve", batch1_param)
#process_batch_correlations("run2_correlations.shelve", batch2_param)
#process_batch_correlations("run3_correlations.shelve", batch3_param)
#compare_corr(batch2_param, batch3_param)
#compare_corr(batch1_param, batch2_param)

#plot_batch_correlations("run1_correlations.shelve", batch1_param)
#plot_batch_correlations("run2_correlations.shelve", batch2_param)
plot_batch_correlations("run3_correlations.shelve", batch3_param)
