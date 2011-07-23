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
# TODO: strip out _old binning code functions
# TODO: convert batch params into ini files
# TODO: make sure all methods here used counts/weights as-saved
# TODO: new methods are fast enough that no need to partition binning plotting,
# and statistics; still use multiprocessing

# note that notes are purely human-readable and the keys do not mean anything
rootdir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/"
run1_notes = {
    "runname": "run1",
    "machine": "sunnyvale",
    "selection_function": "not separable, 1000 catalogs",
    "radio_map": "A cleaned with D",
    "threading": "off",
    "error_1": "job did not complete because of 48-hour limit",
    "error_2": "filename has _sep even though separability was not assumed",
    "error_3": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    "todo": "file filenames"
    }
batch1_param = {
    "path": rootdir + "data_run1",
    "rand:list": {"prefix": "opt_x_radio_mapArand",
                  "suffix": "_noconv_sep",
                  "indices": [0, 1, 2, 3, 4, 6, 7, 15, 16, 17, \
                              18, 19, 20, 21, 42, 43, 44, 45, 59, \
                              60, 61, 62, 63, 64, 97, 98, 99],
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_mapA_noconv_sep",
    "notes": run1_notes
    }

run2_notes = {
    "runname": "run1",
    "machine": "sunnyvale",
    "selection_function": "not separable, 1000 catalogs",
    "radio_map": "A cleaned with D",
    "threading": "on",
    "note": "should be same as run1, but all 100 rand. catalogs finished",
    "error": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    }
batch2_param = {
    "path": rootdir + "data_run2",
    "rand:list": {"prefix": "opt_x_radio_mapArand",
                  "suffix": "_noconv",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_mapA_noconv",
    "notes": run2_notes
    }

run3_notes = {
    "runname": "run3",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "A cleaned with D",
    "threading": "on",
    "error": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    }
batch3_param = {
    "path": rootdir + "data_run3",
    "rand:list": {"prefix": "opt_x_radio_mapArand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_mapA_noconv_sep",
    "notes": run3_notes
    }

run4_notes = {
    "runname": "run4",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD",
    "threading": "on",
    "error": "THESE MAPS ONLY HAD 2 MODES REMOVED",
    }
batch4_param = {
    "path": rootdir + "data_run4",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "notes": run4_notes
    }

run5_notes = {
    "runname": "run5",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "sec A of Kiyo's old 2-way split (signal seen)",
    "threading": "on",
    "note": "this run also includes xcorr with selection function",
    }
truncated = range(100)
truncated.remove(2)
batch5_param = {
    "path": rootdir + "data_run5",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": truncated,
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "selxcorr:file": "optsel_x_radio_combined_noconv_sep",
    "notes": run5_notes
    }

run6_notes = {
    "runname": "run6",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "threading": "on",
    }
batch6_param = {
    "path": rootdir + "data_run6",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "selxcorr:file": "optsel_x_radio_combined_noconv_sep",
    "notes": run6_notes
    }

run7_notes = {
    "runname": "run7",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "notes": "should match run6 exactly, but with new correlate()"
    }
batch7_param = {
    "path": rootdir + "data_run7",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    #"selxcorr:file": "optsel_x_radio_combined_noconv_sep",
    "notes": run7_notes
    }

run8_notes = {
    "runname": "mode_xcorr_run1",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 0-25 modes removed",
    "speedup": "on"
    }
batch8_param = {
    "path": rootdir + "mode_xcorr_run1",
    "mode:list": {"prefix": "opt_x_radio_combined_mode",
                  "suffix": "_noconv_sep",
                  "indices": range(26),
                  "indexfmt": "%d",
                  "id_prefix": "mode"},
    "nullmode:list": {"prefix": "opt_x_radio_combined_nullmode",
                  "suffix": "_noconv_sep",
                  "indices": range(26),
                  "indexfmt": "%d",
                  "id_prefix": "nullmode"},
    "notes": run8_notes
    }

run9_notes = {
    "runname": "test_69old_vs_69new",
    "machine": "prawn",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "new and old treatments up to session 69",
    "speedup": "on"
    }
batch9_param = {
    "path": rootdir + "test_69old_vs_69new",
    "reallyoldway:file": "opt_x_radio_69old_kiyo_noconv_sep",
    "oldway:file": "opt_x_radio_69old_noconv_sep",
    "newway:file": "opt_x_radio_69_noconv_sep",
    "notes": run9_notes
    }

run10_notes = {
    "runname": "run10",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "notes": "should match run7 but using counts as corr bin weights"
    }
batch10_param = {
    "path": rootdir + "xcorr_15modes",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_sep",
    "notes": run10_notes
    }

run11_notes = {
    "runname": "run11",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "convolution": "on",
    "notes": "same as 10, but with convolution of the optical overdensity"
    }
batch11_param = {
    "path": rootdir + "xcorr_withconv",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_conv_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_conv_sep",
    "notes": run11_notes
    }

run12_notes = {
    "runname": "run12",
    "machine": "sunnyvale",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "weighted average of ABCD, 15 modes removed",
    "speedup": "on",
    "meansubtract": "off",
    "notes": "same as 10, but without subtracting means"
    }
batch12_param = {
    "path": rootdir + "xcorr_nomean",
    "rand:list": {"prefix": "opt_x_radio_combined_rand",
                  "suffix": "_noconv_nomean_sep",
                  "indices": range(100),
                  "indexfmt": "%03d",
                  "id_prefix": "rand"},
    "signal:file": "opt_x_radio_combined_noconv_nomean_sep",
    "notes": run12_notes
    }

# /mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/xcorr_15modes_12split/opt_x_radio_15CD_noconv_sep.shelve
run13_notes = {
    "runname": "xcorr_15modes_12split",
    "machine": "prawn",
    "selection_function": "separable, 1000 catalogs",
    "radio_map": "12 ABCD maps",
    "speedup": "on"
    }
batch13_param = {
    "path": rootdir + "xcorr_15modes_12split",
    "xcorr15mAB:file": "opt_x_radio_15AB_noconv_sep",
    "xcorr15mAC:file": "opt_x_radio_15AC_noconv_sep",
    "xcorr15mAD:file": "opt_x_radio_15AD_noconv_sep",
    "xcorr15mBA:file": "opt_x_radio_15BA_noconv_sep",
    "xcorr15mBC:file": "opt_x_radio_15BC_noconv_sep",
    "xcorr15mBD:file": "opt_x_radio_15BD_noconv_sep",
    "xcorr15mCA:file": "opt_x_radio_15CA_noconv_sep",
    "xcorr15mCB:file": "opt_x_radio_15CB_noconv_sep",
    "xcorr15mCD:file": "opt_x_radio_15CD_noconv_sep",
    "xcorr15mDA:file": "opt_x_radio_15DA_noconv_sep",
    "xcorr15mDB:file": "opt_x_radio_15DB_noconv_sep",
    "xcorr15mDC:file": "opt_x_radio_15DC_noconv_sep",
    "notes": run13_notes
    }


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

def make_corr(filename, verbose=False, identifier=None):
    """wrap the plot correlation class which reads correlation object shelve
    files; uses new binning methods in freq-slices"""
    output = {}
    corr_shelve = shelve.open(filename + ".shelve")

    corr = corr_shelve["corr"]
    run_params = corr_shelve["params"]

    try:
        corr_counts = corr_shelve["counts"]
    except KeyError:
        print "WARNING: unable to find counts weighting for correlation"
        corr_counts = None

    if identifier:
        print "binning the correlation function in: " + filename + \
              ".shelve" + " with id=" + identifier
    else:
        print "binning the correlation function in: " + filename + \
              ".shelve"

    if verbose:
        for key in run_params:
            print key + ": " + repr(run_params[key])
        #np.set_printoptions(threshold=np.nan)
        #print corr_shelve["corr"]

    corr[np.isnan(corr)] = 0.
    corr[np.isinf(corr)] = 0.

    lags = sp.array(run_params["lags"])
    frange = run_params["freq"]
    realrange = corr_shelve["freq_axis"]
    corr_2D = fs.rebin_corr_freq_lag(corr, realrange[list(frange)],
                                     weights=corr_counts, return_fbins=True, nfbins=200)
    corr_1D = fs.collapse_correlation_1D(corr_2D[0], corr_2D[2], lags, weights=corr_2D[1])

    output["run_params"] = run_params
    output["lags"] = run_params["lags"]
    output["corr"] = corr
    output["corr_counts"] = corr_counts
    output["freq"] = run_params["freq"]
    output["freq_axis"] = corr_shelve["freq_axis"]
    output["corr1D"] = corr_1D[0]
    output["corr1D_weights"] = corr_1D[1]
    output["corr1D_lags"] = corr_1D[2]
    output["corr2D"] = corr_2D[0]
    output["corr2D_weights"] = corr_2D[1]
    output["corr2D_fbins"] = corr_2D[2]
    if identifier:
        output["identifier"] = identifier

    return output


def plot_corr(shelve_entry, filename, title, coloraxis=None):
    """make plots for a single correlation function (internal)"""

    print "plotting " + filename
    outlog = open(filename + ".txt", 'w')
    run_params = shelve_entry["run_params"]
    outlog.write("run parameters \n" + "-" * 80 + "\n")
    for key in run_params:
        outlog.write(key + ": " + repr(run_params[key]) + "\n")
    outlog.write("\n")

    outlog.write("binned correlation function \n" + "-" * 80 + "\n")
    for (lag, cdat) in zip(shelve_entry["corr1D_lags"], shelve_entry["corr1D"]):
        outlog.write(repr(lag) + repr(cdat) + "\n")

    plot_collapsed(filename + ".png", shelve_entry["corr1D_lags"],
                                 shelve_entry["corr1D"],
                                 cross_power=True, title=title)

    plot_contour(filename + "_contour.png", shelve_entry["corr2D_fbins"],
                 shelve_entry["lags"], shelve_entry["corr2D"],
                 cross_power=True, title=title, coloraxis=coloraxis)


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


def make_shelve_names(batch_param):
    """assemble a list of shelve file names containing the correlation function
    information
    """
    filelist = []
    for item in batch_param:
        token = item.split(":")
        if (len(token) == 2):
            entry_type = token[1]
            ident = token[0]
            if (entry_type == "file"):
                fullpath = batch_param["path"] + "/" + batch_param[item]
                filelist.append((ident, fullpath))
            if (entry_type == "list"):
                list_param = batch_param[item]
                for index in list_param["indices"]:
                    indexstr = list_param["indexfmt"] % index
                    fullpath = batch_param["path"] + "/"
                    fullpath += list_param["prefix"] + indexstr
                    fullpath += list_param["suffix"]
                    listid = list_param["id_prefix"] + indexstr
                    filelist.append((listid, fullpath))

    return filelist


def tuple_list_to_dict(list_in):
    """Simple internal method to convert a list of tuples to a dictionary"""
    dict_out = {}
    for listitem in list_in:
        dict_out[repr(listitem[0])] = listitem[1]
    return dict_out


def compare_corr_one(file_a, file_b, print_params=False):
    corr_a_shelve = shelve.open(file_a)
    corr_b_shelve = shelve.open(file_b)
    corr_a = corr_a_shelve["corr"]
    corr_b = corr_b_shelve["corr"]

    delta = np.max(corr_a - corr_b)
    #print corr_a - corr_b
    difference = False
    if (delta != 0.):
        difference = True
        print "-" * 80
        print "comparing " + file_a, file_b
        print np.max(corr_a), np.max(corr_b), np.max(corr_a - corr_b)
        if print_params:
            run_params_a = corr_a_shelve["params"]
            run_params_b = corr_b_shelve["params"]
            print run_params_a
            print run_params_b
    else:
        print "comparing " + file_a, file_b

    return difference


def compare_corr(batchlist_a, batchlist_b, print_params=False):
    """Compare two sets of correlation function runs"""
    filedict_a = tuple_list_to_dict(make_shelve_names(batchlist_a))
    filedict_b = tuple_list_to_dict(make_shelve_names(batchlist_b))

    common_ids = set(filedict_a.keys()).union(set(filedict_b.keys()))

    difflist = []
    for ident in common_ids:
        try:
            file_a = filedict_a[ident]
        except KeyError:
            #print "dataset 1 does not have index " + ident
            file_a = None

        try:
            file_b = filedict_b[ident]
        except KeyError:
            #print "dataset 2 does not have index " + ident
            file_b = None

        if file_a and file_b:
            difference = compare_corr_one(file_a + ".shelve", 
                                          file_b + ".shelve",
                                          print_params=print_params)
            difflist.append((difference, file_a, file_b))

    return difflist


def wrap_make_corr(runitem):
    """wrapper to the make correlation function for the process pool"""
    (run_id, run_file) = runitem
    return make_corr(run_file, identifier=run_id)


def process_batch_correlations(filename, batch_param):
    """Process a batch of correlation functions"""
    product = shelve.open(filename)
    filelist = make_shelve_names(batch_param)

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = pool.map(wrap_make_corr, filelist)
    # TODO: can pool write this dictionary instead?
    for item in results:
        product[item["identifier"]] = item

    product.close()


def wrap_plot_corr(runitem):
    """wrapper to the correlation function plotter for the process pool"""
    (shelve_entry, filename, title, coloraxis) = runitem
    plot_corr(shelve_entry, filename, title, coloraxis=coloraxis)


def batch_correlations_statistics(filename, batch_param, randtoken="rand"):
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
        rancats[index, :] = shelve_entry["corr1D"]
        #ranaxis = shelve_entry["sep_lags"]
        index += 1

    # TODO: have it be smarter about this rather than assume "signal" exists
    shelve_signal = master["signal"]
    master.close()

    ranstd = np.std(rancats, axis=0)
    ranmean = np.mean(rancats, axis=0)
    rancov = np.corrcoef(rancats, rowvar=0)
    plot_covariance(rancov, "bin-bin_cov.png",
                    axis_labels = shelve_entry["corr1D_lags"])
    print "average binned correlation function and signal \n" + "-" * 80
    output_package = zip(shelve_signal["corr1D_lags"], ranmean,
                                         ranstd, shelve_signal["corr1D"])
    for (lag, cdat, cdaterr, sig) in output_package:
        print lag, cdat, cdaterr, sig

    return output_package


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


def plot_batch_correlations(filename, batch_param, dir_prefix="plots/",
                            color_range=[-0.2, 0.2]):
    """bin a large batch of correlation functions"""
    master = shelve.open(filename)
    filelist = make_shelve_names(batch_param)
    coloraxis = np.linspace(color_range[0], color_range[1], 100, endpoint=True)

    d = os.path.dirname(dir_prefix)
    if not os.path.exists(d):
        print "making the directory: " + dir_prefix
        os.makedirs(d)

    # write out the notes if they exist
    try:
        outlog = open(dir_prefix + "notes.txt", 'w')
        notes = batch_param["notes"]
        for key in notes:
            outlog.write(key + ": " + notes[key] + "\n")
        outlog.write("\n")
    except KeyError:
        print "This batch did not have any notes"

    # pooled plotting
    # TODO: save a .txt file with the run params along with each plot
    runlist = []
    for (run_id, run_file) in filelist:
        shelve_entry = master[run_id]
        run_tag = run_file.split('/')
        run_tag = run_tag[-1]
        runlist.append((shelve_entry, dir_prefix + run_tag,
                        run_tag, coloraxis))

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(wrap_plot_corr, runlist)

    master.close()


# TODO: label lag axes better 
def plot_covariance(matrix_in, filename, axis_labels = [], mask_lower=False):
    if mask_lower:
        mask =  np.tri(matrix_in.shape[0], k=-1)
        matrix_in = np.ma.array(matrix_in, mask=mask) # mask out the lower triangle
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    cmap = cm.get_cmap('jet', 10) 
    cmap.set_bad('w') 
    cax = ax1.imshow(matrix_in, interpolation="nearest", cmap=cmap)
    plt.xlabel("lag (Mpc/h)")
    plt.ylabel("lag (Mpc/h)")
    if (len(axis_labels) > 0):
        print axis_labels
        ax1.set_xticks(range(len(axis_labels)))
        ax1.set_yticks(range(len(axis_labels)))
        ax1.set_xticklabels(axis_labels.astype(np.integer))
        ax1.set_yticklabels(axis_labels.astype(np.integer))
    #ax1.grid(True)
    c = fig.colorbar(cax)
    c.ax.set_ylabel("correlation")
    plt.savefig(filename)


def plot_contour(filename, fbins, lags, corr2D,
                 cross_power=True, title=None, coloraxis=[]):
    a = plt.figure()

    if cross_power:
        corr2D = corr2D * 1.e3
    else:
        corr2D = sp.sign(corr2D) * sp.sqrt(abs(corr2D)) * 1e3

    #a.set_figwidth(a.get_figwidth() / 3.0)
    if len(coloraxis) > 0:
        f = plt.contourf(lags, (fbins) / 1e6, corr2D, coloraxis)
    else:
        f = plt.contourf(lags, (fbins) / 1e6, corr2D)

    f.ax.set_xscale('log')
    f.ax.set_yscale('log')

    plt.axis('scaled')
    plt.xlim((0.05, 0.9))
    plt.ylim((0.8, 100))
    plt.xlabel("angular lag, $\sigma$ (degrees, 34$\cdotp$Mpc/h)")
    plt.ylabel("frequency lag, $\pi$ (MHz, 4.5$\cdotp$Mpc/h)")
    plt.title(title)
    #c = plt.colorbar(f, ticks=coloraxis)
    c = plt.colorbar(f)

    c.ax.set_ylabel("correlation (mK)")
    plt.savefig(filename)


# TODO implement error option
def plot_collapsed(filename, sep_lags, corr1D, errors=None, save_old=False,
                   plot_old=False, cross_power=True, title=None,
                   ylog=True):
    nbins = len(sep_lags)
    a = plt.figure()
    ax = plt.gca()
    if ylog:
        ax.set_yscale("log")
    ax.set_xscale("log")
    elin = 2
    msize = 6

    if cross_power:
        corr1D *= 1e3
    else:
        corr1D = sp.sign(corrf) * sp.sqrt(abs(corrf)) * 1e3

    plt.plot(sep_lags, corr1D, linestyle='None', marker='o',
                         color='b', markersize=msize) 
    if ylog:
        plt.plot(sep_lags, -corr1D, linestyle='None', marker='o',
                             color='r', markersize=msize) 

    # model
    t_lags = sp.arange(0.1, 100, 0.1)
    r0 = 5.5
    rb = 7.0
    t = sp.sqrt(((rb + t_lags) / r0)**(-1.8))
    if cross_power:
        t *= t 

    t = t * 0.15 / t[0]
    f = plt.plot(t_lags, t, marker='None', color='k', linestyle='-')

    if cross_power:
        plt.axis([1.5, 100, 0.0001, 10.])
    else:
        plt.axis([1.5, 100, 0.01, 500.0])
    
    if not ylog:
        plt.axis([1.5, 100, -0.05, 0.12])

    plt.xlabel('lag (Mpc/h)')
    plt.ylabel('correlation (mK)')
    plt.title(title)
    plt.savefig(filename)

if __name__ == '__main__':
    #compare_corr_one("opt_x_radio_mapA_noconv_fast.shelve", 
    #                 "opt_x_radio_mapA_noconv_fastest.shelve", print_params=False)

    #process_batch_correlations("run1_correlations.shelve", batch1_param)
    #process_batch_correlations("run2_correlations.shelve", batch2_param)
    #process_batch_correlations("run3_correlations.shelve", batch3_param)
    #process_batch_correlations("run4_correlations.shelve", batch4_param)
    #process_batch_correlations("run5_correlations.shelve", batch5_param)
    #process_batch_correlations("run6_correlations.shelve", batch6_param)
    #process_batch_correlations("run7_correlations.shelve", batch7_param)
    #process_batch_correlations("run5_correlations_newcorr.shelve", batch5_param)
    #process_batch_correlations("run7_correlations_newcorr.shelve", batch7_param)
    #process_batch_correlations("run8_correlations_modes.shelve", batch8_param)
    #process_batch_correlations("run9_correlations_modes.shelve", batch9_param)
    #process_batch_correlations("run10_correlations_modes.shelve", batch10_param)
    #process_batch_correlations("run11_correlations_modes.shelve", batch11_param)
    #process_batch_correlations("run12_correlations_modes.shelve", batch12_param)
    process_batch_correlations("run13_correlations.shelve", batch13_param)

    #print compare_corr(batch2_param, batch3_param)
    #print compare_corr(batch1_param, batch2_param)
    #print compare_corr(batch6_param, batch7_param)

    #plot_batch_correlations("run1_correlations.shelve", batch1_param,
    #                        dir_prefix="plots/run1b/")
    #plot_batch_correlations("run2_correlations.shelve", batch2_param,
    #                        dir_prefix="plots/run2/")
    #plot_batch_correlations("run3_correlations.shelve", batch3_param,
    #                        dir_prefix="plots/run3/")
    #plot_batch_correlations("run4_correlations.shelve", batch4_param,
    #                        dir_prefix="plots/run4/")
    #plot_batch_correlations("run5_correlations.shelve", batch5_param,
    #                        dir_prefix="plots/run5/",
    #                        color_range=[-0.04, 0.04])
    #plot_batch_correlations("run6_correlations.shelve", batch6_param,
    #                        dir_prefix="plots/run6/",
    #                        color_range=[-0.04, 0.04])
    #plot_batch_correlations("run7_correlations.shelve", batch7_param,
    #                        dir_prefix="plots/run7/",
    #                        color_range=[-0.04, 0.04])
    #plot_batch_correlations("run7_correlations_newcorr.shelve", batch7_param,
    #                        dir_prefix="plots/run7n/",
    #                        color_range=[-10, 10])
    #plot_batch_correlations("run8_correlations_modes.shelve", batch8_param,
    #                        dir_prefix="plots/run8/",
    #                        color_range=[-10, 10])
    #plot_batch_correlations("run9_correlations_modes.shelve", batch9_param,
    #                        dir_prefix="plots/run9/",
    #                        color_range=[-10, 10])
    #plot_batch_correlations("run10_correlations_modes.shelve", batch10_param,
    #                        dir_prefix="plots/run10/",
    #                        color_range=[-10, 10])
    #plot_batch_correlations("run11_correlations_modes.shelve", batch11_param,
    #                        dir_prefix="plots/run11/",
    #                        color_range=[-10, 10])
    #plot_batch_correlations("run12_correlations_modes.shelve", batch12_param,
    #                        dir_prefix="plots/run12/",
    #                        color_range=[-10, 10])
    plot_batch_correlations("run13_correlations.shelve", batch13_param,
                            dir_prefix="plots/run13/",
                            color_range=[-10, 10])

    #batch_correlations_statistics("run1_correlations.shelve", batch1_param)
    #batch_correlations_statistics("run2_correlations.shelve", batch2_param)
    #batch_correlations_statistics("run3_correlations.shelve", batch3_param)
    #batch_correlations_statistics("run4_correlations.shelve", batch4_param)
    #batch_correlations_statistics("run5_correlations.shelve", batch5_param)
    #batch_correlations_statistics("run6_correlations.shelve", batch6_param)
    #batch_correlations_statistics("run7_correlations.shelve", batch7_param)
    #batch_correlations_statistics("run5_correlations_newcorr.shelve", batch5_param)
    #batch_correlations_statistics("run7_correlations_newcorr.shelve", batch7_param)
    #batch_correlations_statistics("run10_correlations_modes.shelve", batch10_param)
    #batch_correlations_statistics("run11_correlations_modes.shelve", batch11_param)
    #batch_correlations_statistics("run12_correlations_modes.shelve", batch12_param)

