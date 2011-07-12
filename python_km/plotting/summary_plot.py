"""Make summary plots of the binned correlation functions for large sets of
random catalogs """
import os
import string
import shelve
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
from correlate import correlation_plots as cp
import multiprocessing
rootdir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/batch_runs/"

# note that notes are purely human-readable and the keys do not mean anything
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


def make_corr(filename, verbose=False, identifier=None):
    """wrap the plot correlation class which reads correlation object shelve
    files"""
    output = {}
    corr_shelve = shelve.open(filename + ".shelve")

    corr = corr_shelve["corr"]
    run_params = corr_shelve["params"]

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


def compare_corr(batchlist_a, batchlist_b, print_params=False):
    """Compare two sets of correlation function runs"""
    filedict_a = tuple_list_to_dict(make_shelve_names(batchlist_a))
    filedict_b = tuple_list_to_dict(make_shelve_names(batchlist_b))

    common_ids = set(filedict_a.keys()).union(set(filedict_b.keys()))

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

        difference = False
        if file_a and file_b:
            corr_a_shelve = shelve.open(file_a)
            corr_b_shelve = shelve.open(file_b)
            corr_a = corr_a_shelve["corr"]
            corr_b = corr_b_shelve["corr"]
            delta = np.max(corr_a - corr_b)
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


if __name__ == '__main__':
    #process_batch_correlations("run1_correlations.shelve", batch1_param)
    #process_batch_correlations("run2_correlations.shelve", batch2_param)
    #process_batch_correlations("run3_correlations.shelve", batch3_param)
    #process_batch_correlations("run4_correlations.shelve", batch4_param)
    #process_batch_correlations("run5_correlations.shelve", batch5_param)

    #print compare_corr(batch2_param, batch3_param)
    #print compare_corr(batch1_param, batch2_param)

    #plot_batch_correlations("run1_correlations.shelve", batch1_param,
    #                        dir_prefix="plots/run1b/")
    #plot_batch_correlations("run2_correlations.shelve", batch2_param,
    #                        dir_prefix="plots/run2/")
    #plot_batch_correlations("run3_correlations.shelve", batch3_param,
    #                        dir_prefix="plots/run3/")
    #plot_batch_correlations("run4_correlations.shelve", batch4_param,
    #                        dir_prefix="plots/run4/")
    plot_batch_correlations("run5_correlations.shelve", batch5_param,
                            dir_prefix="plots/run5/",
                            color_range=[-0.04, 0.04])

    #batch_correlations_statistics("run1_correlations.shelve", batch1_param)
    #batch_correlations_statistics("run2_correlations.shelve", batch2_param)
    #batch_correlations_statistics("run3_correlations.shelve", batch3_param)
    #batch_correlations_statistics("run4_correlations.shelve", batch4_param)
    #batch_correlations_statistics("run5_correlations.shelve", batch5_param)
