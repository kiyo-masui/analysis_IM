"""Make summary plots of the binned correlation functions for large sets of
random catalogs """
import os
import re
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
from core import algebra
import multiprocessing
from kiyopy import parse_ini
# TODO: convert batch params into ini files
# TODO: make sure all methods here used counts/weights as-saved
# TODO: new methods are fast enough that no need to partition binning plotting,
# and statistics; still use multiprocessing


def make_corr(filename, verbose=False, identifier=None, multiplier=1.):
    """wrap the plot correlation class which reads correlation object shelve
    files; uses new binning methods in freq-slices"""
    output = {}
    corr_shelve = shelve.open(filename + ".shelve")

    if (multiplier != 1.):
        print "WARNING: using a multiplier of: " + repr(multiplier)

    corr = corr_shelve["corr"]*multiplier
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
                                     weights=corr_counts, return_fbins=True,
                                     nfbins=200)
    corr_1D = fs.collapse_correlation_1D(corr_2D[0], corr_2D[2], lags,
                                         weights=corr_2D[1])

    output["run_params"] = run_params
    output["lags"] = run_params["lags"]
    # uncomment these only if you need them in the shelve file; makes it huge
    #output["corr"] = corr
    #output["corr_counts"] = corr_counts
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


def plot_corr(shelve_entry, filename, title, coloraxis=None, cross_power=True):
    """make plots for a single correlation function (internal)"""

    print "plotting " + filename
    outlog = open(filename + ".txt", 'w')
    run_params = shelve_entry["run_params"]
    outlog.write("run parameters \n" + "-" * 80 + "\n")
    for key in run_params:
        outlog.write(key + ": " + repr(run_params[key]) + "\n")
    outlog.write("\n")

    outlog.write("binned correlation function \n" + "-" * 80 + "\n")
    for (lag, cdat) in zip(shelve_entry["corr1D_lags"],
                           shelve_entry["corr1D"]):
        outlog.write(repr(lag) + repr(cdat) + "\n")

    plot_collapsed(filename + ".png", shelve_entry["corr1D_lags"],
                                 shelve_entry["corr1D"],
                                 cross_power=cross_power, title=title)

    plot_contour(filename + "_contour.png", shelve_entry["corr2D_fbins"],
                 shelve_entry["lags"], shelve_entry["corr2D"],
                 cross_power=cross_power, title=title, coloraxis=coloraxis)


def make_shelve_names(batch_param, multiplier=1.):
    """assemble a list of shelve file names containing the correlation function
    information; lump a multiplier in with each file
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
                    filelist.append((listid, fullpath, multiplier))

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
    (run_id, run_file, multiplier) = runitem
    return make_corr(run_file, identifier=run_id, multiplier=multiplier)


def process_batch_correlations(filename, batch_param, multiplier=1.):
    """Process a batch of correlation functions"""
    product = shelve.open(filename)
    filelist = make_shelve_names(batch_param, multiplier=multiplier)

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = pool.map(wrap_make_corr, filelist)
    # TODO: can pool write this dictionary instead?
    for item in results:
        product[item["identifier"]] = item

    product.close()


def repair_shelve_files(batch_param, ini_prefix, params_default, param_prefix):
    """Add missing information to shelves"""
    filelist = make_shelve_names(batch_param)
    for (index, filename) in filelist:
        print "repairing: " + filename
        directory = "/".join(filename.split("/")[0:-1]) + "/"
        run_index = re.findall(r'\d+', index)[0]
        ini_file = directory + ini_prefix + run_index + ".ini"
        print ini_file
        params = parse_ini.parse(ini_file, params_default,
                             prefix=param_prefix, feedback=10)

        radio_file1 = params['radio_root1'] + params['radio_data_file1']
        map_radio1 = algebra.make_vect(algebra.load(radio_file1))

        corr_data = shelve.open(filename + ".shelve")
        corr_data["params"] = params
        corr_data["freq_axis"] = map_radio1.get_axis('freq')
        corr_data.close()


def wrap_plot_corr(runitem):
    """wrapper to the correlation function plotter for the process pool"""
    (shelve_entry, filename, title, coloraxis, cross_power) = runitem
    plot_corr(shelve_entry, filename, title,
              coloraxis=coloraxis, cross_power=cross_power)


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
    for (run_id, run_file, multiplier) in randlist:
        print run_file
        shelve_entry = master[run_id]
        rancats[index, :] = shelve_entry["corr1D"]
        #ranaxis = shelve_entry["corr1D_lags"]
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


def fancy_vector(vector, format_string):
    """print a numpy vector with a format string"""
    output = ""
    for entry in vector:
        output += (format_string + " ") % entry

    return output


def batch_compensation_function(filename, batch_param, modetoken="mode"):
    """find the impact of filtering on a run of correlation functions"""
    master = shelve.open(filename)
    filelist = make_shelve_names(batch_param)
    # TODO remove 10 magic number and have it figure this out instead
    nlags = 10

    # make a sublist of calculated correlations for just the random trials
    modelist = []
    for item in filelist:
        if string.find(item[0], modetoken) != -1:
            modelist.append(item)

    print "number of mode subtraction runs: " + repr(len(modelist))
    modecats = np.zeros((len(modelist), nlags))
    index = 0
    for (run_id, run_file, multiplier) in modelist:
        print run_file
        shelve_entry = master[run_id]
        modecats[index, :] = shelve_entry["corr1D"]
        modeaxis = shelve_entry["corr1D_lags"]
        index += 1

    compmode = np.zeros((len(modelist), nlags))
    for modeindex in range(len(modelist)):
        compmode[modeindex, :] = modecats[modeindex, :]/modecats[0, :]

    print modeaxis
    for lagindex in range(nlags):
        lag = modeaxis[lagindex]
        modeloss = compmode[:, lagindex]
        print int(lag), fancy_vector(modeloss, '%5.2g')

def plot_batch_correlations(filename, batch_param, dir_prefix="plots/",
                            color_range=[-0.2, 0.2], cross_power=True):
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
    for (run_id, run_file, multiplier) in filelist:
        shelve_entry = master[run_id]
        run_tag = run_file.split('/')
        run_tag = run_tag[-1]
        runlist.append((shelve_entry, dir_prefix + run_tag,
                        run_tag, coloraxis, cross_power))

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(wrap_plot_corr, runlist)

    master.close()


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
        corr1Dplt = corr1D*1e3
    else:
        corr1Dplt = sp.sign(corr1D) * sp.sqrt(abs(corr1D)) * 1e3

    plt.plot(sep_lags, corr1Dplt, linestyle='None', marker='o',
                         color='b', markersize=msize)
    if ylog:
        plt.plot(sep_lags, -corr1Dplt, linestyle='None', marker='o',
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
        plt.axis([1.5, 100, 0.01, 500.0])
        #plt.axis([1.5, 100, 0.0001, 10.])
    else:
        plt.axis([1.5, 100, 0.01, 500.0])
        #plt.axis([1.5, 100, 0.0001, 10.])

    if not ylog:
        plt.axis([1.5, 100, -0.05, 0.12])

    plt.xlabel('lag (Mpc/h)')
    plt.ylabel('correlation (mK)')
    plt.title(title)
    plt.savefig(filename)

