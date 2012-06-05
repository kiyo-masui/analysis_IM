def repackage_1d_power(pwrspec_compilation):
    case_key = "combination:treatment"
    pwr_cases = dp.unpack_cases(pwrspec_compilation.keys(), case_key, divider=":")

    comb_cases = pwr_cases["combination"]
    # note that the dictionary may not be in any order
    comb_cases.sort()

    summary = {}
    for pwrspec_case in pwrspec_compilation:
        pwrspec_entry = pwrspec_compilation[pwrspec_case]
        parameters = pwrspec_entry[0]
        pwrdata_2d = pwrspec_entry[1][0]
        pwrdata_1d = pwrspec_entry[1][1]
        summary[pwrspec_case] = pwrdata_1d["binavg"]
        bin_left = pwrdata_1d["bin_left"]
        bin_center = pwrdata_1d["bin_center"]
        bin_right = pwrdata_1d["bin_right"]
        histo_counts = pwrdata_1d["counts_histo"]
        #print pwrspec_case, parameters

    num_k = bin_center.shape[0]
    num_pk = len(comb_cases)

    summary_treatment = {}
    for treatment in pwr_cases["treatment"]:
        pwr_treatment = np.zeros((num_k, num_pk))
        for comb, comb_index in zip(comb_cases, range(num_pk)):
            pwrcase = "%s:%s" % (comb, treatment)
            pwr_treatment[:, comb_index] = summary[pwrcase]
        summary_treatment[treatment] = pwr_treatment

    return (bin_left, bin_center, bin_right, histo_counts, summary_treatment)


def convert_2d_to_1d(pwrspec2d_product, logbins=True,
                     bins=None, transfer=None):
    """if bins is not given, just use the x axis"""
    nxbins = len(pwrspec2d_product['bin_x_center'])
    bins_kx = np.zeros(nxbins + 1)
    bins_kx[0: -1] = pwrspec2d_product['bin_x_left']
    bins_kx[-1] = pwrspec2d_product['bin_x_right'][-1]

    nybins = len(pwrspec2d_product['bin_y_center'])
    bins_ky = np.zeros(nybins + 1)
    bins_ky[0: -1] = pwrspec2d_product['bin_y_left']
    bins_ky[-1] = pwrspec2d_product['bin_y_right'][-1]

    if bins is None:
        bins = bins_kx

    entry = {}
    (entry['counts_histo'], entry['binavg']) = \
                        convert_2d_to_1d_driver(pwrspec2d_product['binavg'],
                        pwrspec2d_product['counts_histo'],
                        pwrspec2d_product['bin_x_center'],
                        pwrspec2d_product['bin_y_center'],
                        bins, transfer=transfer)

    bin_left, bin_center, bin_right = binning.bin_edges(bins, log=logbins)
    entry['bin_left'] = bin_left
    entry['bin_center'] = bin_center
    entry['bin_right'] = bin_right

    return entry

def convert_2d_to_1d_driver(pwr_2d, counts_2d, bin_kx, bin_ky, bin_1d,
                            transfer=None):
    """take a 2D power spectrum and the counts matrix (number of modex per k
    cell) and return the binned 1D power spectrum
    pwr_2d is the 2D power
    counts_2d is the counts matrix
    bin_kx is the x-axis
    bin_ky is the x-axis
    bin_1d is the k vector over which to return the result
    """
    # find |k| across the array
    index_array = np.indices(pwr_2d.shape)
    scale_array = np.zeros(index_array.shape)
    scale_array[0, ...] = bin_kx[index_array[0, ...]]
    scale_array[1, ...] = bin_ky[index_array[1, ...]]
    scale_array = np.rollaxis(scale_array, 0, scale_array.ndim)
    radius_array = np.sum(scale_array ** 2., axis=-1) ** 0.5

    radius_flat = radius_array.flatten()
    pwr_2d_flat = pwr_2d.flatten()
    counts_2d_flat = counts_2d.flatten()

    if transfer is not None:
        orig_counts = copy.deepcopy(counts_2d_flat)
        trans_flat = transfer.flatten()
        pwr_2d_flat /= trans_flat
        counts_2d_flat *= trans_flat * trans_flat
        counts_2d_flat[np.isnan(counts_2d_flat)] = 0
        counts_2d_flat[np.isinf(counts_2d_flat)] = 0
        counts_2d_flat[orig_counts == 0] = 0

    count_pwr_prod = counts_2d_flat * pwr_2d_flat
    count_pwr_prod[np.isnan(count_pwr_prod)] = 0
    count_pwr_prod[np.isinf(count_pwr_prod)] = 0
    count_pwr_prod[counts_2d_flat == 0] = 0

    counts_histo = np.histogram(radius_flat, bin_1d,
                                weights=counts_2d_flat)[0]
    binsum_histo = np.histogram(radius_flat, bin_1d,
                                weights=count_pwr_prod)[0]

    binavg = binsum_histo / counts_histo.astype(float)

    return counts_histo, binavg

def summarize_agg_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d,
                          tag, outdir="./plot_data", apply_1d_transfer=None):
    r"""call various power spectral aggregation functions to make 2D, 1D, etc.
    P(k)s to plot.
    """
    fileout = outdir + "/" + tag + "_avg_from2d.dat"
    corr_fileout = outdir + "/" + tag + "_corr_from2d.dat"
    agg_1d_pwrspec_f2d = summarize_1d_agg_pwrspec(pwr_1d_from_2d,
                                fileout,
                                corr_file=corr_fileout,
                                apply_1d_transfer=apply_1d_transfer)

    fileout = outdir + "/" + tag + "_avg.dat"
    corr_fileout = outdir + "/" + tag + "_corr.dat"
    agg_1d_pwrspec = summarize_1d_agg_pwrspec(pwr_1d, fileout,
                                             corr_file=corr_fileout,
                                apply_1d_transfer=apply_1d_transfer)

    fileout = outdir + "/" + tag + "_avg_2d.dat"
    summarize_2d_agg_pwrspec(pwr_2d, fileout, dataname="binavg")

    fileout = outdir + "/" + tag + "_avg_2d_counts.dat"
    summarize_2d_agg_pwrspec(pwr_2d, fileout, dataname="counts_histo")

    return agg_1d_pwrspec_f2d

def summarize_2d_pwrspec(pwr_2d, filename, dataname='binavg', resetnan=0.):
    r"""Write out a single pwrspec
    """
    bin_x_left = np.log10(pwr_2d['bin_x_left'])
    bin_x_center = np.log10(pwr_2d['bin_x_center'])
    bin_x_right = np.log10(pwr_2d['bin_x_right'])
    bin_y_left = np.log10(pwr_2d['bin_y_left'])
    bin_y_center = np.log10(pwr_2d['bin_y_center'])
    bin_y_right = np.log10(pwr_2d['bin_y_right'])

    outfile = open(filename, "w")
    reset_2d = copy.deepcopy(pwr_2d[dataname])
    reset_2d[np.isnan(reset_2d)] = resetnan
    for xind in range(len(bin_x_center)):
        for yind in range(len(bin_y_center)):
            outstr = ("%10.15g " * 7 + "\n") % \
                    (bin_x_left[xind], bin_x_center[xind], bin_x_right[xind], \
                     bin_y_left[yind], bin_y_center[yind], bin_y_right[yind], \
                     reset_2d[xind, yind])
            outfile.write(outstr)

    outfile.close()

def summarize_1d_pwrspec(pwr_1d, filename):
    r"""Write out a 1d power spectrum
    """
    bin_left = pwr_1d['bin_left']
    bin_center = pwr_1d['bin_center']
    bin_right = pwr_1d['bin_right']

    outfile = open(filename, "w")
    for specdata in zip(bin_left, bin_center,
                        bin_right, pwr_1d['binavg']):
        outfile.write(("%10.15g " * 4 + "\n") % specdata)
    outfile.close()

def agg_stat_1d_pwrspec(pwr_1d, apply_1d_transfer=None):
    pwrshp_1d = pwr_1d[0]['binavg'].shape
    n_runs = len(pwr_1d)

    pwrmat_1d = np.zeros((n_runs, pwrshp_1d[0]))
    for index in range(n_runs):
        if apply_1d_transfer is not None:
            print apply_1d_transfer
            pwrmat_1d[index, :] = pwr_1d[index]['binavg'] / apply_1d_transfer
        else:
            pwrmat_1d[index, :] = pwr_1d[index]['binavg']

    mean_1d = np.mean(pwrmat_1d, axis=0)
    std_1d = np.std(pwrmat_1d, axis=0, ddof=1)
    corrmat_1d = np.corrcoef(np.transpose(pwrmat_1d))
    covmat_1d = np.cov(np.transpose(pwrmat_1d))

    return (mean_1d, std_1d, corrmat_1d, covmat_1d)

def summarize_1d_agg_pwrspec(pwr_1d, filename, corr_file=None,
                             apply_1d_transfer=None):
    r"""Summarize the 1D power spectrum from a list of one-dimensional power
    spectrum outputs.
    """
    (mean_1d, std_1d, corrmat_1d, covmat_1d) = agg_stat_1d_pwrspec(pwr_1d,
                                      apply_1d_transfer=apply_1d_transfer)

    # assume that they all have the same binning
    bin_left = pwr_1d[0]['bin_left']
    bin_center = pwr_1d[0]['bin_center']
    bin_right = pwr_1d[0]['bin_right']
    counts_histo = pwr_1d[0]['counts_histo']

    outfile = open(filename, "w")
    for specdata in zip(bin_left, bin_center,
                        bin_right, counts_histo, mean_1d, std_1d):
        outfile.write(("%10.15g " * 6 + "\n") % specdata)
    outfile.close()

    bin_left_lt = np.log10(bin_left)
    bin_center_lt = np.log10(bin_center)
    bin_right_lt = np.log10(bin_right)

    if corr_file is not None:
        outfile = open(corr_file, "w")
        for xind in range(len(bin_center)):
            for yind in range(len(bin_center)):
                outstr = ("%10.15g " * 7 + "\n") % \
                        (bin_left_lt[xind], bin_center_lt[xind], bin_right_lt[xind], \
                         bin_left_lt[yind], bin_center_lt[yind], bin_right_lt[yind], \
                         corrmat_1d[xind, yind])
                outfile.write(outstr)

        outfile.close()

    retval = {}
    retval["mean_1d"] = mean_1d
    retval["std_1d"] = std_1d
    retval["covmat_1d"] = covmat_1d
    retval["bin_left"] = bin_left
    retval["bin_center"] = bin_center
    retval["bin_right"] = bin_right
    retval["counts_histo"] = counts_histo

    return retval

def agg_stat_2d_pwrspec(pwr_2d, dataname='binavg'):
    pwrshp_2d = pwr_2d[0][dataname].shape
    n_runs = len(pwr_2d)

    pwrmat_2d = np.zeros((n_runs, pwrshp_2d[0], pwrshp_2d[1]))
    for index in range(n_runs):
        pwrmat_2d[index, :, :] = pwr_2d[index][dataname]

    mean_2d = np.mean(pwrmat_2d, axis=0)
    std_2d = np.std(pwrmat_2d, axis=0, ddof=1)

    return (mean_2d, std_2d)


def summarize_2d_agg_pwrspec(pwr_2d, filename, dataname='binavg', resetnan=0.):
    r"""Combine a list of 2D power runs and write out
    """
    (mean_2d, std_2d) = agg_stat_2d_pwrspec(pwr_2d, dataname=dataname)

    bin_x_left = np.log10(pwr_2d[0]['bin_x_left'])
    bin_x_center = np.log10(pwr_2d[0]['bin_x_center'])
    bin_x_right = np.log10(pwr_2d[0]['bin_x_right'])
    bin_y_left = np.log10(pwr_2d[0]['bin_y_left'])
    bin_y_center = np.log10(pwr_2d[0]['bin_y_center'])
    bin_y_right = np.log10(pwr_2d[0]['bin_y_right'])

    plot_mean_2d = copy.deepcopy(mean_2d)
    plot_std_2d = copy.deepcopy(std_2d)
    plot_mean_2d[np.isnan(mean_2d)] = resetnan
    plot_std_2d[np.isnan(std_2d)] = resetnan
    outfile = open(filename, "w")
    for xind in range(len(bin_x_center)):
        for yind in range(len(bin_y_center)):
            outstr = ("%10.15g " * 8 + "\n") % \
                    (bin_x_left[xind], bin_x_center[xind], bin_x_right[xind], \
                     bin_y_left[yind], bin_y_center[yind], bin_y_right[yind], \
                     plot_mean_2d[xind, yind], plot_std_2d[xind, yind])
            outfile.write(outstr)

    outfile.close()

    return mean_2d, std_2d


def summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d,
                      tag, outdir="./plot_data"):
    r"""Plot the 1D and 2D power spectra from a run
    """
    fileout = outdir + "/" + tag + "_from2d.dat"
    summarize_1d_pwrspec(pwr_1d_from_2d, fileout)

    fileout = outdir + "/" + tag + ".dat"
    summarize_1d_pwrspec(pwr_1d, fileout)

    fileout = outdir + "/" + tag + "_2d.dat"
    summarize_2d_pwrspec(pwr_2d, fileout, dataname="binavg")

    fileout = outdir + "/" + tag + "_2d_counts.dat"
    summarize_2d_pwrspec(pwr_2d, fileout, dataname="counts_histo")





