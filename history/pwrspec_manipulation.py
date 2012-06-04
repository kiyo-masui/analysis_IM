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


