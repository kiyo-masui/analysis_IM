import numpy as np
import numpy.ma as ma
import shelve
from utils import data_paths as dp
from quadratic_products import pwrspec_estimator as pe
from utils import binning
from plotting import plot_slice
import copy


class PowerSpectrum(object):
    r"""
    All power spectra are in groups:

        autopower: 6 pairs, T treatments

        cross-power: 1 pair signal, N pairs sim; T treatments for all

        cleaning transfer function: N pairs sim, 1 pair signal, (T treatments);
        M pairs template

        beam transfer function: N pairs sim, M pairs template

    TODO:
    add mean counts to summary
    support operations on aggregate (mean + error) power spectra?
    support having a different trans. for each cross-pair (A with B, etc.)
    """

    def __init__(self, filename):
        pwrspec_compilation = shelve.open(filename, "r")
        case_key = "combination:treatment"
        self.pwr_cases = dp.unpack_cases(pwrspec_compilation.keys(),
                                         case_key, divider=":")

        # gather the lists if pair combinations and treatments (cleaning)
        self.comb_cases = self.pwr_cases["combination"]
        self.comb_cases.sort()
        self.treatment_cases = self.pwr_cases["treatment"]
        self.treatment_cases.sort()

        #print self.comb_cases, self.treatment_cases

        # the 2D power can be binned onto 1D
        # this is not done automatically because one can apply transfer
        # functions in 2D, etc. before this rebinning
        self.pwrspec_1d = {}
        self.pwrspec_1d_from_2d = {}
        self.pwrspec_2d = {}

        self.counts_1d = {}
        self.counts_1d_from_2d = {}
        self.gauss_errors_1d_from_2d = {}
        self.counts_2d = {}

        for pwrspec_case in pwrspec_compilation:
            pwrspec_entry = pwrspec_compilation[pwrspec_case]
            self.parameters = pwrspec_entry[0]
            pwrdata_2d = pwrspec_entry[1][0]
            pwrdata_1d = pwrspec_entry[1][1]

            self.pwrspec_1d[pwrspec_case] = pwrdata_1d["binavg"]
            self.counts_1d[pwrspec_case] = pwrdata_1d["counts_histo"]
            self.pwrspec_2d[pwrspec_case] = pwrdata_2d["binavg"]
            self.counts_2d[pwrspec_case] = pwrdata_2d["counts_histo"]

        self.k_1d = {}
        self.k_1d["left"] = pwrdata_1d["bin_left"]
        self.k_1d["center"] = pwrdata_1d["bin_center"]
        self.k_1d["right"] = pwrdata_1d["bin_right"]
        self.k_1d_from_2d = {}

        self.kx_2d = {}
        self.kx_2d["left"] = pwrdata_2d["bin_x_left"]
        self.kx_2d["center"] = pwrdata_2d["bin_x_center"]
        self.kx_2d["right"] = pwrdata_2d["bin_x_right"]
        self.ky_2d = {}
        self.ky_2d["left"] = pwrdata_2d["bin_y_left"]
        self.ky_2d["center"] = pwrdata_2d["bin_y_center"]
        self.ky_2d["right"] = pwrdata_2d["bin_y_right"]

        self.num_kx = pwrdata_2d["bin_x_center"].shape[0]
        self.num_ky = pwrdata_2d["bin_y_center"].shape[0]
        self.num_k_1d = pwrdata_1d["bin_center"].shape[0]
        self.num_k_1d_from_2d = None
        self.num_comb = len(self.comb_cases)
        self.num_treat = len(self.treatment_cases)

        pwrspec_compilation.close()

    def apply_2d_trans_by_treatment(self, transfer_dict):
        r"""here transfer_dict has keys for each treament and the value is the
        2d transfer function. Note that this changes the values saved in the
        class."""
        print "power_spectrum: applying 2D transfer function"
        for treatment in self.treatment_cases:
            transfer = transfer_dict[treatment]
            for comb in self.comb_cases:
                pwrcase = "%s:%s" % (comb, treatment)
                orig_counts = copy.deepcopy(self.counts_2d[pwrcase])

                self.pwrspec_2d[pwrcase] /= transfer
                self.counts_2d[pwrcase] *= transfer * transfer
                #print self.pwrspec_2d[pwrcase].shape, transfer.shape
                #print np.max(transfer)

                nanmask = np.isnan(self.counts_2d[pwrcase])
                infmask = np.isnan(self.counts_2d[pwrcase])
                zeromask = (orig_counts == 0)

                mask = np.where(np.logical_or(
                                np.logical_or(nanmask, infmask), zeromask))

                self.counts_2d[pwrcase][mask] = 0

    def apply_1d_trans_by_treatment(self, transfer_dict):
        r"""here transfer_dict has keys for each treament and the value is the
        1d transfer function. Note that this changes the values saved in the
        class."""
        for treatment in self.treatment_cases:
            transfer = transfer_dict[treatment]
            for comb in self.comb_cases:
                pwrcase = "%s:%s" % (comb, treatment)
                orig_counts = copy.deepcopy(self.counts_1d[pwrcase])

                self.pwrspec_1d[pwrcase] /= transfer
                self.counts_1d[pwrcase] *= transfer * transfer
                # if there is a 2d->1d binning in the class:
                if self.pwrspec_1d_from_2d:
                    self.pwrspec_1d_from_2d[pwrcase] /= transfer
                    self.counts_1d_from_2d[pwrcase] *= transfer * transfer

                nanmask = np.isnan(self.counts_1d[pwrcase])
                infmask = np.isnan(self.counts_1d[pwrcase])
                zeromask = (orig_counts == 0)

                mask = np.where(np.logical_or(
                                np.logical_or(nanmask, infmask), zeromask))

                self.counts_1d[pwrcase][mask] = 0
                if self.pwrspec_1d_from_2d:
                    self.counts_1d_from_2d[pwrcase][mask] = 0

    def convert_2d_to_1d(self, logbins=True,
                         bins=None, weights_2d=None):
        r"""bin the 2D powers onto 1D (the _from_2d variables)
        optionally, weights for each treatment class can be given
        """
        bins_1d = np.zeros(self.num_k_1d + 1)
        bins_1d[0: -1] = self.k_1d["left"]
        bins_1d[-1] = self.k_1d["right"][-1]

        if bins is None:
            bins = bins_1d

        for treatment in self.treatment_cases:
            if weights_2d is not None:
                weights_2d_treat = weights_2d[treatment]
            else:
                weights_2d_treat = None

            for comb in self.comb_cases:
                pwrcase = "%s:%s" % (comb, treatment)
                (self.counts_1d_from_2d[pwrcase], \
                self.gauss_errors_1d_from_2d[pwrcase], \
                self.pwrspec_1d_from_2d[pwrcase]) = \
                                        pe.convert_2d_to_1d_pwrspec(
                                            self.pwrspec_2d[pwrcase],
                                            self.counts_2d[pwrcase],
                                            self.kx_2d["center"],
                                            self.ky_2d["center"],
                                            bins,
                                            weights_2d=weights_2d_treat)

        bin_left, bin_center, bin_right = binning.bin_edges(bins, log=logbins)
        self.k_1d_from_2d["left"] = bin_left
        self.k_1d_from_2d["center"] = bin_center
        self.k_1d_from_2d["right"] = bin_right
        self.num_k_1d_from_2d = bin_center.shape[0]

    def combination_array_1d(self, from_2d=True, debug=False):
        r"""pack the various pair combinations for each treatment -> array
        """

        summary_pwr = {}
        summary_counts = {}
        summary_gerror = {}

        if from_2d:
            num_k = self.num_k_1d_from_2d
        else:
            num_k = self.num_k_1d

        for treatment in self.treatment_cases:
            pwr_treatment = np.zeros((num_k, self.num_comb))
            counts_treatment = np.zeros((num_k, self.num_comb))
            gerror_treatment = np.zeros((num_k, self.num_comb))

            for comb, comb_index in zip(self.comb_cases, range(self.num_comb)):
                pwrcase = "%s:%s" % (comb, treatment)
                if from_2d:
                    if debug:
                        print "combining 2D->1D counts"

                    counts_treatment[:, comb_index] = \
                            self.counts_1d_from_2d[pwrcase]

                    gerror_treatment[:, comb_index] = \
                            self.gauss_errors_1d_from_2d[pwrcase]

                    pwr_treatment[:, comb_index] = \
                            self.pwrspec_1d_from_2d[pwrcase]
                else:
                    if debug:
                        print "WARNING: combining 1D counts"

                    counts_treatment[:, comb_index] = \
                            self.counts_1d[pwrcase]

                    gerror_treatment[:, comb_index] =
                            np.zero_like(counts_treatment[:, comb_index])

                    pwr_treatment[:, comb_index] = \
                            self.pwrspec_1d[pwrcase]

            summary_counts[treatment] = counts_treatment
            summary_gerror[treatment] = gerror_treatment
            summary_pwr[treatment] = pwr_treatment

        return summary_counts, summary_gerror, summary_pwr

    def combination_array_2d(self, counts=False, debug=None):
        r"""pack the various pair combinations for each treatment into an array
        """

        summary_treatment = {}
        for treatment in self.treatment_cases:
            pwr_treatment = np.zeros((self.num_kx, self.num_ky,
                                      self.num_comb))

            for comb, comb_index in zip(self.comb_cases, range(self.num_comb)):
                pwrcase = "%s:%s" % (comb, treatment)
                if counts:
                    pwr_treatment[:, :, comb_index] = self.counts_2d[pwrcase]
                    if debug is not None:
                        outplot_power_file = "%s/counts_%s_%s.png" % \
                                             (debug, treatment, comb)

                else:
                    pwr_treatment[:, :, comb_index] = self.pwrspec_2d[pwrcase]
                    if debug is not None:
                        outplot_power_file = "%s/power_%s_%s.png" % \
                                             (debug, treatment, comb)

                if (debug is not None) and (treatment == "20modes"):
                    plot_slice.simpleplot_2D(outplot_power_file,
                                    np.abs(pwr_treatment[:, :, comb_index]),
                                    np.log10(self.kx_2d['center']),
                                    np.log10(self.ky_2d['center']),
                                    ["logkx", "logky"], 1., "2d array",
                                    "log(abs(array))", logscale=True)

            summary_treatment[treatment] = pwr_treatment

        return summary_treatment

    def agg_stat_1d_pwrspec(self, from_2d=False, use_masked=False):
        r"""TODO: use masked array here"""
        (counts_arr, gerror_arr, comb_arr) = \
            self.combination_array_1d(from_2d=from_2d)

        stat_summary = {}
        for treatment in self.treatment_cases:
            entry = {}
            if use_masked:
                gerror_ma = ma.masked_invalid(gerror_arr[treatment])
                counts_ma = ma.masked_invalid(counts_arr[treatment])
                comb_ma = ma.masked_invalid(comb_arr[treatment])

                if self.num_comb > 1:
                    entry['counts'] = np.ma.mean(counts_ma, axis=1)
                    entry['mean'] = np.ma.mean(comb_ma, axis=1)
                    entry['gauss_std'] = np.ma.mean(gerror_ma, axis=1)
                    entry['std'] = np.ma.std(comb_ma, axis=1, ddof=1)
                    entry['corr'] = np.ma.corrcoef(comb_ma)
                    entry['cov'] = np.ma.cov(comb_ma)
                else:
                    entry['counts'] = counts_ma[:, 0]
                    entry['mean'] = comb_ma[:, 0]
                    entry['gauss_std'] = gerror_ma[:, 0]
                    entry['std'] = 0.
                    entry['corr'] = 0.
                    entry['cov'] = 0.
            else:
                comb_treat = comb_arr[treatment]
                if self.num_comb > 1:
                    entry['counts'] = np.mean(counts_arr[treatment], axis=1)
                    entry['mean'] = np.mean(comb_treat, axis=1)
                    entry['gauss_std'] = np.mean(gerror_arr[treatment], axis=1)
                    entry['std'] = np.std(comb_treat, axis=1, ddof=1)
                    entry['corr'] = np.corrcoef(comb_treat)
                    entry['cov'] = np.cov(comb_treat)
                else:
                    entry['counts'] = counts_arr[treatment][:, 0]
                    entry['mean'] = comb_treat[:, 0]
                    entry['gauss_std'] = gerror_arr[treatment][:, 0]
                    entry['std'] = 0.
                    entry['corr'] = 0.
                    entry['cov'] = 0.

            stat_summary[treatment] = entry

        return stat_summary

    def agg_stat_2d_pwrspec(self, use_masked=False, debug=None):
        r"""TODO: add corr and cov"""
        comb_arr = self.combination_array_2d(counts=False, debug=debug)
        counts_arr = self.combination_array_2d(counts=True, debug=debug)

        stat_summary = {}
        for treatment in self.treatment_cases:
            entry = {}
            if use_masked:
                counts_ma = ma.masked_invalid(counts_arr[treatment])
                comb_ma = ma.masked_invalid(comb_arr[treatment])
                entry['counts'] = np.ma.mean(counts_ma, axis=2)
                entry['mean'] = np.ma.mean(comb_ma, axis=2)
                if self.num_comb > 1:
                    entry['std'] = np.ma.std(comb_ma, axis=2, ddof=1)
                else:
                    entry['std'] = 0.
            else:
                entry['counts'] = np.mean(counts_arr[treatment], axis=2)
                entry['mean'] = np.mean(comb_arr[treatment], axis=2)
                if self.num_comb > 1:
                    entry['std'] = np.std(comb_arr[treatment], axis=2, ddof=1)
                else:
                    entry['std'] = 0.
            stat_summary[treatment] = entry

        return stat_summary

    def summarize_1d_pwrspec(self, pwr_1d, filename):
        r"""Write out a 1d power spectrum
        """
        outfile = open(filename, "w")
        for specdata in zip(self.bin_left, self.bin_center,
                            self.bin_right, self.pwr_1d['binavg']):
            outfile.write(("%10.15g " * 4 + "\n") % specdata)
        outfile.close()

    def summarize_1d_agg_pwrspec(self, treatment, filename, corr_file=None,
                                 from_2d=False):
        # TODO: add counts histo to the summary
        # this could be made more efficient by looping over treatments
        # or calculating values -> self
        stat_summary = self.agg_stat_1d_pwrspec(from_2d=from_2d)
        mean_1d = stat_summary[treatment]["mean"]
        std_1d = stat_summary[treatment]["std"]
        gauss_std_1d = stat_summary[treatment]["gauss_std"]
        corrmat_1d = stat_summary[treatment]["corr"]

        outfile = open(filename, "w")
        for specdata in zip(self.bin_left, self.bin_center,
                            self.bin_right, mean_1d,
                            std_1d, gauss_std_1d):
            outfile.write(("%10.15g " * 6 + "\n") % specdata)
        outfile.close()

        bin_left_lt = np.log10(self.bin_left)
        bin_center_lt = np.log10(self.bin_center)
        bin_right_lt = np.log10(self.bin_right)

        if corr_file is not None:
            outfile = open(corr_file, "w")
            for xind in range(len(bin_center)):
                for yind in range(len(bin_center)):
                    outstr = ("%10.15g " * 7 + "\n") % \
                            (bin_left_lt[xind], bin_center_lt[xind], \
                             bin_right_lt[xind], bin_left_lt[yind], \
                             bin_center_lt[yind], bin_right_lt[yind], \
                             corrmat_1d[xind, yind])
                    outfile.write(outstr)

            outfile.close()

    def summarize_1d_pwrspec_by_treatment(self, basename, from_2d=False):
        for treatment in self.treatment_cases:
            print "WRITE THIS!"

    def summarize_2d_pwrspec(self, pwr_2d, filename, resetnan=0.):
        r"""Write out a single pwrspec
        """
        lbin_x_left = np.log10(self.bin_x_left)
        lbin_x_center = np.log10(self.bin_x_center)
        lbin_x_right = np.log10(self.bin_x_right)
        lbin_y_left = np.log10(self.bin_y_left)
        lbin_y_center = np.log10(self.bin_y_center)
        lbin_y_right = np.log10(self.bin_y_right)

        outfile = open(filename, "w")
        reset_2d = copy.deepcopy(pwr_2d)
        reset_2d[np.isnan(reset_2d)] = resetnan
        for xind in range(len(lbin_x_center)):
            for yind in range(len(lbin_y_center)):
                outstr = ("%10.15g " * 7 + "\n") % \
                        (lbin_x_left[xind], lbin_x_center[xind], \
                         lbin_x_right[xind], lbin_y_left[yind], \
                         lbin_y_center[yind], lbin_y_right[yind], \
                         reset_2d[xind, yind])
                outfile.write(outstr)

        outfile.close()

    def summarize_2d_agg_pwrspec(self, treatment, filename, resetnan=0.):
        r"""Combine a list of 2D power runs and write out
        TODO: add 2D corr
        """
        stat_summary = self.agg_stat_2d_pwrspec()
        mean_2d = stat_summary[treatment]["mean"]
        std_2d = stat_summary[treatment]["std"]
        counts_2d = stat_summary[treatment]["counts"]

        bin_x_left = np.log10(self.kx_2d["left"])
        bin_x_center = np.log10(self.kx_2d["center"])
        bin_x_right = np.log10(self.kx_2d["right"])
        bin_y_left = np.log10(self.ky_2d["left"])
        bin_y_center = np.log10(self.ky_2d["center"])
        bin_y_right = np.log10(self.ky_2d["right"])

        plot_mean_2d = copy.deepcopy(mean_2d)
        plot_std_2d = copy.deepcopy(std_2d)
        plot_mean_2d[np.isnan(mean_2d)] = resetnan
        plot_std_2d[np.isnan(std_2d)] = resetnan
        outfile = open(filename, "w")
        for xind in range(len(bin_x_center)):
            for yind in range(len(bin_y_center)):
                outstr = ("%10.15g " * 9 + "\n") % \
                        (bin_x_left[xind], bin_x_center[xind],
                         bin_x_right[xind], \
                         bin_y_left[yind], bin_y_center[yind],
                         bin_y_right[yind], \
                         counts_2d[xind, yind], plot_mean_2d[xind, yind],
                         plot_std_2d[xind, yind])

                outfile.write(outstr)

        outfile.close()

        return mean_2d, std_2d

    def summarize_pwrspec(self, tag, outdir="./plot_data"):
        r"""Plot the 1D and 2D power spectra from a run
        """
        fileout = outdir + "/" + tag + "_from_2d.dat"
        summarize_1d_pwrspec(pwr_1d_from_2d, fileout)

        fileout = outdir + "/" + tag + ".dat"
        summarize_1d_pwrspec(pwr_1d, fileout)

        fileout = outdir + "/" + tag + "_2d.dat"
        summarize_2d_pwrspec(pwr_2d, fileout, dataname="binavg")

        fileout = outdir + "/" + tag + "_2d_counts.dat"
        summarize_2d_pwrspec(pwr_2d, fileout, dataname="counts_histo")

    def summarize_agg_pwrspec(self, treatment, tag, outdir="./plot_data"):
        r"""Unclear what this will do here"""
        fileout = outdir + "/" + tag + "_avg_from_2d.dat"
        corr_fileout = outdir + "/" + tag + "_corr_from_2d.dat"
        agg_1d_pwrspec_f2d = summarize_1d_agg_pwrspec(pwr_1d_from_2d,
                                    fileout,
                                    corr_file=corr_fileout)

        fileout = outdir + "/" + tag + "_avg.dat"
        corr_fileout = outdir + "/" + tag + "_corr.dat"
        agg_1d_pwrspec = summarize_1d_agg_pwrspec(pwr_1d, fileout,
                                                 corr_file=corr_fileout)

        fileout = outdir + "/" + tag + "_avg_2d.dat"
        summarize_2d_agg_pwrspec(pwr_2d, fileout, dataname="binavg")

        fileout = outdir + "/" + tag + "_avg_2d_counts.dat"
        summarize_2d_agg_pwrspec(pwr_2d, fileout, dataname="counts_histo")

        return agg_1d_pwrspec_f2d

    def print_param(self):
        print self.parameters
