import numpy as np
import shelve
from utils import data_paths as dp

def convert_2d_to_1d_driver(pwr_2d, counts_2d, bin_kx, bin_ky, bin_1d):
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


class PowerSpectrum(object):
    r"""
    All power spectra are in groups:
        autopower: 6 pairs, T treatments
        cross-power: 1 pair signal, N pairs sim; T treatments for all
        cleaning transfer function: N pairs sim, 1 pair signal, (T treatments); M pairs template
        beam transfer function: N pairs sim, M pairs template
    """

    def __init__(self, filename):
        pwrspec_compilation = shelve.open(filename, "r")
        case_key = "combination:treatment"
        self.pwr_cases = dp.unpack_cases(pwrspec_compilation.keys(), case_key, divider=":")

        # gather the lists if pair combinations and treatments (cleaning)
        self.comb_cases = self.pwr_cases["combination"]
        self.comb_cases.sort()
        self.treatment_cases = self.pwr_cases["treatment"]
        self.treatment_cases.sort()

        self.pwrspec_1d = {}
        self.pwrspec_2d = {}
        self.counts_1d = {}
        self.counts_2d = {}
        for pwrspec_case in pwrspec_compilation:
            pwrspec_entry = pwrspec_compilation[pwrspec_case]
            self.parameters = pwrspec_entry[0]
            pwrdata_2d = pwrspec_entry[1][0]
            pwrdata_1d = pwrspec_entry[1][1]

            self.pwrspec_1d[pwrspec_case] = pwrdata_1d["binavg"]
            self.counts_1d[pwrspec_case] = pwrdata_1d["counts_histo"]
            self.bin_left_1d = pwrdata_1d["bin_left"]
            self.bin_center_1d = pwrdata_1d["bin_center"]
            self.bin_right_1d = pwrdata_1d["bin_right"]

            self.pwrspec_2d[pwrspec_case] = pwrdata_2d["binavg"]
            self.counts_2d[pwrspec_case] = pwrdata_2d["counts_histo"]
            self.bin_left_x = pwrdata_2d["bin_x_left"]
            self.bin_center_x = pwrdata_2d["bin_x_center"]
            self.bin_right_x = pwrdata_2d["bin_x_right"]
            self.bin_left_y = pwrdata_2d["bin_y_left"]
            self.bin_center_y = pwrdata_2d["bin_y_center"]
            self.bin_right_y = pwrdata_2d["bin_y_right"]

        self.num_k = self.bin_center_1d.shape[0]
        self.num_comb = len(self.comb_cases)
        self.num_treat = len(self.treatment_cases)

        # the 2D power can be binned onto 1D
        # this is not done automatically because one can apply transfer
        # functions in 2D, etc. before this rebinning
        self.pwrspec_1d_from_2d = {}
        self.counts_1d_from_2d = {}
        self.bin_left_1d_from_2d = None
        self.bin_center_1d_from_2d = None
        self.bin_right_1d_from_2d = None

        pwrspec_compilation.close()

    def apply_2d_trans_by_treatment(self, transfer_dict):
        r"""here transfer_dict has keys for each treament and the value is the
        2d transfer function. Note that this changes the values saved in the
        class."""
        for treatment in self.treatment_cases:
            transfer = transfer_dict[treatment]
            for comb in self.comb_cases:
                pwrcase = "%s:%s" % (comb, treatment)
                orig_counts = copy.deepcopy(self.counts_2d[pwrcase])

                self.pwrspec_2d[pwrcase] /= transfer
                self.counts_2d[pwrcase] *= transfer * transfer

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

    def convert_2d_to_1d(pwrspec2d_product, logbins=True,
                     bins=None, transfer=None):
        r"""bin the 2D powers onto 1D (the _from_2d variables)"""
        bins = self.bin_center_1d
        for treatment in self.treatment_cases:
            for comb in self.comb_cases:
                pwrcase = "%s:%s" % (comb, treatment)
                (self.counts_1d_from_2d[pwrcase], \
                 self.pwrspec_1d_from_2d[pwrcase]) = convert_2d_to_1d_driver(
                                                     self.pwrspec_2d[pwrcase],
                                                     self.counts_2d[pwrcase],
                                                     self.bin_center_x,
                                                     self.bin_center_y,
                                                     bins)

        bin_left, bin_center, bin_right = binning.bin_edges(bins, log=logbins)
        self.bin_left_1d_from_2d = bin_left
        self.bin_center_1d_from_2d = bin_center
        self.bin_right_1d_from_2d = bin_right

    def combination_array(self):
        r"""pack the various pair combinations for each treatment into an array"""
        summary_treatment = {}
        for treatment in self.treatment_cases:
            pwr_treatment = np.zeros((self.num_k, self.num_comb))
            for comb, comb_index in zip(self.comb_cases, range(self.num_comb)):
                pwrcase = "%s:%s" % (comb, treatment)
                pwr_treatment[:, comb_index] = self.pwrspec_1d[pwrcase]
            summary_treatment[treatment] = pwr_treatment

        return summary_treatment

    def print_param(self):
        print self.parameters

    def summarize_agg_pwrspec(self, treatment, tag, outdir="./plot_data"):
        r"""Unclear what this will do here"""
        fileout = outdir + "/" + tag + "_avg_from2d.dat"
        corr_fileout = outdir + "/" + tag + "_corr_from2d.dat"
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
