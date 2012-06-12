import numpy as np
import numpy.ma as ma
from quadratic_products import power_spectrum as ps
import glob
from plotting import plot_slice
import shelve
from kiyopy import parse_ini
from utils import file_tools
# TODO: better interaction between mask and counts

aggregatesummary_init = {
        "directory": "dir",
        "basefile": "file",
        "apply_2d_transfer": None,
        "outfile": "file"
    }

aggregatesummary_prefix = 'as_'

class AggregateSummary(object):
    """take a directory of simulation outputs and merge them into a single
    numpy array/shelve file
    TODO: extend this to include multiple treatments (mode subtraction, etc.)
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          aggregatesummary_init,
                                          prefix=aggregatesummary_prefix)

    def execute(self, processes):
        """produce a list of files to combine and run"""
        outfile = self.params["outfile"]
        filelist = glob.glob('%s/%s*.shelve' % \
                             (self.params['directory'],
                              self.params['basefile']))

        self.produce_summary(filelist, outfile)

    def produce_summary(self, filelist, outfile, debug=False):
        num_sim = len(filelist)
        print "number of simulations to aggregate: %d" % num_sim

        sim_toread = ps.PowerSpectrum(filelist[0])
        # run this here just to get the 2D->1D bin centers
        sim_toread.convert_2d_to_1d()

        k_1d = sim_toread.k_1d
        k_1d_from_2d = sim_toread.k_1d_from_2d
        num_k_1d = sim_toread.num_k_1d
        num_k_1d_from_2d = sim_toread.num_k_1d_from_2d

        (kx_2d, ky_2d) = (sim_toread.kx_2d, sim_toread.ky_2d)
        (num_kx, num_ky) = (sim_toread.num_kx, sim_toread.num_ky)

        print "k_1d bins: %d, %d, kx bins: %d, ky bins: %d" % \
               (num_k_1d, num_k_1d_from_2d, num_kx, num_ky)

        trial_array_1d = np.zeros((num_sim, num_k_1d))
        trial_array_1d_from_2d = np.zeros((num_sim, num_k_1d_from_2d))
        trial_array_2d = np.zeros((num_sim, num_kx, num_ky))
        counts_array_1d = np.zeros((num_sim, num_k_1d))
        counts_array_1d_from_2d = np.zeros((num_sim, num_k_1d_from_2d))
        counts_array_2d = np.zeros((num_sim, num_kx, num_ky))

        if self.params["apply_2d_transfer"] is not None:
            trans_shelve = shelve.open(self.params["apply_2d_transfer"])
            transfer_2d = trans_shelve["transfer_2d"]
            trans_shelve.close()

            transfer_dict = {}
            for treatment in sim_toread.treatment_cases:
                transfer_dict[treatment] = transfer_2d

        for (simfile, index) in zip(filelist, range(num_sim)):
            print simfile, index
            sim_toread = ps.PowerSpectrum(simfile)

            if self.params["apply_2d_transfer"] is not None:
                sim_toread.apply_2d_trans_by_treatment(transfer_dict)

            sim_toread.convert_2d_to_1d()

            agg1d = sim_toread.agg_stat_1d_pwrspec()
            agg1d_from_2d = sim_toread.agg_stat_1d_pwrspec(from_2d=True)
            agg2d = sim_toread.agg_stat_2d_pwrspec()

            trial_array_1d[index, :] = agg1d['0modes']['mean']
            trial_array_1d_from_2d[index, :] = agg1d_from_2d['0modes']['mean']
            trial_array_2d[index, :, :] = agg2d['0modes']['mean']
            counts_array_1d[index, :] = agg1d['0modes']['counts']
            counts_array_1d_from_2d[index, :] = agg1d_from_2d['0modes']['counts']
            counts_array_2d[index, :, :] = agg2d['0modes']['counts']
            if debug:
                nanarray = np.isnan(trial_array_2d[index, :, :])
                print "is nan", np.sum(nanarray)

            del sim_toread

        outshelve = shelve.open(outfile, "n", protocol=-1)
        outshelve["k_1d"] = k_1d
        outshelve["k_1d_from_2d"] = k_1d_from_2d
        outshelve["kx_2d"] = kx_2d
        outshelve["ky_2d"] = ky_2d
        outshelve["pk_1d"] = trial_array_1d
        outshelve["pk_1d_from_2d"] = trial_array_1d_from_2d
        outshelve["pk_2d"] = trial_array_2d
        outshelve["counts_1d"] = counts_array_1d
        outshelve["counts_1d_from_2d"] = counts_array_1d_from_2d
        outshelve["counts_2d"] = counts_array_2d
        outshelve.close()

aggregatestatistics_init = {
        "shelvefile": "file",
        "outputdir": "dir",
    }

aggregatestatistics_prefix = 'ast_'

class AggregateStatistics(object):
    """take the summary shelve and find statistics on it
    TODO: extend this to include multiple treatments (mode subtraction, etc.)
    TODO: have this write plots out to uniform directories
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          aggregatestatistics_init,
                                          prefix=aggregatestatistics_prefix)

        self.summary = shelve.open(self.params["shelvefile"])

    def execute(self, processes):
        """produce a list of files to combine and run"""
        self.aggregate_1d_statistics(from2d = False)
        self.aggregate_1d_statistics(from2d = True)
        self.aggregate_2d_statistics()
        self.summary.close()

    def calc_stat_1d(self, trial_array):
        """take an 1D array of power spectra and find some basic statistics
        on it"""
        stat_1d = {}

        mtrial_array = ma.masked_invalid(trial_array)
        stat_1d['mean'] = np.ma.mean(mtrial_array, axis=0)
        stat_1d['std'] = np.ma.std(mtrial_array, axis=0, ddof=1)
        stat_1d['corr'] = np.ma.corrcoef(np.ma.transpose(mtrial_array), ddof=1)
        stat_1d['cov'] = np.ma.cov(np.ma.transpose(mtrial_array), ddof=1)

        stat_1d['mean'] = np.ma.filled(stat_1d['mean'], fill_value=np.nan)
        stat_1d['std'] = np.ma.filled(stat_1d['std'], fill_value=np.nan)

        return stat_1d

    def aggregate_1d_statistics(self, from2d = False):
        r"""Find the mean, std, cov, corr etc of 1D p(k)'s"""

        if from2d:
            k_1d = self.summary["k_1d_from_2d"]
            trial_array_1d = self.summary["pk_1d_from_2d"]
            counts_array_1d = self.summary["counts_1d_from_2d"]
            outfile = self.params['outputdir'] + "power_1d_from_2d.dat"
        else:
            k_1d = self.summary["k_1d"]
            trial_array_1d = self.summary["pk_1d"]
            counts_array_1d = self.summary["counts_1d"]
            outfile = self.params['outputdir'] + "power_1d.dat"

        stat_1d = self.calc_stat_1d(trial_array_1d)
        counts_1d = self.calc_stat_1d(counts_array_1d)

        if from2d:
            self.summary["pk_1d_from_2d_stat"] = stat_1d
            self.summary["pk_1d_from_2d_counts"] = counts_1d
        else:
            self.summary["pk_1d_stat"] = stat_1d
            self.summary["pk_1d_counts"] = counts_1d

        file_tools.print_multicolumn(k_1d["left"],
                                     k_1d["center"],
                                     k_1d["right"],
                                     counts_1d["mean"],
                                     stat_1d["mean"],
                                     stat_1d["std"],
                                     outfile = outfile)

        logk_1d = np.log10(k_1d['center'])
        outplot_file = self.params['outputdir'] + "sim_corr_1d.png"
        plot_slice.simpleplot_2D(outplot_file, stat_1d['corr'], logk_1d, logk_1d,
                             ["logk", "logk"], 1., "1D corr", "corr")

        outplot_file = self.params['outputdir'] + "sim_cov_1d.png"
        plot_slice.simpleplot_2D(outplot_file, stat_1d['cov'], logk_1d, logk_1d,
                             ["logk", "logk"], 1., "1D covariance", "cov")

    def calc_stat_2d(self, trial_array):
        stat_2d = {}
        mtrial_array = ma.masked_invalid(trial_array)

        stat_2d['mean'] = np.ma.mean(mtrial_array, axis=0)
        stat_2d['std'] = np.ma.std(mtrial_array, axis=0, ddof=1)
        stat_2d['mean'] = np.ma.filled(stat_2d['mean'], fill_value=np.nan)
        stat_2d['std'] = np.ma.filled(stat_2d['std'], fill_value=np.nan)

        num_flat = trial_array.shape[1] * trial_array.shape[2]
        num_sim = trial_array.shape[0]
        stat_2d['flat_axis'] = range(num_flat)
        mtrial_array_flat = np.ma.reshape(mtrial_array, (num_sim, num_flat))
        stat_2d['corr'] = np.corrcoef(np.ma.transpose(mtrial_array_flat), ddof=1)
        stat_2d['cov'] = np.ma.cov(np.ma.transpose(mtrial_array_flat), ddof=1)

        return stat_2d

    def aggregate_2d_statistics(self):
        kx_2d = self.summary["kx_2d"]
        ky_2d = self.summary["ky_2d"]
        trial_array_2d = self.summary["pk_2d"]
        counts_array_2d = self.summary["counts_2d"]
        outfile = self.params['outputdir'] + "power_2d.dat"

        stat_2d = self.calc_stat_2d(trial_array_2d)
        counts_2d = self.calc_stat_2d(counts_array_2d)
        self.summary["pk_2d_stat"] = stat_2d
        self.summary["pk_2d_counts"] = counts_2d

        #mean2d[np.isnan(mean2d)] = 0.
        #mean2d[np.isnan(mean2d)] = 1.e-16
        np.set_printoptions(threshold='nan')

        logkx = np.log10(kx_2d['center'])
        logky = np.log10(ky_2d['center'])

        outplot_file = self.params['outputdir'] + "sim_mean_2d.png"
        plot_slice.simpleplot_2D(outplot_file, stat_2d['mean'], logkx, logky,
                                 ["logkx", "logky"], 1., "2D power", "logP(k)",
                                 logscale=False)

        # can use C or F to do column or row-major
        outplot_file = self.params['outputdir'] + "sim_corr_2d.png"
        plot_slice.simpleplot_2D(outplot_file, stat_2d['corr'],
                                 stat_2d['flat_axis'], stat_2d['flat_axis'],
                                 ["k", "k"], 1., "2D power corr", "corr")

calculatetransfer_init = {
        "shelvefile_in": "file",
        "shelvefile_out": "file",
        "transferfile": "file",
        "outputdir": "dir",
    }

calculatetransfer_prefix = 'atr_'

class CalculateTransfer(object):
    """Calculate a transfer function in 1d/2d
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          calculatetransfer_init,
                                          prefix=calculatetransfer_prefix)

        self.stats_in = shelve.open(self.params["shelvefile_in"])
        self.stats_out = shelve.open(self.params["shelvefile_out"])

    def execute(self, processes):
        """produce a list of files to combine and run"""
        stats_1d_in = self.stats_in["pk_1d_stat"]["mean"]
        stats_1d_out = self.stats_out["pk_1d_stat"]["mean"]
        stats_2d_in = self.stats_in["pk_2d_stat"]["mean"]
        stats_2d_out = self.stats_out["pk_2d_stat"]["mean"]

        counts_1d_in = self.stats_in["pk_1d_counts"]["mean"]
        counts_1d_out = self.stats_out["pk_1d_counts"]["mean"]
        counts_2d_in = self.stats_in["pk_2d_counts"]["mean"]
        counts_2d_out = self.stats_out["pk_2d_counts"]["mean"]

        counts_prod = counts_2d_in * counts_2d_out
        transfer_2d = stats_2d_out / stats_2d_in
        transfer_2d[counts_prod == 0] = 0.

        k_1d = self.stats_in["k_1d"]
        k_1d_from_2d = self.stats_in["k_1d_from_2d"]
        kx_2d = self.stats_in["kx_2d"]
        ky_2d = self.stats_in["ky_2d"]

        outplot_file = self.params['outputdir'] + "transfer_2d.png"
        logkx = np.log10(kx_2d['center'])
        logky = np.log10(ky_2d['center'])
        plot_slice.simpleplot_2D(outplot_file, transfer_2d,
                                 logkx, logky,
                                 ["logkx", "logky"], 1., "2D beam transfer", "T")

        transferfile = shelve.open(self.params["transferfile"])
        transferfile["transfer_2d"] = transfer_2d
        transferfile.close()

        self.stats_in.close()
        self.stats_out.close()
