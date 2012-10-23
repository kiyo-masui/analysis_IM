import numpy as np
import numpy.ma as ma
from quadratic_products import power_spectrum as ps
import glob
import copy
from plotting import plot_slice
import shelve
import h5py
from kiyopy import parse_ini
from utils import file_tools
# TODO: better interaction between mask and counts
# TODO: make object like AggregateSummary that recompiles the physical sims
# TODO: call AggregateStatistics on the physical sims, make plots
# TODO: make transfer from physical -> observed -> observed with beam/meansub

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
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        np.seterr(invalid='raise')

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

        print '%s/%s*.shelve' % \
              (self.params['directory'], self.params['basefile'])

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

        print "AggregateSummary: treatment cases: ", sim_toread.treatment_cases

        if self.params["apply_2d_transfer"] is not None:
            trans_shelve = shelve.open(self.params["apply_2d_transfer"])
            transfer_2d = trans_shelve["transfer_2d"]
            trans_shelve.close()

            transfer_dict = {}
            for treatment in sim_toread.treatment_cases:
                transfer_dict[treatment] = transfer_2d

        result_dict = {}
        for treatment in sim_toread.treatment_cases:
            treatment_dict = {}
            trial_array_1d = np.zeros((num_sim, num_k_1d))
            trial_array_1d_from_2d = np.zeros((num_sim, num_k_1d_from_2d))
            trial_array_2d = np.zeros((num_sim, num_kx, num_ky))
            error_array_1d = np.zeros((num_sim, num_k_1d))
            error_array_1d_from_2d = np.zeros((num_sim, num_k_1d_from_2d))
            counts_array_1d = np.zeros((num_sim, num_k_1d))
            counts_array_1d_from_2d = np.zeros((num_sim, num_k_1d_from_2d))
            counts_array_2d = np.zeros((num_sim, num_kx, num_ky))
            for (simfile, index) in zip(filelist, range(num_sim)):
                print treatment, simfile, index
                sim_toread = ps.PowerSpectrum(simfile)

                if self.params["apply_2d_transfer"] is not None:
                    sim_toread.apply_2d_trans_by_treatment(transfer_dict)

                sim_toread.convert_2d_to_1d()

                agg1d = sim_toread.agg_stat_1d_pwrspec()
                agg1d_from_2d = sim_toread.agg_stat_1d_pwrspec(from_2d=True)
                agg2d = sim_toread.agg_stat_2d_pwrspec()

                # accumulate the means
                trial_array_1d[index, :] = agg1d[treatment]['mean']

                trial_array_1d_from_2d[index, :] = \
                                            agg1d_from_2d[treatment]['mean']

                # accumulate the std across 6 pairs
                error_array_1d[index, :] = agg1d[treatment]['std']

                error_array_1d_from_2d[index, :] = \
                                            agg1d_from_2d[treatment]['std']

                # accumulate the 2D powers
                trial_array_2d[index, :, :] = agg2d[treatment]['mean']

                # accumulate the counts arrays
                counts_array_1d[index, :] = agg1d[treatment]['counts']

                counts_array_1d_from_2d[index, :] = \
                                            agg1d_from_2d[treatment]['counts']

                counts_array_2d[index, :, :] = agg2d[treatment]['counts']

                if debug:
                    nanarray = np.isnan(trial_array_2d[index, :, :])
                    print "is nan", np.sum(nanarray)

                del sim_toread

            # package all the simulations of a given treatment
            treatment_dict["pk_1d"] = copy.deepcopy(trial_array_1d)
            treatment_dict["pk_1d_from_2d"] = \
                            copy.deepcopy(trial_array_1d_from_2d)

            treatment_dict["pkstd_1d"] = copy.deepcopy(error_array_1d)
            treatment_dict["pkstd_1d_from_2d"] = \
                            copy.deepcopy(error_array_1d_from_2d)

            treatment_dict["pk_2d"] = copy.deepcopy(trial_array_2d)

            treatment_dict["counts_1d"] = copy.deepcopy(counts_array_1d)
            treatment_dict["counts_1d_from_2d"] = \
                            copy.deepcopy(counts_array_1d_from_2d)

            treatment_dict["counts_2d"] = copy.deepcopy(counts_array_2d)
            result_dict[treatment] = treatment_dict

        # package all simulations of all treatments into a file
        outshelve = shelve.open(outfile, "n", protocol=-1)
        outshelve["k_1d"] = k_1d
        outshelve["k_1d_from_2d"] = k_1d_from_2d
        outshelve["kx_2d"] = kx_2d
        outshelve["ky_2d"] = ky_2d
        outshelve["results"] = result_dict
        outshelve.close()


aggregatestatistics_init = {
        "shelvefile": "file",
        "outputdir": "dir",
    }

aggregatestatistics_prefix = 'ast_'


class AggregateStatistics(object):
    """take the summary shelve and find statistics on it
    this extends the shelve object assembled in AggregateSummary
    TODO: have this write plots out to uniform directories
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0,
                 make_plot=True):
        self.params = params_dict
        np.seterr(under='raise')
        self.make_plot = make_plot

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          aggregatestatistics_init,
                                          prefix=aggregatestatistics_prefix)

        self.summary = shelve.open(self.params["shelvefile"])
        # get the list of treatments
        self.treatments = self.summary["results"].keys()
        print "AggregateStatistics: treatment cases: ", self.treatments

    def execute(self, processes):
        """produce a list of files to combine and run"""
        stat_results = {}
        for treatment in self.treatments:
            stat_summary = {}
            (stat_1d, error_stat_1d, counts_1d) = \
                    self.aggregate_1d_statistics(treatment, from2d=False)

            stat_summary["pk_1d_stat"] = stat_1d
            stat_summary["pk_1d_errstat"] = error_stat_1d
            stat_summary["pk_1d_counts"] = counts_1d

            (stat_1d, error_stat_1d, counts_1d) = \
                    self.aggregate_1d_statistics(treatment, from2d=True)

            stat_summary["pk_1d_from_2d_stat"] = stat_1d
            stat_summary["pk_1d_from_2d_errstat"] = error_stat_1d
            stat_summary["pk_1d_from_2d_counts"] = counts_1d

            (stat_2d, counts_2d) = self.aggregate_2d_statistics(treatment)

            stat_summary["pk_2d_stat"] = stat_2d
            stat_summary["pk_2d_counts"] = counts_2d

            stat_results[treatment] = copy.deepcopy(stat_summary)

        self.summary["stats"] = stat_results
        print self.summary["stats"]["0modes"].keys()

        self.summary.close()

    def calc_stat_1d(self, trial_array, calc_corr=True):
        """take an 1D array of power spectra and find some basic statistics
        on it"""
        stat_1d = {}

        old_settings = np.seterr(invalid="warn", under="warn")
        mtrial_array = ma.masked_invalid(trial_array)

        stat_1d['mean'] = np.ma.filled(np.ma.mean(mtrial_array, axis=0),
                                       fill_value=np.nan)

        stat_1d['std'] = np.ma.filled(np.ma.std(mtrial_array, axis=0, ddof=1),
                                      fill_value=np.nan)

        if calc_corr:
            mtrial_array_trans = np.ma.transpose(
                                        ma.masked_invalid(trial_array))

            stat_1d['corr'] = np.ma.filled(np.ma.corrcoef(mtrial_array_trans,
                                           ddof=1), fill_value=np.nan)

            stat_1d['cov'] = np.ma.filled(np.ma.cov(mtrial_array_trans,
                                          ddof=1), fill_value=np.nan)

        np.seterr(**old_settings)

        return stat_1d

    def aggregate_1d_statistics(self, treatment, from2d=False):
        r"""Find the mean, std, cov, corr etc of 1D p(k)'s"""
        results = self.summary["results"][treatment]

        if from2d:
            k_1d = self.summary["k_1d_from_2d"]
            trial_array_1d = results["pk_1d_from_2d"]
            error_array_1d = results["pkstd_1d_from_2d"]
            counts_array_1d = results["counts_1d_from_2d"]

            outfile = "%s/power_1d_from_2d_%s.dat" % \
                      (self.params['outputdir'], treatment)
        else:
            k_1d = self.summary["k_1d"]
            trial_array_1d = results["pk_1d"]
            error_array_1d = results["pkstd_1d"]
            counts_array_1d = results["counts_1d"]

            outfile = "%s/power_1d_%s.dat" % \
                      (self.params['outputdir'], treatment)

        stat_1d = self.calc_stat_1d(trial_array_1d)
        error_stat_1d = self.calc_stat_1d(error_array_1d)
        counts_1d = self.calc_stat_1d(counts_array_1d)

        file_tools.print_multicolumn(k_1d["left"],
                                     k_1d["center"],
                                     k_1d["right"],
                                     counts_1d["mean"],
                                     stat_1d["mean"],
                                     stat_1d["std"],
                                     error_stat_1d["mean"],
                                     outfile=outfile)

        logk_1d = np.log10(k_1d['center'])
        if self.make_plot:
            outplot_file = "%s/sim_corr_1d_%s.png" % \
                           (self.params['outputdir'], treatment)
            plot_slice.simpleplot_2D(outplot_file, stat_1d['corr'],
                                     logk_1d, logk_1d,
                                     ["logk", "logk"], 1.,
                                     "1D corr", "corr")

            outplot_file = "%s/sim_cov_1d_%s.png" % \
                           (self.params['outputdir'], treatment)
            plot_slice.simpleplot_2D(outplot_file, stat_1d['cov'],
                                     logk_1d, logk_1d,
                                     ["logk", "logk"], 1.,
                                     "1D covariance", "cov")

        return (stat_1d, error_stat_1d, counts_1d)

    def calc_stat_2d(self, trial_array, calc_corr=False):
        stat_2d = {}
        old_settings = np.seterr(invalid="warn", under="warn")
        mtrial_array = ma.masked_invalid(trial_array)

        stat_2d['mean'] = np.ma.filled(np.ma.mean(mtrial_array, axis=0),
                                       fill_value=np.nan)

        stat_2d['std'] = np.ma.filled(np.ma.std(mtrial_array, axis=0, ddof=1),
                                      fill_value=np.nan)

        num_flat = trial_array.shape[1] * trial_array.shape[2]
        num_sim = trial_array.shape[0]
        stat_2d['flat_axis'] = range(num_flat)

        if calc_corr:
            mtrial_array_flat = np.ma.transpose(np.ma.reshape(mtrial_array,
                                                (num_sim, num_flat)))

            stat_2d['corr'] = np.ma.filled(np.ma.corrcoef(mtrial_array_flat,
                                           ddof=1), fill_value=np.nan)

            stat_2d['cov'] = np.ma.filled(np.ma.cov(mtrial_array_flat,
                                           ddof=1, fill_value=np.nan))

        np.seterr(**old_settings)
        return stat_2d

    def aggregate_2d_statistics(self, treatment):
        kx_2d = self.summary["kx_2d"]
        ky_2d = self.summary["ky_2d"]
        trial_array_2d = self.summary["results"][treatment]["pk_2d"]
        counts_array_2d = self.summary["results"][treatment]["counts_2d"]
        outfile = "%s/power_2d_%s.dat" % \
                  (self.params['outputdir'], treatment)

        stat_2d = self.calc_stat_2d(trial_array_2d)
        counts_2d = self.calc_stat_2d(counts_array_2d)

        #mean2d[np.isnan(mean2d)] = 0.
        #mean2d[np.isnan(mean2d)] = 1.e-16
        np.set_printoptions(threshold='nan')

        logkx = np.log10(kx_2d['center'])
        logky = np.log10(ky_2d['center'])

        if self.make_plot:
            outplot_file = "%s/sim_mean_2d_%s.png" % \
                      (self.params['outputdir'], treatment)
            plot_slice.simpleplot_2D(outplot_file, stat_2d['mean'],
                                     logkx, logky,
                                     ["logkx", "logky"], 1.,
                                     "2D power", "logP(k)",
                                     logscale=False)

            # can use C or F to do column or row-major
            #outplot_file = "%s/sim_corr_2d_%s.png" % \
            #          (self.params['outputdir'], treatment)
            #plot_slice.simpleplot_2D(outplot_file, stat_2d['corr'],
            #                         stat_2d['flat_axis'],
            #                         stat_2d['flat_axis'],
            #                         ["k", "k"], 1.,
            #                         "2D power corr", "corr")

        return (stat_2d, counts_2d)


calculatetransfer_init = {
        "shelvefile_in": "file",
        "shelvefile_out": "file",
        "shelvedatalike_out": "file",
        "transferfile": "file",
        "outputdir": "dir",
    }

calculatetransfer_prefix = 'atr_'


class CalculateTransfer(object):
    """Calculate a transfer function in 1d/2d
    The shelve files in can either be those produced by AggregateSummary above,
    (in the case that they are many sims for each
    TODO: add 1D transfer if we ever need it
    """

    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        np.seterr(invalid='raise')

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          calculatetransfer_init,
                                          prefix=calculatetransfer_prefix)

        print self.params["shelvefile_in"], "->", self.params["shelvefile_out"]

        print "saving sims to data-like file: ", \
              self.params["shelvedatalike_out"]

        self.stats_in = shelve.open(self.params["shelvefile_in"], "r")
        self.stats_dataout = shelve.open(self.params["shelvedatalike_out"])
        self.stats_out = shelve.open(self.params["shelvefile_out"], "w")
        self.treatments_in = self.stats_in["results"].keys()
        self.treatments_out = self.stats_out["results"].keys()

        # If this is measuring the mode cleaning transfer function, it will be
        # with respect to the 0modes removed case. Note that map+sim (with zero
        # modes removed) x sim is quite noisy, so we really don't want to use
        # the zero-mode case of the plussim cleaning runs to estimate this.
        # Instead, use signal-only sims for the 0mode reference.
        print "AggregateStatistics: input treatments: ", self.treatments_in
        print "AggregateStatistics: output treatments: ", self.treatments_out

        if self.treatments_in[0] != "0modes":
            print "Transfer functions must be wrt only 0modes"
            return

    def execute(self, processes):
        """produce a list of files to combine and run"""
        transfer_compilation = {}
        for treatment in self.treatments_out:
            print treatment
            print self.stats_out["results"][treatment].keys()

            stat_in = self.stats_in["stats"]
            stat_out = self.stats_out["stats"]

            stats_1d_in = stat_in["0modes"]["pk_1d_stat"]["mean"]
            stats_1d_out = stat_out[treatment]["pk_1d_stat"]["mean"]
            stats_2d_in = stat_in["0modes"]["pk_2d_stat"]["mean"]
            stats_2d_out = stat_out[treatment]["pk_2d_stat"]["mean"]

            counts_1d_in = stat_in["0modes"]["pk_1d_counts"]["mean"]
            counts_1d_out = stat_out[treatment]["pk_1d_counts"]["mean"]
            counts_2d_in = stat_in["0modes"]["pk_2d_counts"]["mean"]
            counts_2d_out = stat_out[treatment]["pk_2d_counts"]["mean"]

            counts_prod = counts_2d_in * counts_2d_out
            transfer_2d = stats_2d_out / stats_2d_in
            transfer_2d[counts_prod == 0] = 0.

            k_1d = self.stats_in["k_1d"]
            k_1d_from_2d = self.stats_in["k_1d_from_2d"]
            kx_2d = self.stats_in["kx_2d"]
            ky_2d = self.stats_in["ky_2d"]

            transfer_2d_plot = copy.deepcopy(transfer_2d)
            transfer_2d_plot[transfer_2d_plot < 0.] = 0.
            transfer_2d_plot[transfer_2d_plot > 1.] = 1.

            outplot_file = "%s/transfer_2d_%s.png" % \
                           (self.params['outputdir'], treatment)

            logkx = np.log10(kx_2d['center'])
            logky = np.log10(ky_2d['center'])

            plot_slice.simpleplot_2D(outplot_file, transfer_2d_plot,
                                     logkx, logky,
                                     ["logkx", "logky"], 1.,
                                     "2D beam transfer", "T")

            # make a package that looks like the data
            pwrspec2d_product = {}
            pwrspec2d_product["bin_x_left"] = kx_2d["left"]
            pwrspec2d_product["bin_x_center"] = kx_2d["center"]
            pwrspec2d_product["bin_x_right"] = kx_2d["right"]
            pwrspec2d_product["bin_y_left"] = ky_2d["left"]
            pwrspec2d_product["bin_y_center"] = ky_2d["center"]
            pwrspec2d_product["bin_y_right"] = ky_2d["right"]
            pwrspec2d_product["counts_histo"] = counts_2d_out
            pwrspec2d_product["binavg"] = stats_2d_out

            pwrspec1d_product = {}
            pwrspec1d_product["bin_left"] = k_1d["left"]
            pwrspec1d_product["bin_center"] = k_1d["center"]
            pwrspec1d_product["bin_right"] = k_1d["right"]
            pwrspec1d_product["counts_histo"] = counts_1d_out
            pwrspec1d_product["binavg"] = stats_1d_out

            datakey = "data:%s" % treatment
            self.stats_dataout[datakey] = (0, (pwrspec2d_product,
                                               pwrspec1d_product))

            transfer_compilation[treatment] = transfer_2d

        #transferfile = shelve.open(self.params["transferfile"], protocol=0)
        #transferfile["transfer_2d"] = transfer_compilation

        transferfile = h5py.File(self.params["transferfile"], "w")
        for treatment in self.treatments_out:
            transferfile[treatment] = transfer_compilation[treatment]

        transferfile.close()

        self.stats_in.close()
        self.stats_dataout.close()
        self.stats_out.close()
