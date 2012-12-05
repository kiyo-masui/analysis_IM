import numpy as np
import numpy.ma as ma
from quadratic_products import pwrspec_estimator as pe
from quadratic_products import power_spectrum as ps
import glob
import copy
from plotting import plot_slice
import shelve
import h5py
from kiyopy import parse_ini
from utils import file_tools
import shutil
# TODO: better interaction between mask and counts
# TODO: make object like AggregateSummary that recompiles the physical sims
# TODO: call AggregateStatistics on the physical sims, make plots
# TODO: make transfer from physical -> observed -> observed with beam/meansub

aggregatesummary_init = {
        "directory": "dir",
        "basefile": "file",
        "noiseweights_2dto1d": None,
        "fix_weight_treatment": None,
        "apply_2d_beamtransfer": None,
        "apply_2d_modetransfer": None,
        "subtract_pwrspec": None,
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

        data_subtract = None
        if self.params['subtract_pwrspec'] is not None:
            print "agg WARNING: you are subtracting a power spectrum"
            print "from file ", self.params['subtract_pwrspec']
            data_subtract = ps.PowerSpectrum(self.params['subtract_pwrspec'])

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

        # load the transfer functions and 2d->1d noise weights
        transfer_dict = pe.load_transferfunc(
                            self.params["apply_2d_beamtransfer"],
                            self.params["apply_2d_modetransfer"],
                            sim_toread.treatment_cases)

        weights_2d = None
        if self.params["noiseweights_2dto1d"] is not None:
            print "applying 2D noise weights: " + \
                self.params["noiseweights_2dto1d"]

            weightfile = h5py.File(self.params["noiseweights_2dto1d"], "r")
            weights_2d = {}
            for treatment in sim_toread.treatment_cases:
                if self.params["fix_weight_treatment"] is not None:
                    fixed_treatment = self.params["fix_weight_treatment"]
                    print "fixing weight for %s to value at %s" % \
                            (treatment, fixed_treatment)
                    weights_2d[treatment] = weightfile[fixed_treatment].value
                else:
                    weights_2d[treatment] = weightfile[treatment].value

            weightfile.close()

        result_dict = {}
        for treatment in sim_toread.treatment_cases:
            treatment_dict = {}
            treatment_dict["pk_1d"] = np.zeros((num_sim, num_k_1d))
            treatment_dict["pk_1d_from_2d"] = \
                        np.zeros((num_sim, num_k_1d_from_2d))

            treatment_dict["pkstd_1d"] = np.zeros((num_sim, num_k_1d))
            treatment_dict["pkstd_1d_from_2d"] = \
                        np.zeros((num_sim, num_k_1d_from_2d))

            treatment_dict["pk_2d"] = np.zeros((num_sim, num_kx, num_ky))

            treatment_dict["counts_1d"] = np.zeros((num_sim, num_k_1d))
            treatment_dict["counts_1d_from_2d"] = \
                        np.zeros((num_sim, num_k_1d_from_2d))

            treatment_dict["counts_2d"] = np.zeros((num_sim, num_kx, num_ky))
            result_dict[treatment] = treatment_dict

        for (simfile, index) in zip(filelist, range(num_sim)):
            print "processing ", simfile, index
            sim_toread = ps.PowerSpectrum(simfile)

            # optionally subtract some reference spectrum
            if data_subtract is not None:
                print "agg WARNING: you are subtracting a power spectrum"
                # assuming these both have the same treatments
                for pwrspec_case in sim_toread.pwrspec_1d:
                    sim_toread.pwrspec_1d[pwrspec_case] -= \
                        data_subtract.pwrspec_1d[pwrspec_case]

                    sim_toread.pwrspec_2d[pwrspec_case] -= \
                        data_subtract.pwrspec_2d[pwrspec_case]

            sim_toread.apply_2d_trans_by_treatment(transfer_dict)
            sim_toread.convert_2d_to_1d(weights_2d=weights_2d)

            agg1d = sim_toread.agg_stat_1d_pwrspec()
            agg1d_from_2d = sim_toread.agg_stat_1d_pwrspec(from_2d=True)
            agg2d = sim_toread.agg_stat_2d_pwrspec()

            for treatment in sim_toread.treatment_cases:
                result_dict[treatment]["pk_1d"][index, :] = \
                                            agg1d[treatment]['mean']

                result_dict[treatment]["pk_1d_from_2d"][index, :] = \
                                            agg1d_from_2d[treatment]['mean']

                result_dict[treatment]["pkstd_1d"][index, :] = \
                                            agg1d[treatment]['std']

                result_dict[treatment]["pkstd_1d_from_2d"][index, :] = \
                                            agg1d_from_2d[treatment]['std']

                result_dict[treatment]["pk_2d"][index, :, :] = \
                                            agg2d[treatment]['mean']

                result_dict[treatment]["counts_1d"][index, :] = \
                                            agg1d[treatment]['counts']

                result_dict[treatment]["counts_1d_from_2d"][index, :] = \
                                            agg1d_from_2d[treatment]['counts']

                result_dict[treatment]["counts_2d"][index, :, :] = \
                                            agg2d[treatment]['counts']

        # package all simulations of all treatments into a file
        #out_dicttree = shelve.open(outfile, "n", protocol=-1)
        out_dicttree = {}
        out_dicttree["k_1d"] = k_1d
        out_dicttree["k_1d_from_2d"] = k_1d_from_2d
        out_dicttree["kx_2d"] = kx_2d
        out_dicttree["ky_2d"] = ky_2d
        out_dicttree["results"] = result_dict
        file_tools.convert_numpytree_hdf5(out_dicttree, outfile)
        #out_dicttree.close()


aggregatestatistics_init = {
        "aggfile_in": "file.hd5",
        "statfile_out": "data.hd5",
        "outplotdir": "dir"
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

        print "opening: ", self.params["aggfile_in"]
        self.summary = h5py.File(self.params["aggfile_in"], "r")
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

            stat_summary["k_1d"] = self.summary["k_1d"]["center"].value
            stat_summary["k_1d_from_2d"] = \
                    self.summary["k_1d_from_2d"]["center"].value

            stat_summary["kx_2d"] = self.summary["kx_2d"]["center"].value
            stat_summary["ky_2d"] = self.summary["ky_2d"]["center"].value

            stat_results[treatment] = copy.deepcopy(stat_summary)

        # in shelve file model
        #self.summary["stats"] = stat_results
        #print self.summary["stats"]["0modes"].keys()

        self.summary.close()
        file_tools.convert_numpytree_hdf5(stat_results,
                                          self.params["statfile_out"])

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
            trial_array_1d = results["pk_1d_from_2d"].value
            error_array_1d = results["pkstd_1d_from_2d"].value
            counts_array_1d = results["counts_1d_from_2d"].value

            outfile = "%s/power_1d_from_2d_%s.dat" % \
                      (self.params['outplotdir'], treatment)
        else:
            k_1d = self.summary["k_1d"]
            trial_array_1d = results["pk_1d"].value
            error_array_1d = results["pkstd_1d"].value
            counts_array_1d = results["counts_1d"].value

            outfile = "%s/power_1d_%s.dat" % \
                      (self.params['outplotdir'], treatment)

        stat_1d = self.calc_stat_1d(trial_array_1d)
        error_stat_1d = self.calc_stat_1d(error_array_1d)
        counts_1d = self.calc_stat_1d(counts_array_1d)

        file_tools.print_multicolumn(k_1d["left"].value,
                                     k_1d["center"].value,
                                     k_1d["right"].value,
                                     counts_1d["mean"],
                                     stat_1d["mean"],
                                     stat_1d["std"],
                                     error_stat_1d["mean"],
                                     outfile=outfile)

        logk_1d = np.log10(k_1d['center'])
        if self.make_plot:
            outplot_file = "%s/sim_corr_1d_%s.png" % \
                           (self.params['outplotdir'], treatment)
            plot_slice.simpleplot_2D(outplot_file, stat_1d['corr'],
                                     logk_1d, logk_1d,
                                     ["logk", "logk"], 1.,
                                     "1D corr", "corr")

            outplot_file = "%s/sim_cov_1d_%s.png" % \
                           (self.params['outplotdir'], treatment)
            plot_slice.simpleplot_2D(outplot_file, stat_1d['cov'],
                                     logk_1d, logk_1d,
                                     ["logk", "logk"], 1.,
                                     "1D covariance", "cov")

        return (stat_1d, error_stat_1d, counts_1d)

    def calc_stat_2d(self, trial_array, calc_corr=False):
        stat_2d = {}
        old_settings = np.seterr(invalid="warn", under="warn")
        mtrial_array = ma.masked_invalid(trial_array)

        hcheck = []
        for index in range(mtrial_array.shape[0]):
            hcheck.append(np.ma.sum(mtrial_array[index, :, :]))

        hcheck = np.array(hcheck)
        meanval = np.mean(hcheck)
        hcheck = (hcheck - meanval) / meanval
        if not np.all(hcheck == 0.):
            print hcheck

        #print np.apply_over_axes(np.ma.sum, mtrial_array, axes=[1,2])

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
                  (self.params['outplotdir'], treatment)

        stat_2d = self.calc_stat_2d(trial_array_2d)
        counts_2d = self.calc_stat_2d(counts_array_2d)

        #mean2d[np.isnan(mean2d)] = 0.
        #mean2d[np.isnan(mean2d)] = 1.e-16
        np.set_printoptions(threshold='nan')

        logkx = np.log10(kx_2d['center'])
        logky = np.log10(ky_2d['center'])

        if self.make_plot:
            outplot_file = "%s/sim_mean_2d_%s.png" % \
                      (self.params['outplotdir'], treatment)
            plot_slice.simpleplot_2D(outplot_file, stat_2d['mean'],
                                     logkx, logky,
                                     ["logkx", "logky"], 1.,
                                     "2D power", "logP(k)",
                                     logscale=False)

            # can use C or F to do column or row-major
            #outplot_file = "%s/sim_corr_2d_%s.png" % \
            #          (self.params['outplotdir'], treatment)
            #plot_slice.simpleplot_2D(outplot_file, stat_2d['corr'],
            #                         stat_2d['flat_axis'],
            #                         stat_2d['flat_axis'],
            #                         ["k", "k"], 1.,
            #                         "2D power corr", "corr")

        return (stat_2d, counts_2d)


calculatetransfer_init = {
        "powerfile_in": "file",
        "powerfile_out": "file",
        "transferfile": "file",
        "input_multiplier": 1.,
        "outplotdir": "dir",
    }
calculatetransfer_prefix = 'atr_'


class CalculateTransfer(object):
    """Calculate a transfer function in 2d
    """

    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        np.seterr(invalid='raise')

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          calculatetransfer_init,
                                          prefix=calculatetransfer_prefix)

        print self.params["powerfile_in"], "->", self.params["powerfile_out"]

        self.stats_in = h5py.File(self.params["powerfile_in"], "r")
        self.stats_out = h5py.File(self.params["powerfile_out"], "r")
        self.treatments_in = self.stats_in.keys()
        self.treatments_out = self.stats_out.keys()
        self.multiplier = self.params["input_multiplier"]

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
        transferfile = h5py.File(self.params["transferfile"], "w")
        inkey = self.treatments_in[0]
        for treatment in self.treatments_out:

            stats_2d_in = self.stats_in[inkey]["pk_2d_stat"]["mean"].value
            stats_2d_out = self.stats_out[treatment]["pk_2d_stat"]["mean"].value
            counts_2d_in = self.stats_in[inkey]["pk_2d_counts"]["mean"].value
            counts_2d_out = self.stats_out[treatment]["pk_2d_counts"]["mean"].value

            counts_prod = counts_2d_in * counts_2d_out
            transfer_2d = stats_2d_out / (stats_2d_in * self.multiplier)
            transfer_2d[counts_prod == 0] = 0.
            transferfile[treatment] = transfer_2d

            transfer_2d_plot = copy.deepcopy(transfer_2d)
            transfer_2d_plot[transfer_2d_plot < 0.] = 0.
            transfer_2d_plot[transfer_2d_plot > 1.] = 1.

            outplot_file = "%s/transfer_2d_%s.png" % \
                           (self.params['outplotdir'], treatment)

            print "plotting trans to ", outplot_file
            logkx = np.log10(self.stats_in[inkey]["kx_2d"].value)
            logky = np.log10(self.stats_in[inkey]["ky_2d"].value)
            plot_slice.simpleplot_2D(outplot_file, transfer_2d_plot,
                                     logkx, logky,
                                     ["logkx", "logky"], 1.,
                                     "2D beam transfer", "T")

        transferfile.close()
        self.stats_in.close()
        self.stats_out.close()


