from plotting import plot_slice
import numpy as np
from utils import file_tools
from quadratic_products import pwrspec_estimator as pe
from quadratic_products import power_spectrum as ps
from kiyopy import parse_ini
import shelve
import h5py


analyzeautopowerparams_init = {
        "sim_auto_summary": "sim_auto_summary.shelve",
        "sim_xspec_summary": "sim_auto_summary.shelve",
        "data_auto_summary": "data_auto_summary.hd5",
        "data_xspec_summary": "data_auto_summary.hd5",
        "outdir": "./"
               }
analyzeautopowerprefix = 'analyzeautopower_'


class AnalyzeAutopower(object):
    r"""Do some set of final analyses on the autopower"""

    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          analyzeautopowerparams_init,
                                          prefix=analyzeautopowerprefix)

        print self.params
        self.data_auto = h5py.File(self.params["data_auto_summary"], "r")
        self.data_xspec = h5py.File(self.params["data_xspec_summary"], "r")
        self.sim_auto = shelve.open(self.params["sim_auto_summary"], "r")
        self.sim_xspec = shelve.open(self.params["sim_xspec_summary"], "r")

    def execute(self, processes):
        sim_auto_stats = self.sim_auto['stats']['0modes']['pk_1d_from_2d_stat']
        sim_xspec_stats = self.sim_xspec['stats']['0modes']['pk_1d_from_2d_stat']

        sim_xspec_mean = sim_xspec_stats['mean']
        sim_xspec_cov = sim_xspec_stats['cov']
        sim_auto_mean = sim_auto_stats['mean']

        k_vec_sim = self.sim_auto['k_1d_from_2d']['center']
        k_vec = self.data_xspec['0modes']['k_vec'].value

        assert np.array_equal(k_vec_sim, k_vec), "simulation and data unaligned kvecs"
        log_k_vec = np.log10(k_vec)

        for treatment in self.data_xspec:
            est_data_cov = np.zeros_like(sim_xspec_cov)

            data_xspec_mean = self.data_xspec[treatment]['mean'].value
            data_auto_mean = self.data_auto[treatment]['mean'].value
            data_xspec_std = self.data_xspec[treatment]['std'].value
            data_auto_std = self.data_auto[treatment]['std'].value
            print treatment

            outfile = "%s/pk_terms_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, data_xspec_mean,
                                         data_xspec_std, data_auto_mean,
                                         data_auto_std, sim_xspec_mean,
                                         sim_auto_mean,
                                         outfile=outfile)

            data_cov_fac = 7. / 6. * data_xspec_mean * data_xspec_mean
            data_cov_fac += 2. / 3. * data_xspec_mean * data_auto_mean
            data_cov_fac += 1. / 6. * data_auto_mean * data_auto_mean

            sim_cov_fac = 7. / 6. * sim_xspec_mean * sim_xspec_mean
            sim_cov_fac += 2. / 3. * sim_xspec_mean * sim_auto_mean
            sim_cov_fac += 1. / 6. * sim_auto_mean * sim_auto_mean

            cov_fac_ratio = data_cov_fac / sim_cov_fac
            error_multiplier = np.outer(np.sqrt(cov_fac_ratio),
                                        np.sqrt(cov_fac_ratio))

            est_data_cov = sim_xspec_cov * error_multiplier

            np.set_printoptions(threshold=np.nan)
            print est_data_cov

            outfile = "%s/cov_fac_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, data_cov_fac, sim_cov_fac,
                                         outfile=outfile)

            outplot_cov_file = "%s/covariance_%s.png" % \
                                   (self.params['outdir'], treatment)

            plot_slice.simpleplot_2D(outplot_cov_file, est_data_cov,
                                     log_k_vec, log_k_vec,
                                     ["logk", "logk"], 1.,
                                     "1D covariance", "cov")

            new_data_std = np.sqrt(np.diag(est_data_cov))

            outfile = "%s/pk_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, data_xspec_mean,
                                         data_xspec_std, new_data_std,
                                         outfile=outfile)

