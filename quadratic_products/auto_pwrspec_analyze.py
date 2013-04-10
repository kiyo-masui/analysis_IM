from plotting import plot_slice
import numpy as np
import scipy as sp
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


def decorrelate(data_mean, sim_mean, data_cov, use_sqrtm=True, row_dom=True):
    r"""Decorrelate the power spectra
    use_sqrtm uses numpy's sqrtm rather than sqrt through eigh
    """
    whfinite = np.where(np.isfinite(np.diag(data_cov)))[0]
    whslice = slice(np.min(whfinite), np.max(whfinite) + 1)
    rcov = data_cov[whslice, whslice]
    rdata = data_mean[whslice]
    rsim = sim_mean[whslice]

    if use_sqrtm:
        decorr = np.real(sp.linalg.sqrtm(np.linalg.inv(rcov)))
    else:
        eigenval, eigenvec = np.linalg.eigh(np.linalg.inv(rcov))
        decorr = np.dot(eigenvec,
                        np.dot(np.sqrt(np.diag(eigenval)), eigenvec.T))
        #print "CHECK: ", np.linalg.inv(rcov) - np.dot(decorr, decorr)

    decorr_mean = data_mean.copy()
    decorr_cov = data_cov.copy()

    if row_dom:
        decorr_norm = np.sum(decorr, axis=0)
        decorr /= decorr_norm[:,None]

        decorr = np.dot(np.diag(rsim), np.dot(decorr, np.diag(1./rsim)))

        decorr_mean[whslice] = np.dot(decorr, rdata)

        decorr_cov[whslice, whslice] = np.dot(decorr, np.dot(rcov,
                                              decorr.transpose()))

        #print np.dot(decorr, np.dot(rcov, decorr.transpose()))

    else:
        decorr_norm = np.sum(decorr, axis=1)
        decorr /= decorr_norm[None,:]

        decorr = np.dot(np.diag(1./rsim), np.dot(decorr, np.diag(rsim)))

        decorr_mean[whslice] = np.dot(decorr.transpose(), rdata)

        decorr_cov[whslice, whslice] = np.dot(decorr.transpose(),
                                              np.dot(rcov, decorr))

        #print np.dot(decorr.transpose(), np.dot(rcov, decorr))

    return (decorr_mean, decorr_cov, decorr)


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
        self.sim_auto = h5py.File(self.params["sim_auto_summary"], "r")
        self.sim_xspec = h5py.File(self.params["sim_xspec_summary"], "r")

    def execute(self, processes):
        #sim_auto_stats = self.sim_auto['stats']['0modes']['pk_1d_from_2d_stat']
        #sim_xspec_stats = self.sim_xspec['stats']['0modes']['pk_1d_from_2d_stat']
        sim_auto_stats = self.sim_auto['0modes']['pk_1d_from_2d_stat']
        sim_xspec_stats = self.sim_xspec['0modes']['pk_1d_from_2d_stat']

        sim_xspec_mean = sim_xspec_stats['mean'].value
        sim_xspec_cov = sim_xspec_stats['cov'].value
        sim_xspec_corr = sim_xspec_stats['corr'].value
        sim_auto_mean = sim_auto_stats['mean'].value

        #k_vec_sim = self.sim_auto['k_1d_from_2d']['center']
        k_vec_sim = self.sim_auto['0modes']['k_1d_from_2d'].value
        k_vec = self.data_xspec['0modes']['k_vec'].value

        assert np.array_equal(k_vec_sim, k_vec), "simulation and data unaligned kvecs"
        log_k_vec = np.log10(k_vec)

        for treatment in self.data_xspec:
            # the Gaussian error covariance model
            est_data_cov = np.zeros_like(sim_xspec_cov)

            # a covariance model calibrated on the stdev of 6 pairs
            pairs_data_cov = np.zeros_like(sim_xspec_cov)

            data_xspec_mean = self.data_xspec[treatment]['mean'].value
            data_auto_mean = self.data_auto[treatment]['mean'].value
            data_xspec_std = self.data_xspec[treatment]['std'].value / np.sqrt(6.)
            data_auto_std = self.data_auto[treatment]['std'].value / np.sqrt(6.)
            print treatment

            outfile = "%s/pk_terms_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, data_xspec_mean,
                                         data_xspec_std, data_auto_mean,
                                         data_auto_std, sim_xspec_mean,
                                         sim_auto_mean,
                                         outfile=outfile)

            # calibrate the covariance with Gaussian errors
            sample_variance = 7. / 6. * data_xspec_mean * data_xspec_mean
            noise_sample_variance = 2. / 3. * data_xspec_mean * data_auto_mean
            noise_variance = 1. / 6. * data_auto_mean * data_auto_mean
            data_cov_fac = sample_variance + noise_sample_variance + \
                           noise_variance

            sim_cov_fac = 7. / 6. * sim_xspec_mean * sim_xspec_mean
            sim_cov_fac += 2. / 3. * sim_xspec_mean * sim_auto_mean
            sim_cov_fac += 1. / 6. * sim_auto_mean * sim_auto_mean

            # print contibutions
            outfile = "%s/pkerr_contrib_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, np.sqrt(sample_variance),
                                         np.sqrt(noise_sample_variance),
                                         np.sqrt(noise_variance),
                                         outfile=outfile)


            outfile = "%s/cov_fac_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, data_cov_fac, sim_cov_fac,
                                         outfile=outfile)

            # find the Gaussian error model covariance
            cov_fac_ratio = data_cov_fac / sim_cov_fac
            error_multiplier = np.outer(np.sqrt(cov_fac_ratio),
                                        np.sqrt(cov_fac_ratio))

            est_data_cov = sim_xspec_cov * error_multiplier

            np.set_printoptions(threshold=np.nan)
            #print est_data_cov

            # calibrate the covariance with stdev (empirical errors)
            pairscov_fac_ratio = data_xspec_std ** 2.
            pairs_error_multiplier = np.outer(np.sqrt(pairscov_fac_ratio),
                                              np.sqrt(pairscov_fac_ratio))

            pairs_data_cov = sim_xspec_corr * pairs_error_multiplier

            np.set_printoptions(threshold=np.nan)
            #print pairs_data_cov

            # plot the Gaussian error covariance
            outplot_cov_file = "%s/covariance_%s.png" % \
                                   (self.params['outdir'], treatment)

            plot_slice.simpleplot_2D(outplot_cov_file, est_data_cov,
                                     log_k_vec, log_k_vec,
                                     ["logk", "logk"], 1.,
                                     "Gaussian 1D covariance", "cov")

            # plot the pair error covariance
            outplot_cov_file = "%s/pair_covariance_%s.png" % \
                                   (self.params['outdir'], treatment)

            plot_slice.simpleplot_2D(outplot_cov_file, pairs_data_cov,
                                     log_k_vec, log_k_vec,
                                     ["logk", "logk"], 1.,
                                     "1D covariance from pairs", "cov")


            # record the 1D errors
            est_data_std = np.sqrt(np.diag(est_data_cov))
            pairs_data_std = np.sqrt(np.diag(pairs_data_cov))

            outfile = "%s/pk_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, data_xspec_mean,
                                         data_xspec_std, est_data_std,
                                         pairs_data_std,
                                         outfile=outfile)

            # decorrelate
            #(demean, decov, dewindow) = decorrelate(data_xspec_mean, pairs_data_cov)
            (demean, decov, dewindow) = decorrelate(data_xspec_mean,
                                                    sim_xspec_mean, est_data_cov)

            (desim, desimcov, desimwindow) = decorrelate(sim_xspec_mean,
                                                    sim_xspec_mean, est_data_cov)

            outfile = "%s/pk_decorrtest_%s.dat" % (self.params["outdir"], treatment)

            file_tools.print_multicolumn(k_vec, desim, sim_xspec_mean,
                                         outfile=outfile)


            decorr_errors = np.sqrt(np.diag(decov))

            outplot_cov_file = "%s/windows_%s.png" % \
                                   (self.params['outdir'], treatment)

            wvec = range(dewindow.shape[0])
            plot_slice.simpleplot_2D(outplot_cov_file, dewindow,
                                     wvec, wvec,
                                     ["logk", "logk"], 1.,
                                     "1D covariance from pairs", "cov")

            # plot the decorr error covariance
            outplot_cov_file = "%s/decorr_covariance_%s.png" % \
                                   (self.params['outdir'], treatment)

            plot_slice.simpleplot_2D(outplot_cov_file, decov,
                                     log_k_vec, log_k_vec,
                                     ["logk", "logk"], 1.,
                                     "1D covariance from pairs", "cov")

            outfile = "%s/pk_decorr_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, demean, decorr_errors,
                                         outfile=outfile)
