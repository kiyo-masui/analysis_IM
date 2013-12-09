from plotting import plot_slice
import numpy as np
from utils import file_tools
from quadratic_products import pwrspec_estimator as pe
from quadratic_products import power_spectrum as ps
from kiyopy import parse_ini
import shelve
import h5py


autonoiseweightparams_init = {
        "p_noise": "test_map",
        "apply_2d_beamtransfer": None,
        "apply_2d_modetransfer": None,
        "outfile": "noise_weight.hd5"
               }
autonoiseweightprefix = 'autonoiseweight_'


class CompileAutoNoiseweight(object):
    r"""combine AxA, BxB, CxC, DxD into an estimate of the power-spectral noise
    -> inverse variance weight
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          autonoiseweightparams_init,
                                          prefix=autonoiseweightprefix)

    def execute(self, processes):
        pwr_noise = ps.PowerSpectrum(self.params["p_noise"])
        noise2d_agg = pwr_noise.agg_stat_2d_pwrspec()

        transfer_dict = pe.load_transferfunc(
                            self.params["apply_2d_beamtransfer"],
                            self.params["apply_2d_modetransfer"],
                            pwr_noise.treatment_cases)

        weightfile = h5py.File(self.params["outfile"], "w")
        for treatment in pwr_noise.treatment_cases:
            # note that we sum in constant |k| annuli so k^3 in k^3 P(k) is a
            # constant factor
            comb = pwr_noise.comb_cases[0]
            pwrcase = "%s:%s" % (comb, treatment)

            weight = np.abs(pwr_noise.counts_2d[pwrcase])
            weight /= noise2d_agg[treatment]["mean"] ** 2.

            if transfer_dict is not None:
                weight *= transfer_dict[treatment] ** 2.

            # this is a little excessive, but to be safe:
            nanmask = np.isnan(weight)
            infmask = np.isnan(weight)
            zeromask = (pwr_noise.counts_2d[pwrcase] == 0)

            mask = np.where(np.logical_or(
                            np.logical_or(nanmask, infmask), zeromask))

            weight[mask] = 0.

            weightfile[treatment] = weight

        weightfile.close()


autopowerparams_init = {
        "p_map": "test_map",
        "noiseweights_2dto1d": None,
        "apply_2d_beamtransfer": None,
        "apply_2d_modetransfer": None,
        "summaryfile": "./autopower.hd5",
        "kpar_range": None,
        "kperp_range": None,
        "outdir": "./"
               }
autopowerprefix = 'autopower_'


class CompileAutopower(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, autopowerparams_init,
                                          prefix=autopowerprefix)

    def execute(self, processes):
        pwr_map = ps.PowerSpectrum(self.params["p_map"])

        transfer_dict = pe.load_transferfunc(
                            self.params["apply_2d_beamtransfer"],
                            self.params["apply_2d_modetransfer"],
                            pwr_map.treatment_cases)

        pwr_map.apply_2d_trans_by_treatment(transfer_dict)

        # also combine the AxB, etc. into a signal piece that is subtracted
        signal2d_agg = pwr_map.agg_stat_2d_pwrspec()

        kx = pwr_map.kx_2d["center"]
        ky = pwr_map.ky_2d["center"]
        logkx = np.log10(kx)
        logky = np.log10(ky)

        if self.params["kpar_range"] is not None:
            print "restricting k_par to ", self.params["kpar_range"]
            restrict_par = np.where(np.logical_or(
                               ky < self.params["kpar_range"][0],
                               ky > self.params["kpar_range"][1]))
        else:
            restrict_par = None

        if self.params["kperp_range"] is not None:
            print "restricting k_perp to ", self.params["kperp_range"]
            restrict_perp = np.where(np.logical_or(
                               kx < self.params["kperp_range"][0],
                               kx > self.params["kperp_range"][1]))
        else:
            restrict_perp = None

        for treatment in pwr_map.treatment_cases:
            # comb is used for autopower AxB etc.; here just use a placeholder
            comb = pwr_map.comb_cases[0]
            pwrcase = "%s:%s" % (comb, treatment)

            outplot_power_file = "%s/power_2d_%s.png" % \
                                  (self.params['outdir'], treatment)

            plot_slice.simpleplot_2D(outplot_power_file,
                                     np.abs(pwr_map.pwrspec_2d[pwrcase]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "2d power",
                                     "log(abs(2d power))", logscale=True)

            outplot_transfer_file = "%s/transfer_2d_%s.png" % \
                                      (self.params['outdir'], treatment)

            if transfer_dict is not None:
                plot_slice.simpleplot_2D(outplot_transfer_file,
                                         transfer_dict[treatment],
                                         logkx, logky,
                                         ["logkx", "logky"], 1., "2d trans",
                                         "2d trans", logscale=False)
            else:
                plot_slice.simpleplot_2D(outplot_transfer_file,
                                         np.ones_like(pwr_map.counts_2d[pwrcase]),
                                         logkx, logky,
                                         ["logkx", "logky"], 1., "2d trans",
                                         "2d trans", logscale=False)

            outplot_count_file = "%s/countweight_2d_%s.png" % \
                                 (self.params['outdir'], treatment)

            plot_slice.simpleplot_2D(outplot_count_file,
                                     np.abs(pwr_map.counts_2d[pwrcase]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "Log-counts",
                                     "log-counts", logscale=True)


        if self.params["noiseweights_2dto1d"] is not None:
            print "applying 2D noise weights: " + \
                self.params["noiseweights_2dto1d"]

            weightfile = h5py.File(self.params["noiseweights_2dto1d"], "r")
            weights_2d = {}
            for treatment in pwr_map.treatment_cases:
                weights_2d[treatment] = weightfile[treatment].value
                if restrict_perp is not None:
                    weights_2d[treatment][restrict_perp, :] = 0.

                if restrict_par is not None:
                    weights_2d[treatment][:, restrict_par] = 0.

                outplot_weight_file = "%s/noiseweight_2d_%s.png" % \
                                      (self.params['outdir'], treatment)

                plot_slice.simpleplot_2D(outplot_weight_file,
                                     np.abs(weights_2d[treatment]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "Log-weight",
                                     "log-weight", logscale=True)

            weightfile.close()

            pwr_map.convert_2d_to_1d(weights_2d=weights_2d)
        else:
            pwr_map.convert_2d_to_1d()

        # open the output hd5 and attach directories for each treatment
        summary_struct = h5py.File(self.params["summaryfile"], "w")
        summary_by_treatment = {}
        for treatment in pwr_map.treatment_cases:
            summary_by_treatment[treatment] = \
                    summary_struct.create_group(treatment)

        pwr_map_summary = pwr_map.agg_stat_1d_pwrspec(from_2d=True)
        field_list = ["mean", "std", "gauss_std", "cov", "corr"]

        for treatment in pwr_map.treatment_cases:
            # copy data from the summary to the output hd5
            for field in field_list:
                summary_by_treatment[treatment][field] = \
                            pwr_map_summary[treatment][field]

            summary_by_treatment[treatment]["k_vec_left"] = \
                        pwr_map.k_1d_from_2d["left"]

            summary_by_treatment[treatment]["k_vec"] = \
                        pwr_map.k_1d_from_2d["center"]

            summary_by_treatment[treatment]["k_vec_right"] = \
                        pwr_map.k_1d_from_2d["right"]


            outfile = "%s/pk_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(pwr_map.k_1d_from_2d["center"],
                                         pwr_map_summary[treatment]["mean"],
                                         pwr_map_summary[treatment]["std"],
                                         outfile=outfile)

            print "writing to " + outfile

        summary_struct.close()
