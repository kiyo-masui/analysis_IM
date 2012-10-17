from plotting import plot_slice
import numpy as np
from utils import data_paths as dp
from utils import file_tools
from quadratic_products import pwrspec_estimator as pe
from quadratic_products import power_spectrum as ps
from kiyopy import parse_ini
import shelve
import h5py
import glob


crosspowerparams_init = {
        "p_data": "test_map",
        "p_mock": "test_map",
        "apply_2d_beamtransfer": None,
        "apply_2d_modetransfer": None,
        "use_noiseweights_2dto1d": True,
        "kmin_xy": [None, None],
        "outdir": "./"
               }
crosspowerprefix = 'crosspower_'


class CompileCrosspower(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          crosspowerparams_init,
                                          prefix=crosspowerprefix)
            print self.params

    def execute(self, processes):
        print "using data: ", self.params["p_data"]
        print "using mock: ", self.params["p_mock"]
        pwr_data = ps.PowerSpectrum(self.params["p_data"])
        pwr_mock = ps.PowerSpectrum(self.params["p_mock"])

        modetransfer_2d = None
        beamtransfer_2d = None
        if (self.params["apply_2d_beamtransfer"] is not None) or \
           (self.params["apply_2d_modetransfer"] is not None):

            if self.params["apply_2d_modetransfer"] is not None:
                print "Applying 2d transfer from " + \
                      self.params["apply_2d_modetransfer"]
                #trans_shelve = shelve.open(
                #                    self.params["apply_2d_modetransfer"], "r")
                #modetransfer_2d = trans_shelve["transfer_2d"]
                #trans_shelve.close()
                modetransfer_2d = h5py.File(
                                    self.params["apply_2d_modetransfer"], "r")

            if self.params["apply_2d_beamtransfer"] is not None:
                print "Applying 2d transfer from " + \
                      self.params["apply_2d_beamtransfer"]
                #trans_shelve = shelve.open(
                #                    self.params["apply_2d_beamtransfer"], "r")
                #beamtransfer_2d = trans_shelve["transfer_2d"]
                #trans_shelve.close()
                beamtransfer_2d = h5py.File(
                                    self.params["apply_2d_beamtransfer"], "r")

            transfer_dict = {}
            for treatment in pwr_data.treatment_cases:
                if modetransfer_2d is not None:
                    transfer_dict[treatment] = modetransfer_2d[treatment].value
                    if beamtransfer_2d is not None:
                        transfer_dict[treatment] *= \
                                    beamtransfer_2d["0modes"].value
                else:
                    transfer_dict[treatment] = beamtransfer_2d["0modes"].value

            if (modetransfer_2d is not None) or (beamtransfer_2d is not None):
                pwr_mock.apply_2d_trans_by_treatment(transfer_dict)
                pwr_data.apply_2d_trans_by_treatment(transfer_dict)

        # gather the rms of each 2D k-bin over the realizations of mock
        # catalogs; this is used for the 2D->1D weighting
        #mock2d_agg = pwr_mock.agg_stat_2d_pwrspec(debug="./plotheap/")
        mock2d_agg = pwr_mock.agg_stat_2d_pwrspec()

        weights_2d = {}
        # weight by the variance as determined in the mock runs
        for treatment in pwr_mock.treatment_cases:
            weight_for_treatment = 1. / (mock2d_agg[treatment]["std"] * \
                                         mock2d_agg[treatment]["std"])

            kx = pwr_mock.kx_2d['center']
            ky = pwr_mock.ky_2d['center']
            kmin_x = self.params['kmin_xy'][0]
            kmin_y = self.params['kmin_xy'][1]
            if kmin_x is not None:
                restrict = np.where(kx < kmin_x)
                weight_for_treatment[restrict, :] = 0.

            if kmin_y is not None:
                restrict = np.where(ky < kmin_y)
                weight_for_treatment[:, restrict] = 0.

            weight_for_treatment[np.isnan(weight_for_treatment)] = 0.
            weight_for_treatment[np.isinf(weight_for_treatment)] = 0.

            weights_2d[treatment] = weight_for_treatment

            outplot_power_file = "%s/power_2d_%s.png" % \
                                  (self.params['outdir'], treatment)

            outplot_transfer_file = "%s/transfer_2d_%s.png" % \
                                  (self.params['outdir'], treatment)

            outplot_weight_file = "%s/noiseweight_2d_%s.png" % \
                                  (self.params['outdir'], treatment)

            outplot_count_file = "%s/countweight_2d_%s.png" % \
                                 (self.params['outdir'], treatment)

            # comb is used for autopower AxB etc.; here just use a placeholder
            logkx = np.log10(pwr_mock.kx_2d['center'])
            logky = np.log10(pwr_mock.ky_2d['center'])

            comb = pwr_data.comb_cases[0]
            pwrcase = "%s:%s" % (comb, treatment)
            plot_slice.simpleplot_2D(outplot_power_file,
                                     np.abs(pwr_data.pwrspec_2d[pwrcase]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "2d power",
                                     "log(abs(2d power))", logscale=True)

            if (modetransfer_2d is not None) or (beamtransfer_2d is not None):
                plot_slice.simpleplot_2D(outplot_transfer_file,
                                         transfer_dict[treatment],
                                         logkx, logky,
                                         ["logkx", "logky"], 1., "2d trans",
                                         "2d trans", logscale=False)

            #if treatment == "20modes":
            #    np.set_printoptions(threshold=np.nan)
            #    print mock2d_agg[treatment]["std"]
            #    print np.abs(weights_2d[treatment])

            plot_slice.simpleplot_2D(outplot_weight_file,
                                     np.abs(weights_2d[treatment]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "Log-weight",
                                     "log-weight", logscale=True)

            comb = pwr_mock.comb_cases[0]
            pwrcase = "%s:%s" % (comb, treatment)
            plot_slice.simpleplot_2D(outplot_count_file,
                                     np.abs(pwr_mock.counts_2d[pwrcase]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "Log-counts",
                                     "log-counts", logscale=True)

        if self.params["use_noiseweights_2dto1d"]:
            pwr_data.convert_2d_to_1d(weights_2d=weights_2d)
            pwr_mock.convert_2d_to_1d(weights_2d=weights_2d)
        else:
            # use the counts to weight instead by default
            # just be sure the transfer function is applied
            pwr_data.convert_2d_to_1d()
            pwr_mock.convert_2d_to_1d()

        pwr_data_summary = pwr_data.agg_stat_1d_pwrspec(from_2d=True)
        pwr_mock_summary = pwr_mock.agg_stat_1d_pwrspec(from_2d=True)

        # build the hd5 summary file
        sum_filename = "%s/cross_summary.hd5" % self.params["outdir"]
        summaryfile = h5py.File(sum_filename, "w")
        power_1d_out = summaryfile.create_group("power_1d")
        power_1d_cov_out = summaryfile.create_group("power_1d_cov")
        kvec_out = summaryfile.create_group("kvec")

        k_vec = pwr_data.k_1d_from_2d["center"]
        kvec_out["k_1d_left"] = pwr_data.k_1d_from_2d["left"]
        kvec_out["k_1d_center"] = pwr_data.k_1d_from_2d["center"]
        kvec_out["k_1d_right"] = pwr_data.k_1d_from_2d["right"]

        for treatment in pwr_data.treatment_cases:
            mean_data = pwr_data_summary[treatment]["mean"]
            mean_mock = pwr_mock_summary[treatment]["mean"]
            std_mock = pwr_mock_summary[treatment]["std"]
            cov_mock = pwr_mock_summary[treatment]["cov"]
            corr_mock = pwr_mock_summary[treatment]["corr"]

            power_1d_out[treatment] = mean_data
            power_1d_cov_out[treatment] = cov_mock

            logk = np.log10(k_vec)
            outplot_corr_file = "%s/covariance_1d_%s.png" % \
                                 (self.params['outdir'], treatment)

            plot_slice.simpleplot_2D(outplot_corr_file,
                                     corr_mock,
                                     logk, logk,
                                     ["log k", "log k"], 1., "1d corr",
                                     "1d corr", logscale=False)

            outfile = "%s/pk_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, mean_data, mean_mock, std_mock,
                                         outfile=outfile)

            print "writing to " + outfile

        summaryfile.close()



