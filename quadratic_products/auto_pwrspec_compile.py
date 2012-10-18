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


autonoiseweightparams_init = {
        "p_noise": "test_map",
        "apply_2d_beamtransfer": None,
        "apply_2d_modetransfer": None,
        "outfile": "noise_weight.hd5"
               }
autonoiseweightprefix = 'autonoiseweight_'


class CompileAutoNoiseweight(object):
    r"""combine AxA, BxB, CxC, DxD into a noise bias estimate"""
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
                            self.params["apply_2d_modetransfer"])


        weightfile = h5py.File(self.params["outfile"], "w")
        for treatment in pwr_noise.treatment_cases:
            #weight = noise2d_agg[treatment]["mean"] - \
            #         signal2d_agg[treatment]["mean"]
            # note that we sum in constant |k| annuli so k^3 in k^3 P(k) is a
            # constant factor
            comb = pwr_map.comb_cases[0]
            pwrcase = "%s:%s" % (comb, treatment)
            weight = np.abs(pwr_noise.counts_2d[pwrcase]) / 4.
            weight /= noise2d_agg[treatment]["mean"] ** 2.
            if transfer_dict is not None:
                weight *= transfer_dict[treatment] ** 2.

            weightfile[treatment] = weight

        weightfile.close()


fastautopowerparams_init = {
        "p_map": "test_map",
        "p_map_plussim": "test_map",
        "p_cleaned_sim": "test_map",
        "noiseweights_2dto1d": None,
        "apply_2d_transfer": None,
        "outdir": "./"
               }
fastautopowerprefix = 'fastautopower_'


class CompileFastAutopower(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, fastautopowerparams_init,
                                          prefix=fastautopowerprefix)

    def execute(self, processes):
        pwr_map = ps.PowerSpectrum(self.params["p_map"])
        pwr_map_plussim = ps.PowerSpectrum(self.params["p_map_plussim"])
        pwr_cleaned_sim = ps.PowerSpectrum(self.params["p_cleaned_sim"])

        #plussim_basename = self.params["p_map_plussim"].split(".")[0]
        #plussim_list = glob.glob("%s_[0-9]*.shelve" % plussim_basename)
        #plussim_1d_list = []
        #for filename in plussim_list:
        #    print filename
        #    pwr_map_plussim = shelve.open(filename, "r")
        #    plussim_1d_list.append(pe.repackage_1d_power(pwr_map_plussim)[4])
        #    pwr_map_plussim.close()

        #avg_plussim = {}
        #for treatment in pwr_cases["treatment"]:
        #    shape_plussim = plussim_1d_list[0][treatment].shape
        #    ndim_plussim = len(shape_plussim)
        #    num_plussim = len(plussim_1d_list)
        #    shape_avg = shape_plussim + (num_plussim,)
        #    avg_treatment = np.zeros(shape_avg)
        #    for (plussim_realization, ind) in \
        #        zip(plussim_1d_list, range(num_plussim)):
        #        avg_treatment[..., ind] = plussim_realization[treatment]

        #    avg_plussim[treatment] = np.mean(avg_treatment, axis=ndim_plussim)
        #    print avg_plussim[treatment].shape, shape_plussim

        # TODO!!!!!: is this a proper function of treatment?
        if self.params["apply_2d_transfer"] is not None:
            # copy the same beam transfer function for all cases

            trans_shelve = shelve.open(self.params["apply_2d_transfer"])
            transfer_2d = trans_shelve["transfer_2d"]
            trans_shelve.close()

            transfer_dict = {}
            for treatment in pwr_map.treatment_cases:
                transfer_dict[treatment] = transfer_2d

            pwr_map.apply_2d_trans_by_treatment(transfer_dict)
            pwr_map_plussim.apply_2d_trans_by_treatment(transfer_dict)
            pwr_cleaned_sim.apply_2d_trans_by_treatment(transfer_dict)

        # also combine the AxB, etc. into a signal piece that is subtracted
        signal2d_agg = pwr_map.agg_stat_2d_pwrspec()

        for treatment in pwr_map.treatment_cases:
            outplot_power_file = "%s/power_2d_%s.png" % \
                                  (self.params['outdir'], treatment)

            outplot_transfer_file = "%s/transfer_2d_%s.png" % \
                                  (self.params['outdir'], treatment)

            outplot_count_file = "%s/countweight_2d_%s.png" % \
                                 (self.params['outdir'], treatment)

            # comb is used for autopower AxB etc.; here just use a placeholder
            logkx = np.log10(pwr_map.kx_2d['center'])
            logky = np.log10(pwr_map.ky_2d['center'])

            comb = pwr_map.comb_cases[0]
            pwrcase = "%s:%s" % (comb, treatment)

            plot_slice.simpleplot_2D(outplot_power_file,
                                     np.abs(pwr_map.pwrspec_2d[pwrcase]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "2d power",
                                     "log(abs(2d power))", logscale=True)

            if self.params["apply_2d_transfer"] is not None:
                plot_slice.simpleplot_2D(outplot_transfer_file,
                                         transfer_dict[treatment],
                                         logkx, logky,
                                         ["logkx", "logky"], 1., "2d trans",
                                         "2d trans", logscale=False)

            plot_slice.simpleplot_2D(outplot_count_file,
                                     np.abs(pwr_map.counts_2d[pwrcase]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "Log-counts",
                                     "log-counts", logscale=True)

        if self.params["noiseweights_2dto1d"] is not None:
            weightfile = h5py.File(self.params["noiseweights_2dto1d"], "r")
            weights_2d = {}
            for treatment in pwr_map.treatment_cases:
                weights_2d[treatment] = weightfile[treatment].value

                outplot_weight_file = "%s/noiseweight_2d_%s.png" % \
                                      (self.params['outdir'], treatment)


                plot_slice.simpleplot_2D(outplot_weight_file,
                                     np.abs(weights_2d[treatment]),
                                     logkx, logky,
                                     ["logkx", "logky"], 1., "Log-weight",
                                     "log-weight", logscale=True)

            weightfile.close()

            pwr_map.convert_2d_to_1d(weights_2d=weights_2d)
            pwr_map_plussim.convert_2d_to_1d(weights_2d=weights_2d)
            pwr_cleaned_sim.convert_2d_to_1d(weights_2d=weights_2d)
        else:
            pwr_map.convert_2d_to_1d()
            pwr_map_plussim.convert_2d_to_1d()
            pwr_cleaned_sim.convert_2d_to_1d()

        pwr_map_summary = pwr_map.agg_stat_1d_pwrspec(from_2d=True)

        pwr_map_plussim_summary = \
            pwr_map_plussim.agg_stat_1d_pwrspec(from_2d=True)

        pwr_cleaned_sim_summary = \
            pwr_cleaned_sim.agg_stat_1d_pwrspec(from_2d=True)

        reference_pwr = pwr_cleaned_sim_summary["0modes"]["mean"]
        k_vec = pwr_map.k_1d_from_2d["center"]

        for treatment in pwr_map.treatment_cases:
            mean_map_plussim = pwr_map_plussim_summary[treatment]["mean"]
            mean_map = pwr_map_summary[treatment]["mean"]
            std_map = pwr_map_summary[treatment]["std"]
            trans = (mean_map_plussim - mean_map) / reference_pwr

            #mean_map_plussim = np.mean(avg_plussim[treatment], axis=1)
            ##trans = (pwr_map_plussim_1d[treatment]-pwr_map_1d[treatment]) / \
            ##        pwr_cleaned_sim_1d["0modes"]
            #trans = (avg_plussim[treatment]-pwr_map_1d[treatment]) / \
            #        pwr_cleaned_sim_1d["0modes"]
            #corrected_pwr = pwr_map_1d[treatment]/trans
            #trans = np.mean(trans, axis=1)
            #mean_map = np.mean(corrected_pwr, axis=1)
            #std_map = np.mean(corrected_pwr, axis=1)

            outfile = "%s/pk_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, reference_pwr, trans,
                                         mean_map, std_map, outfile=outfile)
            print "writing to " + outfile


autopowerparams_init = {
        "p_map": "test_map",
        "p_noise": "test_map",
        "noiseweights_2dto1d": None,
        "apply_2d_transfer": None,
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

        # TODO!!!!!: is this a proper function of treatment?
        if self.params["apply_2d_transfer"] is not None:
            # copy the same beam transfer function for all cases
            trans_shelve = shelve.open(self.params["apply_2d_transfer"])
            transfer_2d = trans_shelve["transfer_2d"]
            trans_shelve.close()

            transfer_dict = {}
            for treatment in pwr_map.treatment_cases:
                transfer_dict[treatment] = transfer_2d

            pwr_map.apply_2d_trans_by_treatment(transfer_dict)

        # also combine the AxB, etc. into a signal piece that is subtracted
        signal2d_agg = pwr_map.agg_stat_2d_pwrspec()

        if self.params["noiseweights_2dto1d"] is not None:
            weightfile = h5py.File(self.params["noiseweights_2dto1d"], "r")
            weights_2d = {}
            for treatment in pwr_map.treatment_cases:
                weights_2d[treatment] = weightfile[treatment].value

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

        pwr_map_summary = pwr_map.agg_stat_1d_pwrspec(from_2d=True)
        k_vec = pwr_map.k_1d_from_2d["center"]
        for treatment in pwr_map.treatment_cases:
            mean_map = pwr_map_summary[treatment]["mean"]
            std_map = pwr_map_summary[treatment]["std"]

            outfile = "%s/pk_%s.dat" % (self.params["outdir"], treatment)
            file_tools.print_multicolumn(k_vec, mean_map, std_map, outfile=outfile)
            print "writing to " + outfile


