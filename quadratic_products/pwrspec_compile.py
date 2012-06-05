import numpy as np
from utils import data_paths as dp
from quadratic_products import pwrspec_estimator as pe
from quadratic_products import power_spectrum as ps
from kiyopy import parse_ini
import shelve
import glob

params_init = {
        "p_map": "test_map",
        "p_map_plussim": "test_map",
        "p_cleaned_sim": "test_map",
        "outdir": "./"
               }
prefix = 'autopower_'

class CompileAutopower(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

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

        pwr_map_summary = pwr_map.agg_stat_1d_pwrspec()
        pwr_map_plussim_summary = pwr_map_plussim.agg_stat_1d_pwrspec()
        pwr_cleaned_sim_summary = pwr_cleaned_sim.agg_stat_1d_pwrspec()
        reference_pwr = pwr_cleaned_sim_summary["0modes"]["mean"]
        k_vec = pwr_map.bin_center_1d

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

            outfile = open("%s/pk_%s.dat" % (self.params["outdir"], treatment), "w")
            for k, p0, tk, pk_mean, pk_err in zip(k_vec, reference_pwr, trans, mean_map, std_map):
                outfile.write("%10.15g %10.15g %10.15g %10.15g %10.15g\n" % (k, p0, tk, pk_mean, pk_err))

if __name__ == '__main__':
    if len(sys.argv) == 2:
        PwrspecCombinations(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'
