import numpy as np
from utils import data_paths as dp
from kiyopy import parse_ini
import shelve

params_init = {
        "p_map": "test_map",
        "p_map_plussim": "test_map",
        "p_cleaned_sim": "test_map"
               }
prefix = 'autopower_'

def repackage_1d_power(pwrspec_compilation):
    case_key = "combination:treatment"
    pwr_cases = dp.unpack_cases(pwrspec_compilation.keys(), case_key, divider=":")

    comb_cases = pwr_cases["combination"]
    # note that the dictionary may not be in any order
    comb_cases.sort()

    summary = {}
    for pwrspec_case in pwrspec_compilation:
        pwrspec_entry = pwrspec_compilation[pwrspec_case]
        parameters = pwrspec_entry[0]
        pwrdata_2d = pwrspec_entry[1][0]
        pwrdata_1d = pwrspec_entry[1][1]
        summary[pwrspec_case] = pwrdata_1d["binavg"]
        bin_left = pwrdata_1d["bin_left"]
        bin_center = pwrdata_1d["bin_center"]
        bin_right = pwrdata_1d["bin_right"]
        histo_counts = pwrdata_1d["counts_histo"]
        #print pwrspec_case, parameters

    num_k = bin_center.shape[0]
    num_pk = len(comb_cases)

    summary_treatment = {}
    for treatment in pwr_cases["treatment"]:
        pwr_treatment = np.zeros((num_k, num_pk))
        for comb, comb_index in zip(comb_cases, range(num_pk)):
            pwrcase = "%s:%s" % (comb, treatment)
            pwr_treatment[:, comb_index] = summary[pwrcase]
        summary_treatment[treatment] = pwr_treatment

    return (bin_left, bin_center, bin_right, histo_counts, summary_treatment)


class CompileAutopower(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

    def execute(self, processes):
        pwr_map = shelve.open(self.params["p_map"])
        pwr_map_plussim = shelve.open(self.params["p_map_plussim"])
        pwr_cleaned_sim = shelve.open(self.params["p_cleaned_sim"])

        case_key = "combination:treatment"
        pwr_cases = dp.unpack_cases(pwr_map.keys(), case_key, divider=":")

        pwr_map_1d = repackage_1d_power(pwr_map)
        k_vec = pwr_map_1d[1]
        pwr_map_1d = pwr_map_1d[4]

        pwr_map_plussim_1d = repackage_1d_power(pwr_map_plussim)[4]
        pwr_cleaned_sim_1d = repackage_1d_power(pwr_cleaned_sim)[4]

        reference_pwr = np.mean(pwr_cleaned_sim_1d["0modes"], axis=1)

        for treatment in pwr_cases["treatment"]:
            mean_map = np.mean(pwr_map_1d[treatment], axis=1)
            std_map = np.std(pwr_map_1d[treatment], axis=1)
            mean_map_plussim = np.mean(pwr_map_plussim_1d[treatment], axis=1)
            mean_cleaned_sim = np.mean(pwr_cleaned_sim_1d[treatment], axis=1)
            trans = (mean_map_plussim - mean_map) / reference_pwr

            outfile = open("pk_%s.dat" % treatment, "w")
            for k, p0, tk, pk_mean, pk_err in zip(k_vec, reference_pwr, trans, mean_map, std_map):
                outfile.write("%10.15g %10.15g %10.15g %10.15g %10.15g\n" % (k, p0, tk, pk_mean, pk_err))

if __name__ == '__main__':
    if len(sys.argv) == 2:
        PwrspecCombinations(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'
