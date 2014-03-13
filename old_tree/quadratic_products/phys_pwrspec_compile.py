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

physsimparams_init = {
        "pwr_file": "ok.shelve",
        "outdir": "./"
               }
physsimprefix = 'physsim_'


class CompilePhysSim(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, physsimparams_init,
                                          prefix=physsimprefix)

    def execute(self, processes):
        pwr_physsim = ps.PowerSpectrum(self.params["pwr_file"])

        pwr_physsim.convert_2d_to_1d()
        pwr_physsim_summary = pwr_physsim.agg_stat_1d_pwrspec(from_2d=True)
        k_vec = pwr_physsim.k_1d_from_2d["center"]

        treatment = "phys"
        mean_sim = pwr_physsim_summary[treatment]["mean"]
        std_sim = pwr_physsim_summary[treatment]["std"]

        outfile = "%s/pk_%s.dat" % (self.params["outdir"], treatment)
        file_tools.print_multicolumn(k_vec, mean_sim, std_sim, outfile=outfile)
        print "writing to " + outfile

