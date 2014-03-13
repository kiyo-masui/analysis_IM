import numpy as np
from kiyopy import parse_ini
import h5py
from utils import misc as utils

crosssumparams_init = {
        "summary_file": "file",
        "theory_curve": "file",
        "krange": [0.06, 1.],
        "outdir": "./",
        "pk_column": 4
               }
crosssumprefix = 'crosssum_'


class SummarizeCross(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict

        if parameter_file:
            self.params = parse_ini.parse(parameter_file, crosssumparams_init,
                                          prefix=crosssumprefix)

    def execute(self, processes):
        summary_outname = self.params["outdir"] + "theory_fit.dat"
        print "opening: ", self.params["summary_file"]
        print "to summary: ", summary_outname
        cross_summary = h5py.File(self.params["summary_file"], "r")
        # this has: kvec, power_1d, power_1d_cov

        print "using col:%d of %s" % (self.params["pk_column"],
                                      self.params["theory_curve"])

        theory_curve = np.genfromtxt(self.params["theory_curve"])
        theory_curve = theory_curve[:, self.params["pk_column"]]
        #print theory_curve

        bin_center = cross_summary["kvec"]["k_1d_center"].value
        (kmin, kmax) = tuple(self.params["krange"])
        restrict = np.where(np.logical_and(bin_center > kmin,
                                           bin_center < kmax))
        res_slice = slice(min(restrict[0]), max(restrict[0]) + 1)
        #restrict_alt = np.where(restrict)[0][np.newaxis, :]
        #restricted_cov = covmock[restrict_alt][0]

        summary_outfile = open(summary_outname, "w")
        for treatment in cross_summary["power_1d"]:

            pwr = cross_summary["power_1d"][treatment].value
            cov = cross_summary["power_1d_cov"][treatment].value

            amplitude = utils.ampfit(pwr[res_slice],
                                     cov[res_slice, res_slice],
                                     theory_curve[res_slice])

            nmodes = float(treatment.split("modes")[0])
            snr = amplitude['amp'] / amplitude['error']
            outstr = "%d " % nmodes
            outstr += "%(amp)10.5f %(error)10.5f " % amplitude
            outstr += "%(chi2)10.5f %(dof)10.5f %(pte)10.5f " % amplitude
            outstr += "%10.5f" % snr

            print outstr
            summary_outfile.write(outstr + "\n")
