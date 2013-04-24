r"""Code to calculate:
P = ( m_a P_a + m_b P_b ) / (m_c P_c + m_d P_d )

The native output of the power spectrum is (mean_1d, std_1d, covmat_1d)
"""
import copy
import numpy as np
from correlate import compile_pwrspec as cp
from correlate import pwrspec_estimation as pe
from utils import file_tools
from optparse import OptionParser
from kiyopy import parse_ini


def add_pwrspec(pwr_a, pwr_b):
    """assuming errors are uncorrelated"""
    pwr_out = copy.deepcopy(pwr_a)
    pwr_out["mean_1d"] = pwr_a["mean_1d"] + pwr_b["mean_1d"]
    pwr_out["std_1d"] = np.sqrt(pwr_a["std_1d"] * pwr_a["std_1d"] +
                                pwr_b["std_1d"] * pwr_b["std_1d"])
    pwr_out["covmat_1d"] = pwr_a["covmat_1d"] + pwr_b["covmat_1d"]

    return pwr_out


def mul_pwrspec(pwr, mul):
    pwr_out = copy.deepcopy(pwr)
    pwr_out["mean_1d"] = pwr["mean_1d"] * mul
    pwr_out["std_1d"] = pwr["std_1d"] * mul
    pwr_out["covmat_1d"] = pwr["covmat_1d"] * mul * mul

    return pwr_out


# TODO: implement some form of error propagation here... not the best
def divide_pwrspec(pwr_a, pwr_b):
    pwr_out = copy.deepcopy(pwr_a)
    pwr_out["mean_1d"] = pwr_a["mean_1d"]/pwr_b["mean_1d"]
    pwr_out["std_1d"] = np.zeros_like(pwr_a["std_1d"])
    pwr_out["covmat_1d"] = np.zeros_like(pwr_a["covmat_1d"])

    return pwr_out


# TODO implement cross and type=cross option to indicate it
def spectrum_arithmetic(inifile, outdir="./plots/"):
    r"""Perform the operations"""
    params_init = {"pwr_a_ini": None,
                   "pwr_b_ini": None,
                   "pwr_c_ini": None,
                   "pwr_d_ini": None,
                   "mul_a": "1.",
                   "mul_b": "1.",
                   "mul_c": "1.",
                   "mul_d": "1.",
                   "output_tag": "tag identifying the output somehow"}
    prefix="sa_"

    params = parse_ini.parse(inifile, params_init, prefix=prefix)
    print params

    output_root = "%s/%s/" % (outdir, params["output_tag"])
    print output_root

    file_tools.mkparents(output_root)
    parse_ini.write_params(params, output_root + 'params.ini',
                           prefix=prefix)

    if not params["pwr_a_ini"]:
        print "ERROR: at least spectrum A must be given (and C for ratios)"
        return

    pwr_a_all = cp.wrap_batch_gbtpwrspec_data_run(params["pwr_a_ini"],
                                                  generate=False,
                                                  outdir=outdir)

    if params['pwr_b_ini']:
        pwr_b_all = cp.wrap_batch_gbtpwrspec_data_run(params["pwr_b_ini"],
                                                      generate=False,
                                                      outdir=outdir)

    if params['pwr_c_ini']:
        pwr_c_all = cp.wrap_batch_gbtpwrspec_data_run(params["pwr_c_ini"],
                                                      generate=False,
                                                      outdir=outdir)

    if params['pwr_d_ini']:
        pwr_d_all = cp.wrap_batch_gbtpwrspec_data_run(params["pwr_d_ini"],
                                                      generate=False,
                                                      outdir=outdir)

    for treatment in pwr_a_all:
        pwr_a = pwr_a_all[treatment]
        pwr_a = mul_pwrspec(pwr_a, params['mul_a'])
        pwr_numerator = copy.deepcopy(pwr_a)

        if params['pwr_b_ini']:
            pwr_b = pwr_b_all[treatment]
            pwr_b = mul_pwrspec(pwr_b, params['mul_b'])
            pwr_numerator = add_pwrspec(pwr_numerator, pwr_b)

        if params['pwr_c_ini']:
            pwr_c = pwr_c_all[treatment]
            pwr_c = mul_pwrspec(pwr_c, params['mul_c'])
            pwr_denominator = copy.deepcopy(pwr_c)

            if params['pwr_d_ini']:
                pwr_d = pwr_d_all[treatment]
                pwr_d = mul_pwrspec(pwr_d, params['mul_d'])
                pwr_denominator = add_pwrspec(pwr_denominator, pwr_d)

            pwr_final = divide_pwrspec(pwr_numerator, pwr_denominator)
        else:
            pwr_final = pwr_numerator

        filename = "%s/%s_%s.dat" % (output_root, params['output_tag'], treatment)
        outfile = open(filename, "w")
        for specdata in zip(pwr_final['bin_left'], pwr_final['bin_center'],
                            pwr_final['bin_right'], pwr_final['counts_histo'],
                            pwr_final['mean_1d'], pwr_final['std_1d']):
            outfile.write(("%10.15g " * 6 + "\n") % specdata)

        outfile.close()


if __name__ == '__main__':
    parser = OptionParser(usage="usage: %prog [options] filename (-h for help)",
                          version="%prog 1.0")
    parser.add_option("-o", "--outdir", action="store",
                      dest="outdir", default="./plots/",
                      help="directory to write output data to")
    (optparam, inifile) = parser.parse_args()
    optparam = vars(optparam)

    if len(inifile) != 1:
        parser.error("inifile not specified")

    inifile = inifile[0]
    print optparam

    spectrum_arithmetic(inifile, outdir=optparam['outdir'])
