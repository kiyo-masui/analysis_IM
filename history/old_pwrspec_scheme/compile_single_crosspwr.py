r"""Find a single crosspower
"""
import numpy as np
from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from utils import file_tools
from optparse import OptionParser
from kiyopy import parse_ini


def batch_single_crosspwr(left_mapkey,
                          right_simkey, right_weightkey,
                          multiplier=1.,
                          inifile=None, datapath_db=None,
                          outdir="./plots",
                          output_tag=None):
    r"""w_left m_left x w_right m_right"""

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(left_mapkey, "type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    crosspwr_collection = {}
    for treatment in map_cases['treatment']:
        dbkeydict = {}
        dbkeydict['map1_key'] = "%s:map;%s" % (left_mapkey, treatment)
        dbkeydict['map2_key'] = right_simkey
        dbkeydict['noiseinv1_key'] = "%s:weight;%s" % (left_mapkey, treatment)
        dbkeydict['noiseinv2_key'] = right_weightkey
        files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                      datapath_db=datapath_db)

        pwrspec_out_signal = caller.execute(files['map1_key'],
                                            files['map2_key'],
                                            files['noiseinv1_key'],
                                            files['noiseinv2_key'],
                                            inifile=inifile)

        if output_tag:
            pwrspec_out_signal[0]['binavg'] *= multiplier
            pwrspec_out_signal[1]['binavg'] *= multiplier

            pwr_1d_from_2d = pe.convert_2d_to_1d(pwrspec_out_signal[0])

            mtag = "%s_%s" % (output_tag, treatment)
            pe.summarize_pwrspec(pwrspec_out_signal[1],
                                 pwr_1d_from_2d,
                                 pwrspec_out_signal[0],
                                 mtag, outdir=outdir)

            crosspwr_collection[treatment] = (pwrspec_out_signal[1],
                                              pwr_1d_from_2d,
                                              pwrspec_out_signal[0])

    if not output_tag:
        caller.multiprocess_stack()
        return None
    else:
        return crosspwr_collection


def wrap_batch_single_crosspwr(inifile, generate=False, outdir="./plots/"):
    r"""Wrapper to the single crosspwr calculator
    """
    params_init = {"left_mapkey": "some preparation of a map, cleaned",
                   "right_simkey": "a simulation to cross it with",
                   "right_weightkey": "weight to use for that sim",
                   "multiplier": "multiply the 1D and 2D spectra",
                   "spec_ini": "ini file for the spectral estimation",
                   "output_tag": "tag identifying the output somehow"}
    prefix="csc_"

    params = parse_ini.parse(inifile, params_init, prefix=prefix)
    print params

    output_tag = "%s_%s" % (params['left_mapkey'], params['output_tag'])
    output_root = "%s/%s/" % (outdir, output_tag)

    if generate:
        output_tag = None

    print output_root, output_tag
    file_tools.mkparents(output_root)
    parse_ini.write_params(params, output_root + 'params.ini',
                           prefix=prefix)

    datapath_db = data_paths.DataPath()

    return batch_single_crosspwr(params["left_mapkey"],
                                 params["right_simkey"],
                                 params["right_weightkey"],
                                 multiplier=params["multiplier"],
                                 inifile=params["spec_ini"],
                                 datapath_db=datapath_db,
                                 outdir=output_root,
                                 output_tag=output_tag)


if __name__ == '__main__':
    parser = OptionParser(usage="usage: %prog [options] filename (-h for help)",
                          version="%prog 1.0")
    parser.add_option("-g", "--generate", action="store_true",
                      dest="generate", default=False,
                      help="regenerate the cache of quadratic products")
    parser.add_option("-o", "--outdir", action="store",
                      dest="outdir", default="./plots/",
                      help="directory to write output data to")
    (optparam, inifile) = parser.parse_args()
    optparam = vars(optparam)

    if len(inifile) != 1:
        parser.error("inifile not specified")

    inifile = inifile[0]
    print optparam

    wrap_batch_single_crosspwr(inifile,
                               generate=optparam['generate'],
                               outdir=optparam['outdir'])
