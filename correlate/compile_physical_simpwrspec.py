r"""Use the batch job handler to produce the simulation power spectra and then
produce outputs.
"""
from utils import data_paths
from utils import batch_handler
from correlate import pwrspec_estimation as pe
from utils import file_tools


def batch_physical_sim_run(sim_key, inifile=None, datapath_db=None,
                           output_tag=None,
                           outdir="./plots/"):
    """Test the power spectral estimator using simulations"""
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(sim_key, "realization")

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    generate = False if output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    if output_tag:
        output_root = "%s/%s/" % (outdir, output_tag)
        file_tools.mkparents(output_root)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in mock_cases['realization']:
        mapfile = datapath_db.fetch("%s:%s" % (sim_key, index))

        pwrspec_out = caller.execute(mapfile, mapfile, inifile=inifile)

        if output_tag:
            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec_out[0],
                                                      transfer=None))

            pwr_2d.append(pwrspec_out[0])
            pwr_1d.append(pwrspec_out[1])

    if output_tag:
        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d,
                             output_tag, outdir=output_root)

        retval = (pwr_1d, pwr_1d_from_2d, pwr_2d)
    else:
        caller.multiprocess_stack()
        retval = None

    return retval


def wrap_batch_physical_sim_run(inifile, generate=False, outdir="./plots/"):
    r"""Wrapper to the physical sim processing
    """
    params_init = {"simkey": "point to the sims to run",
                   "spec_ini": "ini file for the spectral estimation",
                   "output_tag": "tag identifying the output somehow"}
    prefix="csp_"

    params = parse_ini.parse(inifile, params_init, prefix=prefix)
    print params

    output_tag = "%s_%s" % (params['simkey'], params['output_tag'])
    output_root = "%s/%s/" % (outdir, output_tag)

    if generate:
        output_tag = None

    print output_root
    print output_tag
    file_tools.mkparents(output_root)
    parse_ini.write_params(params, output_root + 'params.ini',
                           prefix=prefix)

    datapath_db = data_paths.DataPath()

    return batch_physical_sim_run(params["simkey"],
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

    wrap_batch_physical_sim_run(inifile,
                                generate=optparam['generate'],
                                outdir=optparam['outdir'])

