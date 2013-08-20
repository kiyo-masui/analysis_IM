r"""Code to estimate power spectra of the GBT data
"""
from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from utils import file_tools
from optparse import OptionParser
from kiyopy import parse_ini
from correlate import compile_crosspwr_transfer as cct


def batch_gbtpwrspec_data_run(map_key, inifile=None, datapath_db=None,
                   output_tag=None, beam_transfer=None,
                   outdir="./plots/",
                   square_1dmodetrans=False,
                   mode_transfer_1d=None,
                   mode_transfer_2d=None):
    r"""Form the pairs of maps1*weight1 x map0*weight0 for calculating the
    auto-power of the GBT data"""

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(map_key, "pair;type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    if output_tag:
        file_tools.mkparents(outdir)

    pwrspec_collection = {}
    for treatment in map_cases['treatment']:
        unique_pairs = data_paths.GBTauto_cross_pairs(map_cases['pair'],
                                                    map_cases['pair'],
                                                    cross_sym="_with_")

        # TODO: make this more elegant
        transfer_2d = None
        if (mode_transfer_2d is not None) and (beam_transfer is None):
            transfer_2d = mode_transfer_2d[treatment][1]

        if (mode_transfer_2d is None) and (beam_transfer is not None):
            transfer_2d = beam_transfer

        if (mode_transfer_2d is not None) and (beam_transfer is not None):
            transfer_2d = mode_transfer_2d[treatment][1] * beam_transfer

        pwr_1d = []
        pwr_2d = []
        pwr_1d_from_2d = []
        for item in unique_pairs:
            dbkeydict = {}
            # NOTE: formerly had "weight" instead of noise_inv
            mapset0 = (map_key, item[0], treatment)
            mapset1 = (map_key, item[1], treatment)
            dbkeydict['map1_key'] = "%s:%s;map;%s" % mapset0
            dbkeydict['map2_key'] = "%s:%s;map;%s" % mapset1
            dbkeydict['noiseinv1_key'] = "%s:%s;noise_inv;%s" % mapset0
            dbkeydict['noiseinv2_key'] = "%s:%s;noise_inv;%s" % mapset1
            files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                      datapath_db=datapath_db)

            pwrspec_out = caller.execute(files['map1_key'],
                                         files['map2_key'],
                                         files['noiseinv1_key'],
                                         files['noiseinv2_key'],
                                         inifile=inifile)

            if output_tag:
                pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec_out[0],
                                      transfer=transfer_2d))

                pwr_2d.append(pwrspec_out[0])
                pwr_1d.append(pwrspec_out[1])

        if output_tag:
            mtag = output_tag + "_%s" % treatment
            if mode_transfer_1d is not None:
                transfunc = mode_transfer_1d[treatment][0]
                if square_1dmodetrans:
                    transfunc *= transfunc
            else:
                transfunc = None

            agg_pwrspec = pe.summarize_agg_pwrspec(pwr_1d,
                                                   pwr_1d_from_2d, pwr_2d, mtag,
                                                   outdir=outdir,
                                                   apply_1d_transfer=transfunc)

            # (mean_1d, std_1d, covmat_1d)
            pwrspec_collection[treatment] = agg_pwrspec

    if output_tag:
        return pwrspec_collection
    else:
        caller.multiprocess_stack()
        return None


def call_data_autopower(basemaps, treatments, inifile=None, generate=False,
                        outdir="./plots/", mode_transfer_1d=None,
                        mode_transfer_2d=None, beam_transfer=None, alttag=None):
    r"""Call a chunk of batch data runs for e.g. different map products
    """
    datapath_db = data_paths.DataPath()

    for base in basemaps:
        for treatment in treatments:
            output_tag = base + treatment
            if alttag:
                output_tag += "_" + output_tag

            output_root = "%s/%s/" % (outdir, output_tag)

            if generate:
                output_tag = None

            batch_gbtpwrspec_data_run(mapname,
                                 inifile=inifile,
                                 datapath_db=datapath_db,
                                 output_tag=output_tag,
                                 outdir=output_root,
                                 beam_transfer=beam_transfer,
                                 mode_transfer_1d=mode_transfer_1d,
                                 mode_transfer_2d=mode_transfer_2d)


def wrap_batch_gbtpwrspec_data_run(inifile, generate=False,
                                    outdir="./plots/"):
    r"""Wrapper to the GBT x GBT calculation"""
    params_init = {"gbt_mapkey": "cleaned GBT map",
                   "mode_transfer_1d_ini": "ini file -> 1d trans. function",
                   "mode_transfer_2d_ini": "ini file -> 2d trans. function",
                   "beam_transfer_ini": "ini file -> 2d beam trans. function",
                   "square_1dmodetrans": False,
                   "spec_ini": "ini file for the spectral estimation",
                   "output_tag": "tag identifying the output somehow"}
    prefix="cp_"

    params = parse_ini.parse(inifile, params_init, prefix=prefix)
    print params

    output_tag = "%s_%s" % (params['gbt_mapkey'], params['output_tag'])
    output_root = "%s/%s/" % (outdir, output_tag)

    if generate:
        output_tag = None

    print output_root
    print output_tag
    file_tools.mkparents(output_root)
    parse_ini.write_params(params, output_root + 'params.ini',
                           prefix=prefix)

    datapath_db = data_paths.DataPath()

    mode_transfer_1d=None
    if params["mode_transfer_1d_ini"]:
        mode_transfer_1d = cct.wrap_batch_crosspwr_transfer(
                                            params["mode_transfer_1d_ini"],
                                            generate=generate,
                                            outdir=outdir)

    return batch_gbtpwrspec_data_run(params["gbt_mapkey"],
                         inifile=params["spec_ini"],
                         datapath_db=datapath_db,
                         outdir=output_root,
                         output_tag=output_tag,
                         beam_transfer=None,
                         square_1dmodetrans = params["square_1dmodetrans"],
                         mode_transfer_1d=mode_transfer_1d,
                         mode_transfer_2d=None)


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

    wrap_batch_gbtpwrspec_data_run(inifile,
                                   generate=optparam['generate'],
                                   outdir=optparam['outdir'])
