r"""Code to estimate power spectra of the GBT data
"""
from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler


def batch_data_run(map_key, inifile=None, datapath_db=None,
                   usecache_output_tag=None, beam_transfer=None,
                   outdir="./plot_data_v2/",
                   mode_transfer_1d=None,
                   mode_transfer_2d=None):
    r"""Form the pairs of maps1*weight1 x map0*weight0 for calculating the
    auto-power of the GBT data"""

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(map_key, "pair;type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if usecache_output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

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

            if usecache_output_tag:
                pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec_out[0],
                                      transfer=transfer_2d))

                pwr_2d.append(pwrspec_out[0])
                pwr_1d.append(pwrspec_out[1])

                mtag = usecache_output_tag + "_%s" % treatment
                if mode_transfer_1d is not None:
                    transfunc = mode_transfer_1d[treatment][0]
                else:
                    transfunc = None

                pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, mtag,
                                     outdir=outdir,
                                     apply_1d_transfer=transfunc)

    if usecache_output_tag:
        # TODO: could output the P(k)'s for the various treatments
        return None
    else:
        caller.multiprocess_stack()
        return None


def call_data_autopower(basemaps, treatments, inifile=None, generate=False,
                        outdir="./plot_data_v2/", mode_transfer_1d=None,
                        mode_transfer_2d=None, beam_transfer=None):
    r"""Call a chunk of batch data runs for e.g. different map products
    """
    datapath_db = data_paths.DataPath()
    # TODO: put transfer functions in the naming tags

    for base in basemaps:
        for treatment in treatments:
            mapname = base + treatment

            usecache_output_tag = None
            if not generate:
                usecache_output_tag = mapname

            batch_data_run(mapname,
                           inifile=inifile,
                           datapath_db=datapath_db,
                           usecache_output_tag=usecache_output_tag,
                           beam_transfer=beam_transfer,
                           outdir=outdir,
                           mode_transfer_1d=mode_transfer_1d,
                           mode_transfer_2d=mode_transfer_2d)


def batch_wigglez_automock_run(mock_key, sel_key,
                               inifile=None, datapath_db=None,
                               usecache_output_tag=None):
    r"""TODO: make this work; wrote this but never really needed it yet
    """
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if usecache_output_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    for index in mock_cases['realization']:
        dbkeydict = {}
        dbkeydict['map1_key'] = "%s:%s" % (mock_key, index)
        dbkeydict['map2_key'] = "%s:%s" % (mock_key, index)
        dbkeydict['noiseinv1_key'] = sel_key
        dbkeydict['noiseinv2_key'] = sel_key
        files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                    datapath_db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    caller.multiprocess_stack()
