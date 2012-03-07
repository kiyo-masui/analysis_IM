r"""Use the batch job handler to produce the simulation power spectra and then
produce outputs.
"""
from utils import data_paths
from utils import batch_handler
from correlate import pwrspec_estimation as pe


def batch_physical_sim_run(sim_key, inifile=None, datapath_db=None,
                           usecache_output_tag=None, transfer=None,
                           outdir="./plot_data_v2/"):
    """Test the power spectral estimator using simulations"""
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(sim_key, "realization")

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    generate = False if usecache_output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in mock_cases['realization']:
        mapfile = datapath_db.fetch("%s:%s" % (sim_key, index))

        pwrspec_out = caller.execute(mapfile, mapfile, inifile=inifile)

        if usecache_output_tag:
            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec_out[0],
                                                  transfer=transfer))

            pwr_2d.append(pwrspec_out[0])
            pwr_1d.append(pwrspec_out[1])

    if usecache_output_tag:
        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d,
                             usecache_output_tag, outdir=outdir)

        retval = (pwr_1d, pwr_1d_from_2d, pwr_2d)
    else:
        caller.multiprocess_stack()
        retval = None

    return retval


def batch_sim_run(simleft_key, simright_key,
                  weightleft_key, weightright_key,
                  inifile=None, datapath_db=None,
                  outdir="./plot_data_v2/",
                  usecache_output_tag=None, transfer=None):
    r"""
    typical weight matrix:
    db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv
    """

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(simleft_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if usecache_output_tag else True
    if generate:
        print "REGENERATING the power spectrum result cache: "

    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in mock_cases['realization']:
        dbkeydict = {}
        dbkeydict['map1_key'] = "%s:%s" % (simleft_key, index)
        dbkeydict['map2_key'] = "%s:%s" % (simright_key, index)
        dbkeydict['noiseinv1_key'] = weightleft_key
        dbkeydict['noiseinv2_key'] = weightright_key
        files = data_paths.convert_dbkeydict_to_filedict(dbkeydict,
                                                    datapath_db=datapath_db)

        pwrspec_out = caller.execute(files['map1_key'], files['map2_key'],
                                     files['noiseinv1_key'],
                                     files['noiseinv2_key'],
                                     inifile=inifile)

        if usecache_output_tag:
            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwrspec_out[0],
                                                  transfer=transfer))

            pwr_2d.append(pwrspec_out[0])
            pwr_1d.append(pwrspec_out[1])

    if usecache_output_tag:
        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d,
                             usecache_output_tag, outdir=outdir)
        retval = (pwr_1d, pwr_1d_from_2d, pwr_2d)
    else:
        caller.multiprocess_stack()
        retval = None

    return retval


def call_sim_autopower(basesims, treatments, weight, inifile=None,
                       inifile_phys=None, generate=False,
                       outdir="./plot_data_v2/"):
    r"""run all of the auto-power theory cases
    """
    datapath_db = data_paths.DataPath()

    for base in basesims:
        mapname = "%s_physical" % base

        usecache_output_tag = None
        if not generate:
            usecache_output_tag = mapname

        batch_physical_sim_run(mapname,
                               inifile=inifile_phys,
                               datapath_db=datapath_db,
                               outdir=outdir,
                               usecache_output_tag=usecache_output_tag)

        for treatment in treatments:
            mapname = base + treatment

            usecache_output_tag = None
            if not generate:
                usecache_output_tag = mapname

            batch_sim_run(mapname, mapname,
                          weight, weight, inifile=inifile,
                          datapath_db=datapath_db,
                          outdir=outdir,
                          usecache_output_tag=usecache_output_tag)
