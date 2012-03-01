import numpy as np
import math
from utils import data_paths
from utils import batch_handler
from kiyopy import parse_ini
from correlate import batch_quadratic as bq
from correlate import pwrspec_estimation as pe

def batch_physical_sim_run(sim_key, inifile=None, datapath_db=None,
                           compile_tag=None):
    """Test the power spectral estimator using simulations"""
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(sim_key, "realization")

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    generate = False if compile_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in mock_cases['realization']:
        map1_file = datapath_db.fetch("%s:%s" % (sim_key, index))
        map2_file = datapath_db.fetch("%s:%s" % (sim_key, index))

        pwr2d_run, pwr1d_run = caller.execute(map1_file, map2_file,
                                              inifile=inifile)

        if compile_tag:
            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                                  transfer=transfer))

            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)

    if compile_tag:
        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, compile_tag,
                             outdir="./plot_data_v2")

        retval = (pwr_1d, pwr_1d_from_2d, pwr_2d)
    else:
        caller.multiprocess_stack()
        retval = None

    return retval


def batch_sim_run(simleft_key, simright_key,
                  weightleft_key, weightright_key,
                  inifile=None, datapath_db=None,
                  compile_tag=None):
    r"""
    typical weight matrix:
    db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv
    """

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(simleft_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if compile_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in mock_cases['realization']:
        input = {}
        input['map1_key'] = "%s:%s" % (simleft_key, index)
        input['map2_key'] = "%s:%s" % (simright_key, index)
        input['noiseinv1_key'] = weightleft_key
        input['noiseinv2_key'] = weightright_key
        files = bq.convert_keydict_to_filedict(input, db=datapath_db)

        pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                              files['map2_key'],
                                              files['noiseinv1_key'],
                                              files['noiseinv2_key'],
                                              inifile=inifile)
        if compile_tag:
            pwr_1d_from_2d.append(pe.convert_2d_to_1d(pwr2d_run,
                                                  transfer=transfer))

            pwr_2d.append(pwr2d_run)
            pwr_1d.append(pwr1d_run)


    if compile_tag:
        pe.summarize_pwrspec(pwr_1d, pwr_1d_from_2d, pwr_2d, compile_tag,
                             outdir="./plot_data_v2")
        retval = (pwr_1d, pwr_1d_from_2d, pwr_2d)
    else:
        caller.multiprocess_stack()
        retval = None

    return retval
