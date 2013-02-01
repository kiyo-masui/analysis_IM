from utils import data_paths
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from utils import file_tools


def batch_wigglez_automock_run(mock_key, sel_key,
                               inifile=None, datapath_db=None,
                               output_tag=None):
    r"""TODO: make this work; wrote this but never really needed it yet
    """
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    mock_cases = datapath_db.fileset_cases(mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if output_tag else True
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

