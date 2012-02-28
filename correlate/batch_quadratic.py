r"""Code to run large batches of quadratic estimators on combinations of
data/sims
"""
import numpy as np
import math
from utils import data_paths
from correlate import map_pair as mp
from correlate import pwrspec_estimation as pe
from utils import batch_handler
from kiyopy import parse_ini


def call_xspec_run(map1_key, map2_key,
                   noiseinv1_key, noiseinv2_key,
                   inifile=None):
    r"""a free-standing function which calls the xspec analysis
    """
    params_init = {
        "unitless": True,
        "return_3d": False,
        "truncate": False,
        "window": None,
        "refinement": 2,
        "pad": 5,
        "order": 2,
        "freq_list": tuple(range(256)),
        "bins": [0.00765314, 2.49977141, 35]
                   }
    prefix = 'xs_'

    params = parse_ini.parse(inifile, params_init, prefix=prefix)
    if inifile is None:
        print "WARNING: no ini file for pwrspec estimation"

    # initialize and calculate the xspec
    simpair = mp.MapPair(map1_key, map2_key,
                         noiseinv1_key, noiseinv2_key,
                         params['freq_list'])

    bparam = params['bins']
    bins = np.logspace(math.log10(bparam[0]),
                       math.log10(bparam[1]),
                       num = bparam[2], endpoint=True)

    retval = simpair.pwrspec_summary(window=params['window'],
                                     unitless=params['unitless'],
                                     bins=bins,
                                     truncate=params['truncate'],
                                     refinement=params['refinement'],
                                     pad=params['pad'],
                                     order=params['order'],
                                     return_3d=params['return_3d'])

    return retval


def call_phys_space_run(cube1_file, cube2_file,
                        inifile=None):
    """Directly call the power spectral estimation on some physical vol"""
    params_init = {
        "unitless": True,
        "return_3d": False,
        "truncate": False,
        "window": "blackman",
        "bins": [0.00765314, 2.49977141, 35]
                   }
    prefix = 'xs_'

    params = parse_ini.parse(inifile, params_init, prefix=prefix)
    if inifile is None:
        print "WARNING: no ini file for pwrspec estimation"

    bparam = params['bins']
    bins = np.logspace(math.log10(bparam[0]),
                       math.log10(bparam[1]),
                       num = bparam[2], endpoint=True)

    retval = pe.calculate_xspec_file(cube1_file, cube2_file, bins,
                    weight1_file=None, weight2_file=None,
                    truncate=params['truncate'], window=params['window'],
                    return_3d=params['return_3d'], unitless=params['unitless'])

    return retval


def batch_physical_sim_run(sim_key, inifile=None):
    """Test the power spectral estimator using simulations"""
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath
    mock_cases = datapath_db.fileset_cases(sim_key, "realization")
    fileset = datapath_db.fileset_cases(sim_key, "realization")

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for index in mock_cases['realization']:
        map1_file = datapath_db.fetch("%s:%s" % (sim_key, index))
        map2_file = datapath_db.fetch("%s:%s" % (sim_key, index))

        caller.execute(map1_file, map2_file,
                       inifile=inifile)

    caller.multiprocess_stack()


def convert_keydict_to_filedict(dbkeydict, db=None):
    if db is None:
        db = data_paths.DataPath()

    filedict = {}
    for name in dbkeydict:
        filedict[name] = db.fetch(dbkeydict[name])

    return filedict


def batch_sim_run(simleft_key, simright_key,
                  weightleft_key, weightright_key,
                  inifile=None):
    r"""
    typical weight matrix:
    db:GBT_15hr_map_cleaned_0mode:A_with_B;noise_inv"""
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    mock_cases = datapath_db.fileset_cases(simleft_key, "realization")

    for index in mock_cases['realization']:
        input = {}
        input['map1_key'] = "%s:%s" % (simleft_key, index)
        input['map2_key'] = "%s:%s" % (simright_key, index)
        input['noiseinv1_key'] = "%s" % weightleft_key
        input['noiseinv2_key'] = "%s" % weightright_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    caller.multiprocess_stack()


def batch_GBTxwigglez_data_run(gbt_map_key, wigglez_map_key,
                               wigglez_mock_key, wigglez_selection_key,
                               inifile=None):
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    map_cases = datapath_db.fileset_cases(gbt_map_key, "type;treatment")
    mock_cases = datapath_db.fileset_cases(wigglez_mock_key, "realization")

    for treatment in map_cases['treatment']:
        input = {}
        input['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
        input['map2_key'] = "%s" % wigglez_map_key
        input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
        input['noiseinv2_key'] = "%s" % wigglez_selection_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    for treatment in map_cases['treatment']:
        for index in mock_cases['realization']:
            input = {}
            input['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
            input['map2_key'] = "%s:%d" % (wigglez_mock_key, index)
            input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
            input['noiseinv2_key'] = "%s" % wigglez_selection_key
            files = convert_keydict_to_filedict(input, db=datapath_db)

            caller.execute(files['map1_key'], files['map2_key'],
                           files['noiseinv1_key'], files['noiseinv2_key'],
                           inifile=inifile)

    caller.multiprocess_stack()


def batch_GBTxwigglez_trans_run(sim_key, sim_wigglez,
                                base_sim_GBT, gbt_map_key,
                                wigglez_selection_key,
                                inifile=None):
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    map_cases = datapath_db.fileset_cases(sim_key, "type;treatment")

    map1_key = "db:%s:0" % base_sim_GBT
    map2_key = "db:%s:0" % sim_wigglez
    noiseinv1_key = "db:%s:weight;0modes" % (gbt_map_key)
    noiseinv2_key = "db:%s" % wigglez_selection_key

    caller.execute(map1_key, map2_key,
                   noiseinv1_key, noiseinv2_key,
                   inifile=inifile)

    for treatment in map_cases['treatment']:
        map1_key = "db:%s:map;%s" % (sim_key, treatment)
        map2_key = "db:%s:0" % sim_wigglez
        noiseinv1_key = "db:%s:weight;%s" % (gbt_map_key, treatment)
        noiseinv2_key = "db:%s" % wigglez_selection_key

        caller.execute(map1_key, map2_key,
                       noiseinv1_key, noiseinv2_key,
                       inifile=inifile)

    caller.multiprocess_stack()


def batch_one_sided_trans_run(modeloss_simkey, sim_key,
                              modeloss_weight_root,
                              inifile=None):
    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    # TODO: this is not a generic loop based on new db structure
    map1_key = "db:%s:0" % sim_key
    map2_key = "db:%s:0" % sim_key
    noiseinv1_key = "db:%s:weight;0modes" % (modeloss_weight_root)
    noiseinv2_key = "db:%s:weight;0modes" % (modeloss_weight_root)

    caller.execute(map1_key, map2_key,
                   noiseinv1_key, noiseinv2_key,
                   inifile=inifile)

    map_cases = datapath_db.fileset_cases(sim_key, "type;treatment")
    for treatment in map_cases['treatment']:
        map1_key = "db:%s:map;%s" % (modeloss_simkey, treatment)
        map2_key = "db:%s:0" % sim_key

        noiseinv1_key = "db:%s:weight;%s" % (modeloss_weight_root, \
                                                 treatment)

        noiseinv2_key = "db:%s:ones;%s" % (modeloss_weight_root, \
                                               treatment)

        caller.execute(map1_key, map2_key,
                       noiseinv1_key, noiseinv2_key,
                       inifile=inifile)

    caller.multiprocess_stack()


def batch_wigglez_automock_run(mock_key, sel_key,
                               inifile=None):

    datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    mock_cases = datapath_db.fileset_cases(mock_key, "realization")
    for index in mock_cases['realization']:
        map1_key = "db:%s:%s" % (mock_key, index)
        map2_key = "db:%s:%s" % (mock_key, index)
        noiseinv1_key = "db:%s" % sel_key
        noiseinv2_key = "db:%s" % sel_key

        caller.execute(map1_key, map2_key,
                       noiseinv1_key, noiseinv2_key,
                       inifile=inifile)

    caller.multiprocess_stack()


def batch_data_run(map_key, inifile=None):
    datapath_db = data_paths.DataPath()
    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    map_cases = datapath_db.fileset_cases(map_key, "pair;type;treatment")
    for treatment in map_cases['treatment']:
        uniq_pairs = data_paths.GBTauto_cross_pairs(map_cases['pair'],
                                                    map_cases['pair'],
                                                    cross_sym="_with_")

        for item in uniq_pairs:
            map1_key = "db:%s:%s;map;%s" % (map_key, item[0], treatment)
            map2_key = "db:%s:%s;map;%s" % (map_key, item[1], treatment)
            noiseinv1_key = "db:%s:%s;weight;%s" % (map_key, item[0], treatment)
            noiseinv2_key = "db:%s:%s;weight;%s" % (map_key, item[1], treatment)

            caller.execute(map1_key, map2_key,
                           noiseinv1_key, noiseinv2_key,
                           inifile=inifile)

    caller.multiprocess_stack()
