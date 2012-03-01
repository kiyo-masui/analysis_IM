r"""Code to run large batches of quadratic estimators on combinations of
data/sims

The code is run in two steps:
1. with compile_tag=None (default) calculate and cache all the power spectra
2. with compile_tag=some identifier for the outputs, read the caches, write a
summary file, and return the power spectra and bins

The reason for this structure is that it is expensive to compute the hundreds
of power spectra that go into a plot, but one also needs to try out many
different ways of combining the outputs to form transfer functions, etc.,
different 2D->1D binnings, etc.

The bin edges for the output power spectra are in
        bin_left = pwr_1d[0]['bin_left']
        bin_center = pwr_1d[0]['bin_center']
        bin_right = pwr_1d[0]['bin_right']

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
                         params['freq_list'], avoid_db=True)

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


def batch_physical_sim_run(sim_key, inifile=None, datapath_db=None,
                           compile_tag=None):
    """Test the power spectral estimator using simulations"""
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    mock_cases = datapath_db.fileset_cases(sim_key, "realization")

    funcname = "correlate.batch_quadratic.call_phys_space_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

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


def convert_keydict_to_filedict(dbkeydict, db=None):
    if db is None:
        db = data_paths.DataPath()

    filedict = {}
    for name in dbkeydict:
        filedict[name] = db.fetch(dbkeydict[name])

    return filedict


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

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    mock_cases = datapath_db.fileset_cases(simleft_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    pwr_1d = []
    pwr_1d_from_2d = []
    pwr_2d = []
    for index in mock_cases['realization']:
        input = {}
        input['map1_key'] = "%s:%s" % (simleft_key, index)
        input['map2_key'] = "%s:%s" % (simright_key, index)
        input['noiseinv1_key'] = weightleft_key
        input['noiseinv2_key'] = weightright_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

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

def batch_GBTxwigglez_data_run(gbt_map_key, wigglez_map_key,
                               wigglez_mock_key, wigglez_selection_key,
                               inifile=None, datapath_db=None):
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    map_cases = datapath_db.fileset_cases(gbt_map_key, "type;treatment")
    mock_cases = datapath_db.fileset_cases(wigglez_mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for treatment in map_cases['treatment']:
        input = {}
        input['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
        input['map2_key'] = wigglez_map_key
        input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
        input['noiseinv2_key'] = wigglez_selection_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    for treatment in map_cases['treatment']:
        for index in mock_cases['realization']:
            input = {}
            input['map1_key'] = "%s:map;%s" % (gbt_map_key, treatment)
            input['map2_key'] = "%s:%s" % (wigglez_mock_key, index)
            input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
            input['noiseinv2_key'] = wigglez_selection_key
            files = convert_keydict_to_filedict(input, db=datapath_db)

            caller.execute(files['map1_key'], files['map2_key'],
                           files['noiseinv1_key'], files['noiseinv2_key'],
                           inifile=inifile)

    caller.multiprocess_stack()


def batch_GBTxwigglez_trans_run(sim_key, sim_wigglez,
                                base_sim_GBT, gbt_map_key,
                                wigglez_selection_key,
                                inifile=None, datapath_db=None):
    r"""
    Assume that the 0'th realization sim is used in the cleaning sims
    """
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    map_cases = datapath_db.fileset_cases(sim_key, "type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    input = {}
    input['map1_key'] = "%s:0" % base_sim_GBT
    input['map2_key'] = "%s:0" % sim_wigglez
    input['noiseinv1_key'] = "%s:weight;0modes" % (gbt_map_key)
    input['noiseinv2_key'] = wigglez_selection_key
    files = convert_keydict_to_filedict(input, db=datapath_db)

    caller.execute(files['map1_key'], files['map2_key'],
                   files['noiseinv1_key'], files['noiseinv2_key'],
                   inifile=inifile)

    for treatment in map_cases['treatment']:
        input = {}
        input['map1_key'] = "%s:map;%s" % (sim_key, treatment)
        input['map2_key'] = "%s:0" % sim_wigglez
        input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
        input['noiseinv2_key'] = wigglez_selection_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    caller.multiprocess_stack()


def batch_one_sided_trans_run(modeloss_simkey, sim_key,
                              modeloss_weight_root,
                              inifile=None, datapath_db=None):
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    map_cases = datapath_db.fileset_cases(sim_key, "type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    input = {}
    input['map1_key'] = "%s:0" % sim_key
    input['map2_key'] = "%s:0" % sim_key
    input['noiseinv1_key'] = "%s:weight;0modes" % (modeloss_weight_root)
    input['noiseinv2_key'] = "%s:weight;0modes" % (modeloss_weight_root)
    files = convert_keydict_to_filedict(input, db=datapath_db)

    caller.execute(files['map1_key'], files['map2_key'],
                   files['noiseinv1_key'], files['noiseinv2_key'],
                   inifile=inifile)

    for treatment in map_cases['treatment']:
        input = {}
        input['map1_key'] = "%s:map;%s" % (modeloss_simkey, treatment)
        input['map2_key'] = "%s:0" % sim_key

        input['noiseinv1_key'] = "%s:weight;%s" % (modeloss_weight_root, \
                                                   treatment)

        input['noiseinv2_key'] = "%s:ones;%s" % (modeloss_weight_root, \
                                                 treatment)
        files = convert_keydict_to_filedict(input, db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    caller.multiprocess_stack()


def batch_wigglez_automock_run(mock_key, sel_key,
                               inifile=None, datapath_db=None):

    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    mock_cases = datapath_db.fileset_cases(mock_key, "realization")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for index in mock_cases['realization']:
        input = {}
        input['map1_key'] = "%s:%s" % (mock_key, index)
        input['map2_key'] = "%s:%s" % (mock_key, index)
        input['noiseinv1_key'] = sel_key
        input['noiseinv2_key'] = sel_key
        files = convert_keydict_to_filedict(input, db=datapath_db)

        caller.execute(files['map1_key'], files['map2_key'],
                       files['noiseinv1_key'], files['noiseinv2_key'],
                       inifile=inifile)

    caller.multiprocess_stack()


def batch_data_run(map_key, inifile=None, datapath_db=None):
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    outpath = datapath_db.fetch("quadratic_batch_data")
    print "writing to: " + outpath

    map_cases = datapath_db.fileset_cases(map_key, "pair;type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    caller = batch_handler.MemoizeBatch(funcname, outpath,
                                        generate=True, verbose=True)

    for treatment in map_cases['treatment']:
        uniq_pairs = data_paths.GBTauto_cross_pairs(map_cases['pair'],
                                                    map_cases['pair'],
                                                    cross_sym="_with_")

        for item in uniq_pairs:
            input = {}
            input['map1_key'] = "%s:%s;map;%s" % (map_key, item[0], treatment)
            input['map2_key'] = "%s:%s;map;%s" % (map_key, item[1], treatment)
            input['noiseinv1_key'] = "%s:%s;weight;%s" % (map_key, item[0], treatment)
            input['noiseinv2_key'] = "%s:%s;weight;%s" % (map_key, item[1], treatment)
            files = convert_keydict_to_filedict(input, db=datapath_db)

            caller.execute(files['map1_key'], files['map2_key'],
                           files['noiseinv1_key'], files['noiseinv2_key'],
                           inifile=inifile)

    caller.multiprocess_stack()


def theory_power_spectrum(redshift_filekey, bin_centers,
                          unitless=True, fileout="theory.dat"):
    r"""TODO: make this work ... needs some minor changes
    probably move elsewhere
    """
    zfilename = datapath_db.fetch(redshift_filekey, intend_read=True,
                                      pick='1')

    zspace_cube = algebra.make_vect(algebra.load(zfilename))
    simobj = corr21cm.Corr21cm.like_kiyo_map(zspace_cube)
    pwrspec_input = simobj.get_pwrspec(bin_center)
    if unitless:
        pwrspec_input *= bin_center ** 3. / 2. / math.pi / math.pi

    outfile = open(fileout, "w")
    for specdata in zip(bin_left, bin_center, bin_right, pwrspec_input):
        outfile.write(("%10.15g " * 4 + "\n") % specdata)

    outfile.close()



