import numpy as np
import math
from utils import data_paths
from utils import batch_handler
from kiyopy import parse_ini
from correlate import batch_quadratic as bq
from correlate import pwrspec_estimation as pe


def batch_GBTxwigglez_trans_run(sim_key, sim_wigglez,
                                base_sim_GBT, gbt_map_key,
                                wigglez_selection_key,
                                inifile=None, datapath_db=None
                                compile_tag=None):
    r"""
    Assume that the 0'th realization sim is used in the cleaning sims
    """
    if datapath_db is None:
        datapath_db = data_paths.DataPath()

    cache_path = datapath_db.fetch("quadratic_batch_data")

    map_cases = datapath_db.fileset_cases(sim_key, "type;treatment")

    funcname = "correlate.batch_quadratic.call_xspec_run"
    generate = False if compile_tag else True
    caller = batch_handler.MemoizeBatch(funcname, cache_path,
                                        generate=generate, verbose=True)

    input = {}
    input['map1_key'] = "%s:0" % base_sim_GBT
    input['map2_key'] = "%s:0" % sim_wigglez
    input['noiseinv1_key'] = "%s:weight;0modes" % (gbt_map_key)
    input['noiseinv2_key'] = wigglez_selection_key
    files = bq.convert_keydict_to_filedict(input, db=datapath_db)

    pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                          files['map2_key'],
                                          files['noiseinv1_key'],
                                          files['noiseinv2_key'],
                                          inifile=inifile)

    if compile_tag:
        pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run)

        pwr_1d_zero = pwr1d_run['binavg']
        pwr_2d_zero = pwr2d_run['binavg']
        pwr_1d_from_2d_zero = pwr1d_from_2d['binavg']

    transfer_functions = {}
    for treatment in map_cases['treatment']:
        input = {}
        input['map1_key'] = "%s:map;%s" % (sim_key, treatment)
        input['map2_key'] = "%s:0" % sim_wigglez
        input['noiseinv1_key'] = "%s:weight;%s" % (gbt_map_key, treatment)
        input['noiseinv2_key'] = wigglez_selection_key
        files = bq.convert_keydict_to_filedict(input, db=datapath_db)

        pwr2d_run, pwr1d_run = caller.execute(files['map1_key'],
                                              files['map2_key'],
                                              files['noiseinv1_key'],
                                              files['noiseinv2_key'],
                                              inifile=inifile)

        if compile_tag:
            pwr1d_from_2d = pe.convert_2d_to_1d(pwr2d_run)

            pwr_1d = pwr1d_run['binavg']
            pwr_1d_from_2d = pwr1d_from_2d['binavg']
            trans1d_mode = pwr_1d/pwr_1d_zero
            trans1d_from2d_mode = pwr_1d_from_2d/pwr_1d_from_2d_zero
            trans2d_mode = pwr2d_run['binavg']/pwr_2d_zero

            transfer_functions[mode_num] = (trans1d_mode,
                                            trans1d_from2d_mode,
                                            trans2d_mode)

            # assume that they all have the same binning
            bin_left = pwr1d_run['bin_left']
            bin_center = pwr1d_run['bin_center']
            bin_right = pwr1d_run['bin_right']
            counts_histo = pwr1d_run['counts_histo']

            filename = "./plot_data/%s_%dmodes_trans.dat" % (tag, mode_num)
            outfile = open(filename, "w")
            for specdata in zip(bin_left, bin_center,
                            bin_right, counts_histo, pwr_1d, trans1d_mode,
                            trans1d_from2d_mode):
                outfile.write(("%10.15g " * 7 + "\n") % specdata)
            outfile.close()

    if compile_tag:
        return transfer_functions
    else:
        caller.multiprocess_stack()
        return None
