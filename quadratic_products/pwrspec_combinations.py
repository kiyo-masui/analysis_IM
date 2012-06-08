import os
import math
import numpy as np
from utils import aggregate_outputs
from kiyopy import parse_ini
from quadratic_products import pwrspec_estimator as pe
from foreground_clean import map_pair as mp
from utils import data_paths as dp

def pwrspec_caller(map1_key, map2_key,
                   noiseinv1_key, noiseinv2_key,
                   params):

    mappair = mp.MapPair(map1_key, map2_key,
                         noiseinv1_key, noiseinv2_key,
                         params['freq_list'],
                         input_filenames=True)

    bparam = params['bins']
    bins = np.logspace(math.log10(bparam[0]),
                       math.log10(bparam[1]),
                       num=bparam[2], endpoint=True)

    if params["degrade_resolution"]:
        mappair.degrade_resolution()

    if params["factorizable_noise"]:
        mappair.make_noise_factorizable()

    if params["meansub"]:
        mappair.subtract_weighted_mean()

    retval = mappair.pwrspec_summary(window=params['window'],
                                     unitless=params['unitless'],
                                     bins=bins,
                                     truncate=params['truncate'],
                                     refinement=params['refinement'],
                                     pad=params['pad'],
                                     order=params['order'],
                                     return_3d=params['return_3d'])

    return retval


gbtdataautopower_init = {
        "map_key": "test_map",
        "outfile": "test_file.shelve",
        "unitless": True,
        "return_3d": False,
        "truncate": False,
        "window": None,
        "degrade_resolution": False,
        "factorizable_noise": False,
        "meansub": False,
        "refinement": 2,
        "pad": 5,
        "order": 2,
        "freq_list": tuple(range(256)),
        "bins": [0.00765314, 2.49977141, 35]
               }
gbtdataautopower_prefix = 'xs_'

class GbtDataAutopower(object):
    r"""Handler to call the power spectral estimation for different
    combinations of auto, cross-powers
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          gbtdataautopower_init,
                                          prefix=gbtdataautopower_prefix)

        self.freq_list = np.array(self.params['freq_list'], dtype=int)


    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations.pwrspec_caller"
        caller = aggregate_outputs.AggregateOutputs(funcname)

        map_key = self.params['map_key']
        map_cases = self.datapath_db.fileset_cases(map_key,
                                                   "pair;type;treatment")

        unique_pairs = dp.GBTauto_cross_pairs(map_cases['pair'],
                                              map_cases['pair'],
                                              cross_sym="_with_")

        for treatment in map_cases['treatment']:
            for item in unique_pairs:
                dbkeydict = {}
                mapset0 = (map_key, item[0], treatment)
                mapset1 = (map_key, item[1], treatment)
                dbkeydict['map1_key'] = "%s:%s;map;%s" % mapset0
                dbkeydict['map2_key'] = "%s:%s;map;%s" % mapset1
                dbkeydict['noiseinv1_key'] = "%s:%s;noise_inv;%s" % mapset0
                dbkeydict['noiseinv2_key'] = "%s:%s;noise_inv;%s" % mapset1
                files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                                         datapath_db=self.datapath_db)

                execute_key = "%s:%s" % (item[0], treatment)
                caller.execute(files['map1_key'],
                               files['map2_key'],
                               files['noiseinv1_key'],
                               files['noiseinv2_key'],
                               self.params,
                               execute_key=execute_key)


        caller.multiprocess_stack(self.params["outfile"], debug=False)


batchsimautopower_init = {
        "map_key": "test_map",
        "sim_file": "sim.npy",
        "outfile": "test_file.shelve",
        "unitless": True,
        "return_3d": False,
        "truncate": False,
        "window": None,
        "degrade_resolution": False,
        "factorizable_noise": False,
        "meansub": False,
        "refinement": 2,
        "pad": 5,
        "order": 2,
        "freq_list": tuple(range(256)),
        "bins": [0.00765314, 2.49977141, 35]
               }
batchsimautopower_prefix = 'bs_'

class BatchSimAutopower(object):
    r"""Handler to call the power spectral estimation for different
    combinations of auto, cross-powers
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          batchsimautopower_init,
                                          prefix=batchsimautopower_prefix)

        self.freq_list = np.array(self.params['freq_list'], dtype=int)


    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations.pwrspec_caller"
        caller = aggregate_outputs.AggregateOutputs(funcname)

        map_key = self.params['map_key']
        sim_file = self.params['sim_file']
        map_cases = self.datapath_db.fileset_cases(map_key,
                                                   "pair;type;treatment")

        unique_pairs = dp.GBTauto_cross_pairs(map_cases['pair'],
                                              map_cases['pair'],
                                              cross_sym="_with_")

        treatment = "0modes"

        for item in unique_pairs:
            dbkeydict = {}
            mapset0 = (map_key, item[0], treatment)
            mapset1 = (map_key, item[1], treatment)
            dbkeydict['map1_key'] = "%s:%s;map;%s" % mapset0
            dbkeydict['map2_key'] = "%s:%s;map;%s" % mapset1
            dbkeydict['noiseinv1_key'] = "%s:%s;noise_inv;%s" % mapset0
            dbkeydict['noiseinv2_key'] = "%s:%s;noise_inv;%s" % mapset1
            files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                                     datapath_db=self.datapath_db)

            execute_key = "%s:%s" % (item[0], treatment)
            caller.execute(sim_file,
                           sim_file,
                           files['noiseinv1_key'],
                           files['noiseinv2_key'],
                           self.params,
                           execute_key=execute_key)

        caller.multiprocess_stack(self.params["outfile"], debug=False)
