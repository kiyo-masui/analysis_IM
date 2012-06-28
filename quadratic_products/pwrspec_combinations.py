import os
import math
import numpy as np
import shutil
from utils import aggregate_outputs
from kiyopy import parse_ini
from quadratic_products import pwrspec_estimator as pe
from foreground_clean import map_pair as mp
from utils import data_paths as dp


def pwrspec_caller(map1_key, map2_key,
                   noiseinv1_key, noiseinv2_key,
                   params):
    r"""Call the cross-power estimator for the real data"""

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


def phys_pwrspec_caller(cube1_file, cube2_file, params):
    r"""Call the cross-power estimator on simulations in physical
    coordinates"""
    bparam = params['bins']
    bins = np.logspace(math.log10(bparam[0]),
                       math.log10(bparam[1]),
                       num=bparam[2], endpoint=True)

    retval = pe.calculate_xspec_file(cube1_file, cube2_file, bins,
                                     weight1_file=None, weight2_file=None,
                                     unitless=params['unitless'],
                                     window=params['window'],
                                     truncate=params['truncate'],
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
    r"""Handle GBT x GBT
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


gbtdatanoisepower_init = {
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
gbtdatanoisepower_prefix = 'ns_'

class GbtDataNoisePower(object):
    r"""Handle GBT x GBT; AxA, BxB, CxC, etc.
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          gbtdatanoisepower_init,
                                          prefix=gbtdatanoisepower_prefix)

        self.freq_list = np.array(self.params['freq_list'], dtype=int)


    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations.pwrspec_caller"
        caller = aggregate_outputs.AggregateOutputs(funcname)

        map_key = self.params['map_key']
        map_cases = self.datapath_db.fileset_cases(map_key,
                                                   "pair;type;treatment")

        # We use the same cleaned maps and the weights are the same for
        # various pairs, e.g. A_with_B is the same as A_with_C, etc. because
        # the mode cleaning does not impact the weighting functions
        noise_pairs = {"A_with_A": "A_with_B",
                       "B_with_B": "B_with_A",
                       "C_with_C": "C_with_A",
                       "D_with_D": "D_with_A"
                      }

        for treatment in map_cases['treatment']:
            for item in noise_pairs:
                dbkeydict = {}
                mapset0 = (map_key, noise_pairs[item], treatment)
                mapset1 = (map_key, noise_pairs[item], treatment)
                dbkeydict['map1_key'] = "%s:%s;map;%s" % mapset0
                dbkeydict['map2_key'] = "%s:%s;map;%s" % mapset1
                dbkeydict['noiseinv1_key'] = "%s:%s;noise_inv;%s" % mapset0
                dbkeydict['noiseinv2_key'] = "%s:%s;noise_inv;%s" % mapset1
                files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                                         datapath_db=self.datapath_db)

                execute_key = "%s:%s" % (item, treatment)
                caller.execute(files['map1_key'],
                               files['map2_key'],
                               files['noiseinv1_key'],
                               files['noiseinv2_key'],
                               self.params,
                               execute_key=execute_key)


        caller.multiprocess_stack(self.params["outfile"], debug=False)


crosspower_init = {
        "map_key": "test_map",
        "wigglez_key": "WiggleZ_15hr_binned_data",
        "wigglez_sel_key": "WiggleZ_15hr_montecarlo",
        "wigglez_mock_key": "WiggleZ_15hr_mock",
        "outfile_mock": "test_file.shelve",
        "outfile_data": "test_file.shelve",
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
crosspower_prefix = 'wxs_'

class WiggleZxGBT(object):
    r"""Handle GBT x WiggleZ
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          crosspower_init,
                                          prefix=crosspower_prefix)

        print self.params
        self.freq_list = np.array(self.params['freq_list'], dtype=int)


    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations.pwrspec_caller"
        caller_data = aggregate_outputs.AggregateOutputs(funcname)
        caller_mock = aggregate_outputs.AggregateOutputs(funcname)

        wigglez_key = self.params['wigglez_key']
        sel_key = self.params['wigglez_sel_key']
        mock_key = self.params['wigglez_mock_key']
        mock_files = self.datapath_db.fetch(mock_key)

        map_key = self.params['map_key']
        map_cases = self.datapath_db.fileset_cases(map_key,
                                                   "type;treatment")

        for treatment in map_cases['treatment']:
            dbkeydict = {}
            dbkeydict['map1_key'] = "%s:map;%s" % (map_key, treatment)
            dbkeydict['map2_key'] = wigglez_key

            dbkeydict['noiseinv1_key'] = "%s:weight;%s" % \
                                             (map_key, treatment)

            dbkeydict['noiseinv2_key'] = sel_key
            files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                        datapath_db=self.datapath_db)

            execute_key = "data:%s" % treatment
            #print files, execute_key

            caller_data.execute(files['map1_key'],
                           files['map2_key'],
                           files['noiseinv1_key'],
                           files['noiseinv2_key'],
                           self.params,
                           execute_key=execute_key)

        caller_data.multiprocess_stack(self.params["outfile_data"],
                                       debug=False, ncpu=18)

        for treatment in map_cases['treatment']:
            for item in mock_files[0]:
                dbkeydict = {}
                dbkeydict['map1_key'] = "%s:map;%s" % (map_key, treatment)
                dbkeydict['map2_key'] = "%s:%s" % (mock_key, item)

                dbkeydict['noiseinv1_key'] = "%s:weight;%s" % \
                                             (map_key, treatment)

                dbkeydict['noiseinv2_key'] = sel_key
                files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                                         datapath_db=self.datapath_db)
                execute_key = "mock%s:%s" % (item, treatment)
                #print files, execute_key

                caller_mock.execute(files['map1_key'],
                               files['map2_key'],
                               files['noiseinv1_key'],
                               files['noiseinv2_key'],
                               self.params,
                               execute_key=execute_key)


        caller_mock.multiprocess_stack(self.params["outfile_mock"], debug=False, ncpu=18)


crosspowersim_init = {
        "map_key": "test_map",
        "tack_on": None,
        "wigglez_sim_file": "WiggleZ_15hr_binned_data",
        "wigglez_sel_key": "WiggleZ_15hr_montecarlo",
        "outfile_data": "test_file.shelve",
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
crosspowersim_prefix = 'wxss_'

class WiggleZxGBT_modesim(object):
    r"""Handle GBT x WiggleZ mode subtraction simulations
    TODO: consider merging this with WiggleZxGBT
    NOTE: in this case, it is better to parallelize over pipeline instance
    (across many sims) rather than just parallelize the spectral est. -- here
    ncpu=1 by default.
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          crosspowersim_init,
                                          prefix=crosspowersim_prefix)

        print self.params
        self.freq_list = np.array(self.params['freq_list'], dtype=int)


    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations.pwrspec_caller"
        caller_data = aggregate_outputs.AggregateOutputs(funcname)

        wigglez_sim_file = self.params['wigglez_sim_file']
        sel_key = self.params['wigglez_sel_key']
        map_key = self.params['map_key']
        tack_on = self.params['tack_on']

        map_cases = self.datapath_db.fileset_cases(map_key,
                                                   "type;treatment")

        for treatment in map_cases['treatment']:
            sim_mapkey = "%s:map;%s" % (map_key, treatment)
            sim_mapfile = self.datapath_db.fetch(sim_mapkey, tack_on=tack_on)
            sim_noisekey = "%s:weight;%s" % (map_key, treatment)
            sim_noisefile = self.datapath_db.fetch(sim_noisekey, tack_on=tack_on)
            wigglez_selfile = self.datapath_db.fetch(sel_key)

            execute_key = "sim:%s" % treatment

            caller_data.execute(sim_mapfile,
                                wigglez_sim_file,
                                sim_noisefile,
                                wigglez_selfile,
                                self.params,
                                execute_key=execute_key)

        caller_data.multiprocess_stack(self.params["outfile_data"],
                                       debug=False, ncpu=1)


# this does not actually do any batch processing, but just wraps a single
# quadratic estimator output in the same packaging/pipeline interaction
batchsimcrosspower_init = {
        "map_key": "test_map",
        "sim_file": "sim.npy",
        "wigglez_sim_file": "sim.npy",
        "wigglez_sel_key": "WiggleZ_15hr_montecarlo",
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
        "order": 1,
        "freq_list": tuple(range(256)),
        "bins": [0.00765314, 2.49977141, 35]
               }
batchsimcrosspower_prefix = 'bxs_'

class BatchSimCrosspower(object):
    r"""Take cross-power relevant for finding the cross beam transfer function
    """
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          batchsimcrosspower_init,
                                          prefix=batchsimcrosspower_prefix)

        self.freq_list = np.array(self.params['freq_list'], dtype=int)


    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations.pwrspec_caller"
        caller = aggregate_outputs.AggregateOutputs(funcname)

        dbkeydict = {}
        dbkeydict['noiseinv1_key'] = "%s:weight;0modes" % \
                                     self.params['map_key']

        dbkeydict['noiseinv2_key'] = self.params['wigglez_sel_key']
        files = dp.convert_dbkeydict_to_filedict(dbkeydict,
                                                 datapath_db=self.datapath_db)

        execute_key = "sim:0modes"
        caller.execute(self.params['sim_file'],
                       self.params['wigglez_sim_file'],
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

batchphysicalsim_init = {
        "sim_key": "sim_stuff",
        "outfile": "test_file.shelve",
        "unitless": True,
        "return_3d": False,
        "truncate": False,
        "window": None,
        "bins": [0.00765314, 2.49977141, 35]
               }
batchphysicalsim_prefix = 'bps_'

class BatchPhysicalSim(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          batchphysicalsim_init,
                                          prefix=batchphysicalsim_prefix)

    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations"
        funcname += ".phys_pwrspec_caller"
        caller = aggregate_outputs.AggregateOutputs(funcname)

        sim_files = self.datapath_db.fetch(self.params['sim_key'])

        for index in sim_files[0]:
            execute_key = "%s:phys" % index
            caller.execute(sim_files[1][index],
                           sim_files[1][index],
                           self.params,
                           execute_key=execute_key)

        caller.multiprocess_stack(self.params["outfile"], debug=False)

# this does not actually do any batch processing, but just wraps a single
# quadratic estimator output in the same packaging/pipeline interaction
singlephysicalsim_init = {
        "sim_file_left": "sim_stuff",
        "sim_file_right": "sim_stuff",
        "outfile": "test_file.shelve",
        "unitless": True,
        "return_3d": False,
        "truncate": False,
        "window": None,
        "bins": [0.00765314, 2.49977141, 35]
               }
singlephysicalsim_prefix = 'sps_'

class SinglePhysicalSim(object):
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          singlephysicalsim_init,
                                          prefix=singlephysicalsim_prefix)

    def execute(self, processes):
        funcname = "quadratic_products.pwrspec_combinations"
        funcname += ".phys_pwrspec_caller"
        caller = aggregate_outputs.AggregateOutputs(funcname)

        execute_key = "sim:phys"
        caller.execute(self.params['sim_file_left'],
                       self.params['sim_file_right'],
                       self.params, execute_key=execute_key)

        caller.multiprocess_stack(self.params["outfile"], debug=False)


cleanup_init = {
    "path_key": "ok",
    "tack_on": None
    }

cleanup_prefix = "clean_"

class CleanupCleanedMaps(object):
    r"""delete cleaned maps to save space (~4 GB each)"""
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        self.datapath_db = dp.DataPath()

        if parameter_file:
            self.params = parse_ini.parse(parameter_file,
                                          cleanup_init,
                                          prefix=cleanup_prefix)

    def execute(self, processes):
        directory = self.datapath_db.fetch(self.params['path_key'],
                                           tack_on=self.params['tack_on'])
        print ">>> DELETE DIRECTORY: %s" % directory
        try:
            shutil.rmtree(directory)
        except:
            print "This directory does not exist"
