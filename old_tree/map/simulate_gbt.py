"""Make a fake set of data that look like a real template but are thermal noise
and simulated signal"""
import numpy as np
import core.algebra as algebra
from numpy import random
import struct
from kiyopy import parse_ini
from utils import batch_handler
from utils import data_paths
from map import simulate_thermal
from map import simulate_gbt_signal
import os

params_init = {
               'template_key': "database to map",
               'output_key': "database to out map",
               'tack_on': None,
               'total_integration': 100.,
               'scenario': 'str',
               'refinement': 2.,
               'multiplier': 1.,
               'seed': -1
               }
prefix = 'sgbt_'


class SimulateGbt(object):
    r"""make a fake dataset based on a real one"""

    @batch_handler.log_timing
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

        self.template_key = self.params['template_key']
        self.output_key = self.params['output_key']
        self.total_integration = self.params['total_integration']
        self.scenario = self.params['scenario']
        self.refinement = self.params['refinement']
        self.multiplier = self.params['multiplier']
        self.tack_on = self.params['tack_on']

        # set the random seed
        if (self.params['seed'] < 0):
            print "no seed given; generating one (are you sure?)"
            # The usual seed is not fine enough for parallel jobs
            randsource = open("/dev/random", "rb")
            self.seed = struct.unpack("I", randsource.read(4))[0]
            #self.seed = abs(long(outfile_physical.__hash__()))
        else:
            self.seed = self.params['seed']

        random.seed(self.seed)

        self.datapath_db = data_paths.DataPath()

        self.input_weight_maps = self.return_maplist(self.template_key,
                                                     "noise_weight")

        self.output_weight_maps = self.return_maplist(self.output_key,
                                                      "noise_weight",
                                                      tack_on=self.tack_on)

        self.output_maps = self.return_maplist(self.output_key,
                                               "clean_map",
                                               tack_on=self.tack_on)

        self.output_delta_thermal = []
        self.output_thermal = []
        for mapfile in self.output_maps:
            basename = os.path.splitext(mapfile)[0]
            self.output_delta_thermal.append(basename + "_deltat.npy")
            self.output_thermal.append(basename + "_thermal.npy")

        self.output_signal = "gaussian_signal_simulation.npy"

        print "input weight maps: ", self.input_weight_maps
        print "output weight maps: ", self.output_weight_maps
        self.output_root = os.path.dirname(self.output_weight_maps[0])
        self.output_root += "/"
        print "output directory: ", self.output_root
        if not os.path.isdir(self.output_root):
            os.mkdir(self.output_root)


    def return_maplist(self, db_key, map_type, tack_on=None,
                       ignore=['firstpass']):
        template_flist = self.datapath_db.fetch(db_key, tack_on=tack_on)
        map_combinations = data_paths.unpack_cases(template_flist[0],
                                                   "sec;map_type")

        print map_combinations
        assert map_type in map_combinations['map_type'], \
               "no weight maps!"

        maplist = []
        seclist = map_combinations['sec']
        if ignore is not None:
            seclist = [tag for tag in seclist if tag not in ignore]

        seclist.sort()

        for sec in seclist:
            mapname = "%s;%s" % (sec, map_type)
            maplist.append(template_flist[1][mapname])

        return maplist

    def execute(self, processes):
        self.execute_thermalsim()
        self.execute_signalsim()
        self.execute_assembledir()

    @batch_handler.log_timing
    def execute_thermalsim(self):
        # note that we want each section to have a different seed
        # starting from the base random seed and incrementing
        seed_inc = self.seed
        for (weight_file, delta_file, thermal_file) in \
                zip(self.input_weight_maps,
                    self.output_delta_thermal,
                    self.output_thermal):

            params_init = {
               'output_file': thermal_file,
               'delta_temp_file': delta_file,
               'weight_file': weight_file,
               'total_integration': self.total_integration,
               'max_stdev': np.inf,
               'seed': seed_inc
               }

            print "parameters for thermal noise sim: ", params_init
            simthermal = simulate_thermal.SimulateSingleThermal(
                                                params_dict=params_init)
            simthermal.execute(0)
            seed_inc += 1

    @batch_handler.log_timing
    def execute_signalsim(self):
        params_init = {
               'output_root': self.output_root,
               'template_file': self.input_weight_maps[0],
               'outfile_physical': None,
               'outfile_raw': None,
               'outfile_delta': None,
               'outfile_beam': self.output_signal,
               'outfile_meansub': None,
               'outfile_degrade': None,
               'scenario': self.scenario,
               'seed': self.seed,
               'refinement': self.refinement,
               'weightfile': self.input_weight_maps[0]
               }
        print "parameters for signal sim: ", params_init
        sim_signal = simulate_gbt_signal.SimulateGbtSignal(
                                                params_dict=params_init)

        sim_signal.execute(0)

    @batch_handler.log_timing
    def execute_assembledir(self):
        # link the weights through to the simulation directory
        for (weight_file_in, weight_file_out) in \
                zip(self.input_weight_maps, self.output_weight_maps):
            os.symlink(weight_file_in, weight_file_out)
            os.symlink(weight_file_in + ".meta", weight_file_out + ".meta")

        signalfile = self.output_root + self.output_signal
        signalmap = algebra.make_vect(algebra.load(signalfile))
        signalmap *= self.multiplier

        # now load the signal simulation add thermal noise and save
        for (thermal_file, mapfile) in \
                zip(self.output_thermal, self.output_maps):
            thermalmap = algebra.make_vect(algebra.load(thermal_file))
            algebra.save(mapfile, signalmap + thermalmap)

