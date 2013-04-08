"""Make simulated thermal noise realizations
"""
import numpy as np
import core.algebra as algebra
from numpy import random
import struct
from kiyopy import parse_ini
from utils import batch_handler

# given a noise file and total integration time, generate a noise map
# for areas with no integration time, use max_variance (in K)
params_init = {
               'output_file': "./thermal_noise.npy",
               'delta_temp_file': "./delta_temp_noise.npy",
               'weight_file': "./ok.npy",
               'total_integration': 100.,
               'max_stdev': np.inf,
               'seed': -1
               }
prefix = 'sst_'

class SimulateSingleThermal(object):
    r"""make a single simulation of thermal noise"""

    @batch_handler.log_timing
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

        self.output_file = self.params['output_file']
        self.delta_temp_file = self.params['delta_temp_file']
        self.total_integration = self.params['total_integration']
        self.weight_map = algebra.make_vect(
                                algebra.load(self.params['weight_file']))

        self.max_stdev = self.params['max_stdev']

        # set the random seed
        if (self.params['seed'] < 0):
            # The usual seed is not fine enough for parallel jobs
            randsource = open("/dev/random", "rb")
            self.seed = struct.unpack("I", randsource.read(4))[0]
            #self.seed = abs(long(outfile_physical.__hash__()))
        else:
            self.seed = self.params['seed']

        random.seed(self.seed)

    def noise_model(self, freq):
        r"""simple noise model: linear from 30K at 700 MHz to 20K at 800 MHz"""
        return 20 + (900. - freq)/200.*(10)

    @batch_handler.log_timing
    def execute(self, processes):
        freq_list = self.weight_map.get_axis("freq")
        delta_freq = np.roll(freq_list, 1) - freq_list
        bandwidth = delta_freq[1]
        # this may be too picky
        assert bandwidth == np.mean(delta_freq[1:]), "bad freq. axis"
        print "assuming uniform bandwidth: ", bandwidth

        integration_time = self.total_integration * 3600.

        freq_weight = np.apply_over_axes(np.sum, self.weight_map, [1, 2])
        best_freq = np.argmax(freq_weight.flatten())
        print "best freq: ", freq_list[best_freq]

        # normalize such that sum of the best slice is the total integration
        # time
        norm = integration_time / np.sum(self.weight_map[best_freq, :, :])
        self.weight_map *= norm

        sys_temp = self.noise_model(freq_list[best_freq]/1.e6)
        print "system temp at best: ", sys_temp

        #print self.weight_map[best_freq, :, :]
        delta_temp_map = sys_temp / np.sqrt(bandwidth * self.weight_map)
        delta_temp_map[np.isinf(delta_temp_map)] = self.max_stdev

        # print the delta T for the best slice in mK
        #print delta_temp_map[best_freq, :, :] * 1000.

        algebra.save(self.delta_temp_file, delta_temp_map)

        thermal_noise = np.random.randn(*delta_temp_map.shape)
        thermal_noise *= delta_temp_map
        thermal_noise = algebra.make_vect(thermal_noise, axis_names=('freq', 'ra', 'dec'))
        thermal_noise.copy_axis_info(delta_temp_map)

        #print thermal_noise[best_freq, :, :] * 1000.

        algebra.save(self.output_file, thermal_noise)

