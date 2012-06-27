"""Make simulated signal realizations in the GBT survey volume"""
import scipy as sp
import numpy as np
from simulations import corr21cm
#from simulations import foregroundsck
#from simulations import pointsource
import core.algebra as algebra
import map.beam as beam
from core import constants as cc
import multiprocessing
from utils import data_paths
import sys
from numpy import random
import struct
from kiyopy import parse_ini
from utils import file_tools
from utils import units
from utils import batch_handler
import os


# TODO: confirm ra=x, dec=y (thetax = 5, thetay = 3 in 15hr)
def realize_simulation(template_map, scenario=None, seed=None, refinement=2):
    """do basic handling to call Richard's simulation code
    here we use 300 h km/s from WiggleZ for streaming dispersion
    Notes on foreground calls for the future (also copy param member func.):
    syn = foregroundsck.Synchrotron()
    synfield = syn.getfield()
    ps = pointsource.DiMatteo()
    psf = ps.getfield()
    """
    if seed is not None:
        random.seed(seed)

    if scenario == "nostr":
        print "running dd+vv and no streaming case"
        simobj = corr21cm.Corr21cm.like_kiyo_map(template_map)
        maps = simobj.get_kiyo_field_physical(refinement=refinement)

    else:
        if scenario == "str":
            print "running dd+vv and streaming simulation"
            simobj = corr21cm.Corr21cm.like_kiyo_map(template_map,
                                                     sigma_v=300.*0.72)
            maps = simobj.get_kiyo_field_physical(refinement=refinement)

        if scenario == "ideal":
            print "running dd-only and no mean simulation"
            simobj = corr21cm.Corr21cm.like_kiyo_map(template_map)
            maps = simobj.get_kiyo_field_physical(refinement=refinement,
                                density_only=True,
                                no_mean=True,
                                no_evolution=True)

    return maps


def make_simulation_set(template_file, output_root, outfile_physical=None,
                        outfile_raw=None, outfile_beam=None,
                        outfile_beam_plus_data=None,
                        verbose=True, scenario=None, seed=None,
                        refinement=2):
    """Produce simulated GBT data volumes of three types:
    (from dimensions of a given template file)
    0. the simulation in a physical volume
    1. the raw simulation in z
    2. the simulation convolved by the instrumental beam
    3. the simulation convolved by the instrumental beam plus the template
    """
    template_map = algebra.make_vect(algebra.load(template_file))

    # The usual seed is not fine enough for parallel jobs
    if not seed:
        randsource = open("/dev/random", "rb")
        seed = struct.unpack("I", randsource.read(4))[0]
        #seed = abs(long(outfile_physical.__hash__()))

    (gbtsim, gbtphys, physdim) = realize_simulation(template_map,
                                                    scenario=scenario,
                                                    seed=seed,
                                                    refinement=refinement)

    phys_map = algebra.make_vect(gbtphys, axis_names=('freq', 'ra', 'dec'))
    pshp = phys_map.shape

    # define the axes of the physical map; several alternatives are commented
    info = {}
    info['axes'] = ('freq', 'ra', 'dec')
    info['type'] = 'vect'
    info['freq_delta'] = abs(physdim[0] - physdim[1]) / float(pshp[0] - 1)
    info['freq_centre'] = physdim[0] + info['freq_delta'] * float(pshp[0] // 2)
    #        'freq_centre': abs(physdim[0] + physdim[1]) / 2.,

    info['ra_delta'] = abs(physdim[2]) / float(pshp[1] - 1)
    #info['ra_centre'] = info['ra_delta'] * float(pshp[1] // 2)
    #        'ra_centre': abs(physdim[2]) / 2.,
    info['ra_centre'] = 0.

    info['dec_delta'] = abs(physdim[3]) / float(pshp[2] - 1)
    #info['dec_centre'] = info['dec_delta'] * float(pshp[2] // 2)
    #        'dec_centre': abs(physdim[3]) / 2.,
    info['dec_centre'] = 0.

    phys_map.info = info

    sim_map = algebra.make_vect(gbtsim, axis_names=('freq', 'ra', 'dec'))
    sim_map.copy_axis_info(template_map)

    if outfile_raw:
        algebra.save(output_root + outfile_raw, sim_map)

    if outfile_physical:
        algebra.save(output_root + outfile_physical, phys_map)

    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    freq_data *= 1.0e6

    beamobj = beam.GaussianBeam(beam_data, freq_data)
    sim_map_withbeam = beamobj.apply(sim_map)

    if outfile_beam:
        algebra.save(output_root + outfile_beam, sim_map_withbeam)

    sim_map_withbeam += template_map
    if outfile_beam_plus_data:
        algebra.save(output_root + outfile_beam_plus_data, sim_map_withbeam)


def generate_delta_sim(input_file, output_file):
    r"""make the map with the temperature divided out (delta)"""
    print "reading %s -> %s (dividing by T_b(z))" % (input_file, output_file)

    simmap = algebra.make_vect(algebra.load(input_file))
    freq_axis = simmap.get_axis('freq') / 1.e6
    z_axis = units.nu21 / freq_axis - 1.0

    simobj = corr21cm.Corr21cm()
    T_b = simobj.T_b(z_axis)*1e-3

    simmap /= T_b[:, np.newaxis, np.newaxis]

    print "saving to" + output_file
    algebra.save(output_file, simmap)


def generate_proc_sim(input_file, weightfile, output_file,
                      meansub=False, degrade=False):
    r"""make the maps with various combinations of beam conv/meansub"""
    print "%s -> %s (beam, etc.)" % (input_file, output_file)
    simmap = algebra.make_vect(algebra.load(input_file))

    if degrade:
        print "performing common resolution convolution"
        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
        freq_data *= 1.0e6
        beam_diff = sp.sqrt(max(1.1 * beam_data) ** 2 - (beam_data) ** 2)
        common_resolution = beam.GaussianBeam(beam_diff, freq_data)
        # Convolve to a common resolution.
        simmap = common_resolution.apply(simmap)

    if meansub:
        print "performing mean subtraction"
        noise_inv = algebra.make_vect(algebra.load(weightfile))
        means = sp.sum(sp.sum(noise_inv * simmap, -1), -1)
        means /= sp.sum(sp.sum(noise_inv, -1), -1)
        means.shape += (1, 1)
        simmap -= means
        # the weights will be zero in some places
        simmap[noise_inv < 1.e-20] = 0.

    # extra sanity?
    simmap[np.isinf(simmap)] = 0.
    simmap[np.isnan(simmap)] = 0.

    print "saving to" + output_file
    algebra.save(output_file, simmap)


params_init = {
               'output_root': "./",
               'template_file': "ok.npy",
               'outfile_physical': "ok.npy",
               'outfile_raw': "ok.npy",
               'outfile_beam': "ok.npy",
               'outfile_beam_plus_data': "ok.npy",
               'scenario': 'str',
               'seed': -1,
               'refinement': 2
               }
prefix = 'sg_'

class SimulateGbt():
    r"""Class to handle signal-only sim ini files"""
    @batch_handler.log_timing
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

        if not os.path.isdir(self.params['output_root']):
            os.mkdir(self.params['output_root'])

    @batch_handler.log_timing
    def execute(self, processes):
        make_simulation_set(self.params['template_file'],
                            self.params['output_root'],
                            self.params['outfile_physical'],
                            self.params['outfile_raw'],
                            self.params['outfile_beam'],
                            self.params['outfile_beam_plus_data'],
                            verbose=True,
                            scenario=self.params['scenario'],
                            seed=self.params['seed'],
                            refinement=self.params['refinement'])
