"""Make simulated signal realizations in the GBT survey volume
Notes on foreground calls for the future (also copy param member func.):
    from simulations import foregroundsck
    from simulations import pointsource
    syn = foregroundsck.Synchrotron()
    synfield = syn.getfield()
    ps = pointsource.DiMatteo()
    psf = ps.getfield()
"""
import os
import copy
import numpy as np
from simulations import corr21cm
import core.algebra as algebra
import map.beam as beam
from numpy import random
import struct
from kiyopy import parse_ini
from utils import units
from utils import batch_handler
from mpi4py import MPI


params_init = {
               'output_root': "./",
               'template_file': "ok.npy",
               'outfile_physical': None,
               'outfile_raw': None,
               'outfile_delta': None,
               'outfile_beam': None,
               'outfile_meansub': None,
               'outfile_degrade': None,
               'degrade_factor' : 1.1,
               'scenario': 'str',
               'seed': -1,
               'refinement': 2,
               'weightfile': None,
               'simnum' : None,
               }
prefix = 'sg_'

class SimulateGbt(object):
    r"""Class to handle signal-only sim ini files"""

    @batch_handler.log_timing
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

        if not os.path.isdir(self.params['output_root']):
            os.mkdir(self.params['output_root'])

        self.refinement = self.params['refinement']
        self.scenario = self.params['scenario']
        self.template_file = self.params['template_file']
        self.output_root = self.params['output_root']
        # here we use 300 h km/s from WiggleZ for streaming dispersion
        self.streaming_dispersion = 300.*0.72

        self.template_map = algebra.make_vect(
                                algebra.load(self.template_file))

        # determine the beam model
        self.beam_data = np.array([0.316148488246, 0.306805630985,
                                   0.293729620792, 0.281176247549,
                                   0.270856788455, 0.26745856078,
                                   0.258910010848, 0.249188429031])

        self.freq_data = np.array([695, 725, 755, 785, 815, 845, 875, 905],
                                 dtype=float)
        self.freq_data *= 1.0e6

        # set the random seed
        if (self.params['seed'] < 0):
            # The usual seed is not fine enough for parallel jobs
            randsource = open("/dev/random", "rb")
            self.seed = struct.unpack("I", randsource.read(4))[0]
            #self.seed = abs(long(outfile_physical.__hash__()))
        else:
            self.seed = self.params['seed']

        random.seed(self.seed)

        # register any maps that need to be produced
        self.sim_map_phys = None
        self.sim_map = None
        self.sim_map_delta = None
        self.sim_map_withbeam = None
        self.sim_map_meansub = None
        self.sim_map_degrade = None

    def mpiexecute(self, processes):
        r"""The MPI method """
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        params = self.params

        if params['simnum'] != None:
            n_sim = params['simnum']
        else: 
            print 'MPI need total simulation number. '
            n_sim = 1
            exit()

        outfile_raw_temp = params['outfile_raw']
        outfile_physical_temp = params['outfile_physical']
        outfile_delta_temp = params['outfile_delta']
        outfile_beam_temp = params['outfile_beam']
        outfile_meansub_temp = params['outfile_meansub']
        outfile_degrade_temp = params['outfile_degrade']

        comm.barrier()

        if rank<n_sim:
            sim_list = range(rank, n_sim, size)
            print "RANK %d : "%rank, sim_list
            print 
            for sim in sim_list:
                if outfile_raw_temp:
                    self.params['outfile_raw'] = outfile_raw_temp%sim

                if outfile_physical_temp:
                    self.params['outfile_physical'] = outfile_physical_temp%sim

                if outfile_delta_temp:
                    self.params['outfile_delta'] = outfile_delta_temp%sim

                if outfile_beam_temp:
                    self.params['outfile_beam'] = outfile_beam_temp%sim

                if outfile_meansub_temp:
                    self.params['outfile_meansub'] = outfile_meansub_temp%sim

                if outfile_degrade_temp:
                    self.params['outfile_degrade'] = outfile_degrade_temp%sim
                print "RANK %d : creat sim_%03d maps"%(rank, sim)
                print 
                self.execute(processes)
                print "RANK %d : creat sim_%03d maps finished"%(rank, sim)

        comm.barrier()


    @batch_handler.log_timing
    def execute(self, processes):
        # this generates the raw physical and observation space sims
        self.realize_simulation()

        if self.params['outfile_raw']:
            algebra.save(self.output_root + self.params['outfile_raw'],
                         self.sim_map)

        if self.params['outfile_physical']:
            algebra.save(self.output_root + self.params['outfile_physical'],
                         self.sim_map_phys)

        if self.params['outfile_delta']:
            self.make_delta_sim()
            algebra.save(self.output_root + self.params['outfile_delta'],
                         self.sim_map_delta)

        if self.params['outfile_beam']:
            self.convolve_by_beam()
            algebra.save(self.output_root + self.params['outfile_beam'],
                         self.sim_map_withbeam)

        if self.params['outfile_meansub']:
            self.subtract_mean()
            algebra.save(self.output_root + self.params['outfile_meansub'],
                         self.sim_map_meansub)

        if self.params['outfile_degrade']:
            self.degrade_to_common_res()
            algebra.save(self.output_root + self.params['outfile_degrade'],
                         self.sim_map_degrade)

    @batch_handler.log_timing
    def realize_simulation(self):
        """do basic handling to call Richard's simulation code
        this produces self.sim_map and self.sim_map_phys
        """
        if self.scenario == "nostr":
            print "running dd+vv and no streaming case"
            simobj = corr21cm.Corr21cm.like_kiyo_map(self.template_map)
            maps = simobj.get_kiyo_field_physical(refinement=self.refinement)

        else:
            if self.scenario == "str":
                print "running dd+vv and streaming simulation"
                simobj = corr21cm.Corr21cm.like_kiyo_map(self.template_map,
                                           sigma_v=self.streaming_dispersion)

                maps = simobj.get_kiyo_field_physical(refinement=self.refinement)

            if self.scenario == "ideal":
                print "running dd-only and no mean simulation"
                simobj = corr21cm.Corr21cm.like_kiyo_map(self.template_map)
                maps = simobj.get_kiyo_field_physical(
                                            refinement=self.refinement,
                                            density_only=True,
                                            no_mean=True,
                                            no_evolution=True)

        (gbtsim, gbtphys, physdim) = maps

        # process the physical-space map
        self.sim_map_phys = algebra.make_vect(gbtphys, axis_names=('freq', 'ra', 'dec'))
        pshp = self.sim_map_phys.shape

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

        self.sim_map_phys.info = info

        # process the map in observation coordinates
        self.sim_map = algebra.make_vect(gbtsim, axis_names=('freq', 'ra', 'dec'))
        self.sim_map.copy_axis_info(self.template_map)

    @batch_handler.log_timing
    def make_delta_sim(self):
        r"""this produces self.sim_map_delta"""
        freq_axis = self.sim_map.get_axis('freq') / 1.e6
        z_axis = units.nu21 / freq_axis - 1.0

        simobj = corr21cm.Corr21cm()
        T_b = simobj.T_b(z_axis) * 1e-3

        self.sim_map_delta = copy.deepcopy(self.sim_map)
        self.sim_map_delta /= T_b[:, np.newaxis, np.newaxis]

    @batch_handler.log_timing
    def convolve_by_beam(self):
        r"""this produces self.sim_map_withbeam"""
        beamobj = beam.GaussianBeam(self.beam_data, self.freq_data)
        self.sim_map_withbeam = beamobj.apply(self.sim_map)

    @batch_handler.log_timing
    def degrade_to_common_res(self):
        r"""this produces self.sim_map_degrade"""
        # this depends on having simulations with the means subtracted
        if self.sim_map_meansub is None:
            self.subtract_mean()

        degrade_factor = self.params['degrade_factor']
        print "degrading the resolution to a %2.1f common beam"%degrade_factor
        beam_diff = np.sqrt(max(degrade_factor*self.beam_data)**2-(self.beam_data)**2)
        common_resolution = beam.GaussianBeam(beam_diff, self.freq_data)
        # Convolve to a common resolution.
        self.sim_map_degrade = common_resolution.apply(self.sim_map_meansub)
        #sim_map[np.isinf(sim_map)] = 0.
        #sim_map[np.isnan(sim_map)] = 0.

    @batch_handler.log_timing
    def subtract_mean(self):
        r"""this produces self.sim_map_meansub"""
        # this depends on having simulations convolved by the beam
        if self.sim_map_withbeam is None:
            self.convolve_by_beam()

        self.sim_map_meansub = copy.deepcopy(self.sim_map_withbeam)
        print "sim meansub using: " + self.params['weightfile']
        noise_inv = algebra.make_vect(algebra.load(self.params['weightfile']))
        means = np.sum(np.sum(noise_inv * self.sim_map_meansub, -1), -1)
        means /= np.sum(np.sum(noise_inv, -1), -1)
        means.shape += (1, 1)
        self.sim_map_meansub -= means
        # the weights will be zero in some places
        self.sim_map_meansub[noise_inv < 1.e-20] = 0.
        #self.sim_map_meansub[np.isinf(self.sim_map)] = 0.
        #self.sim_map_meansub[np.isnan(self.sim_map)] = 0.
