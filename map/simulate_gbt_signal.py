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
from utils import data_paths
from mpi4py import MPI


params_init = {
               'output_root': "./",
               'template_file': "ok.npy",
               'outfile_physical': None,
               'outfile_raw': None,
               'outfile_delta': None,
               'outfile_optsim': None,
               'outfile_beam': None,
               'outfile_meansub': None,
               'outfile_degrade': None,
               'scenario': 'str',
               'seed': -1,
               'seed_file': '',
               'refinement': 2,
               'selection_file': None,
               'optcatalog_file': None,
               'weightfile': None,
               'simnum' : None
               }
prefix = 'sg_'

class SimulateGbtSignal(object):
    r"""Class to handle signal-only sim ini files"""

    #@batch_handler.log_timing
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = params_dict
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)


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
        outfile_optsim_temp = params['outfile_optsim']
        outfile_beam_temp = params['outfile_beam']
        outfile_meansub_temp = params['outfile_meansub']
        outfile_degrade_temp = params['outfile_degrade']

        seed_list = np.loadtxt(params['seed_file']).astype('int')

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

                if outfile_optsim_temp:
                    self.params['outfile_optsim'] = outfile_optsim_temp%sim

                if outfile_beam_temp:
                    self.params['outfile_beam'] = outfile_beam_temp%sim

                if outfile_meansub_temp:
                    self.params['outfile_meansub'] = outfile_meansub_temp%sim

                if outfile_degrade_temp:
                    self.params['outfile_degrade'] = outfile_degrade_temp%sim
                print "RANK %d : creat sim_%03d maps"%(rank, sim)
                print 

                ## set the random seed
                self.seed = seed_list[sim]
                random.seed(self.seed)

                self.execute(processes, rank)
                print "RANK %d : creat sim_%03d maps finished"%(rank, sim)

        comm.barrier()


    #@batch_handler.log_timing
    def execute(self, processes, rank=0):

        #if not os.path.isdir(self.params['output_root']):
        if not os.path.exists(self.params['output_root']):
            os.mkdir(self.params['output_root'])

        self.refinement = self.params['refinement']
        self.scenario = self.params['scenario']
        self.template_file = self.params['template_file']
        self.output_root = self.params['output_root']
        # here we use 300 h km/s from WiggleZ for streaming dispersion
        self.streaming_dispersion = 300.*0.72

        self.template_map = algebra.make_vect( algebra.load(self.template_file))
        #self.datapath_db = data_paths.DataPath()
        #self.template_map = self.datapath_db.fetch_multi(self.template_file)

        # determine the beam model
        self.beam_data = np.array([0.316148488246, 0.306805630985,
                                   0.293729620792, 0.281176247549,
                                   0.270856788455, 0.26745856078,
                                   0.258910010848, 0.249188429031])

        self.freq_data = np.array([695, 725, 755, 785, 815, 845, 875, 905],
                                 dtype=float)
        self.freq_data *= 1.0e6

        # register any maps that need to be produced
        self.sim_map_phys = None
        self.sim_map = None
        self.sim_map_delta = None
        self.sim_map_optsim = None
        self.sim_map_withbeam = None
        self.sim_map_meansub = None
        self.sim_map_degrade = None

        ## set the random seed
        #if (self.params['seed'] < 0):
        #    # The usual seed is not fine enough for parallel jobs
        #    randsource = open("/dev/random", "rb")
        #    self.seed = struct.unpack("I", randsource.read(4))[0]
        #    randsource.close()
        #    #self.seed = abs(long(outfile_physical.__hash__()))
        #else:
        #    self.seed = self.params['seed']
        #print "Rank %d: set random seed: "%rank, self.seed
        #random.seed(self.seed)

        # this generates the raw physical and observation space sims
        self.realize_simulation(rank)

        if self.params['outfile_raw']:
            filename = self.output_root + self.params['outfile_raw']
            print "saving raw sim to ", filename
            algebra.save(filename, self.sim_map)

        if self.params['outfile_physical']:
            filename = self.output_root + self.params['outfile_physical']
            print "saving physical-space sim to ", filename
            algebra.save(filename, self.sim_map_phys)

        if self.params['outfile_delta']:
            self.make_delta_sim()

            filename = self.output_root + self.params['outfile_delta']
            print "saving overdensity sim to ", filename
            algebra.save(filename, self.sim_map_delta)

        if self.params['outfile_optsim']:
            self.make_opt_sim()

            filename = self.output_root + self.params['outfile_optsim']
            print "saving optical galaxy sim to ", filename
            algebra.save(filename, self.sim_map_optsim)

        if self.params['outfile_beam']:
            self.convolve_by_beam()

            filename = self.output_root + self.params['outfile_beam']
            print "saving beam-convolved sim to ", filename
            algebra.save(filename, self.sim_map_withbeam)

        if self.params['outfile_meansub']:
            self.subtract_mean()

            filename = self.output_root + self.params['outfile_meansub']
            print "saving beam-convolved, meansub sim to ", filename
            algebra.save(filename, self.sim_map_meansub)

        if self.params['outfile_degrade']:
            self.degrade_to_common_res()

            filename = self.output_root + self.params['outfile_degrade']
            print "saving beam-convolved, meansub, degrade sim to ", filename
            algebra.save(filename, self.sim_map_degrade)

    #@batch_handler.log_timing
    def realize_simulation(self, rank=0):
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

        print "Rank %d"%rank + "set scenario: " + self.scenario

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
        print "Rank %d: creat sim map"%rank

    #@batch_handler.log_timing
    def make_delta_sim(self):
        r"""this produces self.sim_map_delta"""
        print "making sim in units of overdensity"
        freq_axis = self.sim_map.get_axis('freq') / 1.e6
        z_axis = units.nu21 / freq_axis - 1.0

        simobj = corr21cm.Corr21cm()
        T_b = simobj.T_b(z_axis) * 1e-3

        self.sim_map_delta = copy.deepcopy(self.sim_map)
        self.sim_map_delta /= T_b[:, np.newaxis, np.newaxis]

    #@batch_handler.log_timing
    def make_opt_sim(self):
        r"""this produces self.sim_map_optsim"""

        print "making sim of optically-selected galaxies"
        selection_function = \
                self.datapath_db.fetch_multi(self.params['selection_file'])

        poisson_vect = np.vectorize(np.random.poisson)

        print self.sim_map_delta

        mean_num_gal = (self.sim_map_delta + 1.) * selection_function
        print mean_num_gal

        self.sim_map_optsim = poisson_vect(mean_num_gal)
        self.sim_map_optsim = \
            algebra.make_vect(self.sim_map_optsim.astype(float), axis_names=('freq', 'ra', 'dec'))

        if self.params['optcatalog_file']:
            optical_catalog = \
                self.datapath_db.fetch_multi(self.params['optcatalog_file'])

            # convert from delta to N
            optical_catalog = (optical_catalog + 1.) * selection_function

            print np.sum(optical_catalog), np.sum(self.sim_map_optsim)

        self.sim_map_optsim = self.sim_map_optsim / selection_function - 1.

        self.sim_map_optsim.copy_axis_info(self.sim_map_delta)

    #@batch_handler.log_timing
    def convolve_by_beam(self):
        r"""this produces self.sim_map_withbeam"""
        print "convolving simulation by beam"
        beamobj = beam.GaussianBeam(self.beam_data, self.freq_data)
        self.sim_map_withbeam = beamobj.apply(self.sim_map)

    #@batch_handler.log_timing
    def degrade_to_common_res(self):
        r"""this produces self.sim_map_degrade"""
        print "degrading to common resolution"
        # this depends on having simulations with the means subtracted
        if self.sim_map_meansub is None:
            self.subtract_mean()

        beam_diff = np.sqrt(max(1.1 * self.beam_data) ** 2 - (self.beam_data) ** 2)
        common_resolution = beam.GaussianBeam(beam_diff, self.freq_data)
        # Convolve to a common resolution.
        self.sim_map_degrade = common_resolution.apply(self.sim_map_meansub)
        #sim_map[np.isinf(sim_map)] = 0.
        #sim_map[np.isnan(sim_map)] = 0.

    #@batch_handler.log_timing
    def subtract_mean(self):
        r"""this produces self.sim_map_meansub"""
        print "subtracting mean from simulation"
        # this depends on having simulations convolved by the beam
        if self.sim_map_withbeam is None:
            self.convolve_by_beam()

        self.sim_map_meansub = copy.deepcopy(self.sim_map_withbeam)
        print "sim meansub using: " + self.params['weightfile']
        #noise_inv = self.datapath_db.fetch_multi(self.params['weightfile'])
        noise_inv = algebra.make_vect(algebra.load(self.params['weightfile']))
        means = np.sum(np.sum(noise_inv * self.sim_map_meansub, -1), -1)
        means /= np.sum(np.sum(noise_inv, -1), -1)
        means.shape += (1, 1)
        self.sim_map_meansub -= means
        # the weights will be zero in some places
        self.sim_map_meansub[noise_inv < 1.e-20] = 0.
        #self.sim_map_meansub[np.isinf(self.sim_map)] = 0.
        #self.sim_map_meansub[np.isnan(self.sim_map)] = 0.
