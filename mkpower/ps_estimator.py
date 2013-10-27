#! /usr/bin/env python

import scipy as sp
import numpy as np
from numpy.fft import *
import scipy.linalg as linalg
import multiprocessing as mp

from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from scipy import integrate
from math import *
from sys import *
import matplotlib.pyplot as plt
from mpi4py import MPI
import os

import functions


pi = np.pi
deg2rad = pi/180.

params_init = {
    'processes' : 1,

    'sim_root'  : './',
    'ssm_root'  : './',
    'gbt_root'  : './',
    'opt_root'  : './',

    'sim_numb'  : 100,
    'sim_fact'  : 0,

    'cut_list'  : [],

    'boxshape'  : (128,128,128),
    'discrete'  : 3,

    'kbin_num'  : 20,
    'kbin_min'  : -1.,
    'kbin_max'  : -1.,

    'ps_root'   : './',
    'ps_name'   : 'ps',

    # 'None': return error, 'auto': , 'cros': , 'wigglez', '2df'
    'ps_type'   : 'None', 
    'ps_mode'   : None,

    'est_transfer' : False,
    'est_powerspc' : False,
    'est_gausserr' : False,
    'est_noiseerr' : False,
    'est_powersim' : False,
    'est_powershn' : False,

    'do_2dpower' : True,
    'do_1dpower' : True,
}
prefix = 'pse_'

class PowerSpectrumEstimator(object):
    """Calculate The Power Spectrum"""
    def __init__(self, parameter_file_or_dict=None, feedback=1):
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, 
                                      params_init, 
                                      prefix=prefix, 
                                      feedback=feedback)

        self.feedback=feedback

    def mpiexecute(self, nprocesses=1):
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        self.execute(nprocesses=1, comm=comm, rank=rank, size=size)

    def execute(self, nprocesses=1, comm=None, rank=0, size=1):

        if self.params['est_transfer']:
            self.estimate(comm, rank, size, 'rf')
            self.estimate(comm, rank, size, 'tr')
        if self.params['est_powerspc']:
            self.estimate(comm, rank, size, 'ps')
        if self.params['est_powershn']:
            self.estimate(comm, rank, size, 'sn')
        if self.params['est_gausserr']:
            self.estimate(comm, rank, size, 'ns')
        if self.params['est_noiseerr']:
            self.estimate(comm, rank, size, 'ne')
        if self.params['est_powersim']:
            self.estimate(comm, rank, size, 'si')

    def estimate(self, comm, rank, size, step):
        params = self.params

        if rank == 0:
            imap_list, nmap_list, tabs_list = self.prepare_maps(params, step)
        else:
            imap_list = []
            nmap_list = []
            tabs_list = []

        if comm != None:
            print 'broadcast map list'
            comm.barrier()
            imap_list = comm.bcast(imap_list, root=0)
            nmap_list = comm.bcast(nmap_list, root=0)
            tabs_list = comm.bcast(tabs_list, root=0)
            comm.barrier()

        pair_numb = len(imap_list)
        if rank < pair_numb:
            imap_list = imap_list[rank::size]
            nmap_list = nmap_list[rank::size]
            tabs_list = tabs_list[rank::size]
            print "rank %03d have %d maps to process"%(rank, len(imap_list))

            if step == 'ne':
                factors = ['square', 'square']
            elif params['sim_fact'] != 0:
                if step == 'rf':   factors = [params['sim_fact'], params['sim_fact']]
                elif step == 'tr': factors = [1., params['sim_fact']]
                else:  factors = None
            else: factors = None
            ps_2d, kn_2d, ps_1d, kn_1d =\
                self.process_maps(params, imap_list, nmap_list, 
                                  tabs_list, step, factors)

        if comm != None:
            if rank < pair_numb:
                if rank != 0:
                    print "rank %d: sending data"%rank
                    comm.ssend(ps_2d, dest=0, tag=11)
                    comm.ssend(kn_2d, dest=0, tag=12)
                    comm.ssend(ps_1d, dest=0, tag=13)
                    comm.ssend(kn_1d, dest=0, tag=14)
                    print "rank %d: data sent"%rank
                else:
                    if size > pair_numb:
                        active_rank_num = pair_numb
                    else:
                        active_rank_num = size
                    for i in range(1, active_rank_num):
                        print "rank %d: receiving data from rank %d"%(rank, i)
                        ps_2d = np.concatenate([ps_2d, comm.recv(source=i, 
                                                                 tag=11)],
                                                                 axis=0)
                        kn_2d = np.concatenate([ps_2d, comm.recv(source=i, 
                                                                 tag=12)], 
                                                                 axis=0)
                        ps_1d = np.concatenate([ps_1d, comm.recv(source=i, 
                                                                 tag=13)],
                                                                 axis=0)
                        kn_1d = np.concatenate([ps_1d, comm.recv(source=i, 
                                                                 tag=14)], 
                                                                 axis=0)
                        print "rank %d: received data from rank %d"%(rank, i)
                    for i in range(active_rank_num, size):
                        message = comm.recv(source=i, tag=15)
            else:
                comm.ssend(1, dest=0, tag=15)

            comm.barrier()

        if rank == 0:
            self.save(params, ps_2d, kn_2d, ps_1d, kn_1d, step)

        if comm != None:
            comm.barrier()

    def get_kbin_edges(self, params):
        k_e = np.logspace(np.log10(params['kbin_min']), 
                          np.log10(params['kbin_max']), 
                          num=params['kbin_num'])
        k_s = k_e[-1]/k_e[-2]
        k_e = np.append( k_e, k_e[-1]*k_s ) / (np.sqrt(k_s))
        return k_e

    def apply_cut_list(self, nmap):
        if len(self.params) != 0:
            nmap[self.params['cut_list']] = 0.
        return nmap

    def ps_estimate(self, params, imap_pair, nmap_pair, factors=None):
        
        print "Estimate power using map from"
        print len(imap_pair[0])*'-' +\
              "\n [%s\n *%s]\n x[%s\n *%s]\n"%\
              ( nmap_pair[0], imap_pair[0], nmap_pair[1], imap_pair[1]) +\
              len(imap_pair[0])*'-' + '\n'

        imap1 = algebra.make_vect(algebra.load(imap_pair[0]))
        imap2 = algebra.make_vect(algebra.load(imap_pair[1]))

        if factors != None:
            if factors[0] == 'square':
                imap1 = imap1**2
            else:
                imap1 *= factors[0]
            if factors[1] == 'square':
                imap2 = imap2**2
            else:
                imap2 *= factors[1]
        if nmap_pair[0] != None:
            nmap1 = algebra.make_vect(algebra.load(nmap_pair[0]))
        else:
            print "No weighting given, use 1"
            nmap1 = np.ones_like(imap1)
        if nmap_pair[1] != None:
            nmap2 = algebra.make_vect(algebra.load(nmap_pair[1]))
        else:
            print "No weighting given, use 1"
            nmap2 = np.ones_like(imap2)
        
        self.apply_cut_list(nmap1)
        self.apply_cut_list(nmap2)

        k_edges_p = self.get_kbin_edges(params)
        k_edges_v = self.get_kbin_edges(params)

        ps_box = functions.BOX(imap1, imap2, nmap1, nmap2)
        ps_box.mapping_to_xyz()
        ps_box.estimate_ps_3d()
        ps_box.convert_ps_to_unitless()
        ps_box.convert_3dps_to_2dps(k_edges_p, k_edges_v)
        ps_box.convert_3dps_to_1dps(k_edges_p)

        return ps_box.ps_2d, ps_box.kn_2d, ps_box.ps_1d, ps_box.kn_1d

    def process_maps(self, params, imap_list, nmap_list, tabs_list, step, 
                     factors=None):

        ps_2d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        kn_2d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        ps_1d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        kn_1d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        for i in range(len(imap_list)):
            ps_2d[i], kn_2d[i], ps_1d[i], kn_1d[i] =\
                self.ps_estimate(params, imap_list[i], nmap_list[i], factors)
            if len(tabs_list) != 0:
                self.save_section(params, ps_2d[i], kn_2d[i], step, tabs_list[i])

        return ps_2d, kn_2d, ps_1d, kn_1d

    def prepare_maps(self, params, map_set):
        '''
        map_set: 'rf' : prepare maps for reference calculation, 
                 'tr' : prepare maps for transfer function calculation,
                 'ps' : prepare maps for power spectrum calculation,
                 'si' : prepare maps for simulation.
        '''

        ps_mode = params['ps_mode']
        imap_list = []
        nmap_list = []
        tabs_list = []

        if map_set == 'rf' or (map_set == 'si' and params['ps_type'] == 'cros'):
            print 'prepare maps for reference calculation'
            imaps_a = functions.get_mapdict(params['sim_root'], selection='raw')
            #imaps_a = functions.get_mapdict(params['sim_root'], selection='beammeansub')
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            imaps_b = functions.get_mapdict(params['sim_root'], selection='delta')
            #imaps_b = functions.get_mapdict(params['sim_root'], selection='raw')
            #nmaps_b = functions.get_mapdict(params['gbt_root'])
            nmaps_b = functions.get_mapdict(params['opt_root'], selection='selection')

            for i in range(params['sim_numb']):
                #imaps_a = functions.get_mapdict(params['ssm_root']%i)
                imap_list.append([imaps_a[1]['%d'%i], 
                                  #imaps_a[1]['map;0modes'], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  #nmaps_b[1]['weight;%dmodes'%ps_mode]])
                                  nmaps_b])
        elif map_set == 'tr':
            print 'prepare maps for transfer function calculation'
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            imaps_b = functions.get_mapdict(params['sim_root'], selection='delta')
            #imaps_b = functions.get_mapdict(params['sim_root'], selection='raw')
            #imaps_b = functions.get_mapdict(params['sim_root'], selection='beammeansub')
            #nmaps_b = functions.get_mapdict(params['gbt_root'])
            nmaps_b = functions.get_mapdict(params['opt_root'], selection='selection')

            for i in range(params['sim_numb']):
                imaps_a = functions.get_mapdict(params['ssm_root']%i)
                imap_list.append([imaps_a[1]['map;%dmodes'%ps_mode], 
                                  imaps_b[1]['%d'%i]])
                                  #imaps_b[1]['0']])
                                  #imaps_a[1]['map;%dmodes'%ps_mode]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  #nmaps_b[1]['weight;%dmodes'%ps_mode]])
                                  nmaps_b])
        elif map_set == 'ps' and params['ps_type'] == 'cros':
            #imaps_a = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/'\
            #        + 'test_allbeams_27n30_10by7_clean_map_I_1315.npy'
            #nmaps_a = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/'\
            #        + 'test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
            imaps_a = functions.get_mapdict(params['gbt_root'])
            nmaps_a = functions.get_mapdict(params['gbt_root'])
            imap = functions.get_mapdict(params['opt_root'], selection='real')
            nmap = functions.get_mapdict(params['opt_root'], selection='sele')
            #imaps_b = '/home/ycli/data/2df/map/'\
            #        + 'test_allbeams_27n30_10by7_clean_map_I_1315/'\
            #        + 'real_map_2df_delta.npy'
            #nmaps_b = '/home/ycli/data/2df/map/'\
            #        + 'test_allbeams_27n30_10by7_clean_map_I_1315/'\
            #        + 'sele_map_2df.npy'
            imap_list = [[imaps_a[1]['map;%dmodes'%ps_mode],    imap],]
            nmap_list = [[nmaps_a[1]['weight;%dmodes'%ps_mode], nmap],]
        elif map_set == 'ps' and params['ps_type'] == 'auto':
            print 'prepare maps for power spectrum calculation'
            imaps = functions.get_mapdict(params['gbt_root'])
            #section = ['A', 'B', 'C', 'D']
            section = ['A', 'B', 'C']
            for i in range(len(section)):
                for j in range(i+1, len(section)):
                    s1 = '%s_with_%s'%(section[i], section[j])
                    s2 = '%s_with_%s'%(section[j], section[i])
                    imap_list.append([imaps[1]['%s;map;%dmodes'%(s1, ps_mode)],
                                      imaps[1]['%s;map;%dmodes'%(s2, ps_mode)]])
                    nmap_list.append([imaps[1]['%s;noise_inv;%dmodes'%(s1, ps_mode)],
                                      imaps[1]['%s;noise_inv;%dmodes'%(s2, ps_mode)]])
                    tabs_list.append('%s%sx%s%s'%(section[i], section[j], 
                                                  section[j], section[i]))
        elif map_set == 'ns' and params['ps_type'] == 'auto':
            imaps = functions.get_mapdict(params['gbt_root'])
            #section = ['A', 'B', 'C', 'D']
            section = ['A', 'B', 'C']
            for i in range(len(section)):
                for j in range(len(section)):
                    if i == j:
                        continue
                    s = '%s_with_%s'%(section[i], section[j])
                    imap_list.append([imaps[1]['%s;map;%dmodes'%(s, ps_mode)],
                                      imaps[1]['%s;map;%dmodes'%(s, ps_mode)]])
                    nmap_list.append([imaps[1]['%s;noise_inv;%dmodes'%(s, ps_mode)],
                                      imaps[1]['%s;noise_inv;%dmodes'%(s, ps_mode)]])
                    tabs_list.append('%s%s'%(section[i], section[j]))

        elif map_set == 'si' and params['ps_type'] == 'auto':
            print 'prepare maps for simulation calculation'
            imaps_a = functions.get_mapdict(params['sim_root'], selection='raw')
            #imaps_a = functions.get_mapdict(params['sim_root'], selection='beammeansub')
            #imaps_a = functions.get_mapdict(params['sim_root'], selection='degradebeam')
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            imaps_b = functions.get_mapdict(params['sim_root'], selection='raw')
            #imaps_b = functions.get_mapdict(params['sim_root'], selection='beammeansub')
            nmaps_b = functions.get_mapdict(params['gbt_root'])

            for i in range(params['sim_numb']):
                imap_list.append([imaps_a[1]['%d'%i], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  nmaps_b[1]['weight;%dmodes'%ps_mode]])
        elif map_set == 'ne' and params['ps_type'] == 'auto':
            print 'prepare maps for noise error calculation'
            imaps = functions.get_mapdict(params['gbt_root'])
            #section = ['A', 'B', 'C', 'D']
            section = ['A', 'B', 'C']
            for i in range(len(section)):
                for j in range(len(section)):
                    if i == j:
                        continue
                    s = '%s_with_%s'%(section[i], section[j])
                    imap_list.append([imaps[1]['%s;noise_diag;%dmodes'%(s, ps_mode)],
                                      imaps[1]['%s;noise_diag;%dmodes'%(s, ps_mode)]])
                    nmap_list.append([None, None])
                    tabs_list.append('%s%s'%(section[i], section[j]))

        elif map_set == 'ps' and params['ps_type'] == '2df':
            print 'prepare maps of real catalog for 2df power sepctrum'
            imap = functions.get_mapdict(params['opt_root'], selection='real')
            nmap = functions.get_mapdict(params['opt_root'], selection='sele')
            imap_list.append([imap, imap])
            nmap_list.append([nmap, nmap])
        elif map_set == 'sn' and params['ps_type'] == '2df':
            print 'prepare maps of mock catalog for 2df power sepctrum short noise'
            imap = functions.get_mapdict(params['opt_root'], selection='mock')
            nmap = functions.get_mapdict(params['opt_root'], selection='sele')
            #for i in range(len(imap[0])):
            for i in range(100):
                imap_list.append([imap[1]['%d'%i],imap[1]['%d'%i]])
                nmap_list.append([nmap, nmap])
        elif map_set == 'sn' and params['ps_type'] == 'cros':
            print 'prepare maps of mock catalog for 2df power sepctrum short noise'
            #imaps_a = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/'\
            #        + 'test_allbeams_27n30_10by7_clean_map_I_1315.npy'
            #nmaps_a = '/home/ycli/data/parkes/maps_by_chris/pks_27n30_10by7/'\
            #        + 'test_allbeams_27n30_10by7_noise_weight_I_1315.npy'
            imaps_a = functions.get_mapdict(params['gbt_root'])
            nmaps_a = functions.get_mapdict(params['gbt_root'])
            imap = functions.get_mapdict(params['opt_root'], selection='mock')
            nmap = functions.get_mapdict(params['opt_root'], selection='sele')
            for i in range(len(imap[0])):
                imap_list.append([imaps_a[1]['map;%dmodes'%ps_mode], 
                                  imap[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  nmap])
        elif map_set == 'ps' and params['ps_type'] == 'wigglez':
            print 'prepare maps of mock catalog for wigglez power sepctrum short noise'
            imap = functions.get_mapdict(params['opt_root'], selection='data')
            nmap = functions.get_mapdict(params['opt_root'].replace('binned_delta', 'binned'), 
                                         selection='separable')
            imap_list.append([imap, imap])
            nmap_list.append([nmap, nmap])
        elif map_set == 'sn' and params['ps_type'] == 'wigglez':
            print 'prepare maps of mock catalog for wigglez power sepctrum short noise'
            imap = functions.get_mapdict(params['opt_root'])
            nmap = functions.get_mapdict(params['opt_root'].replace('binned_delta', 'binned'), 
                                         selection='separable')
            for i in range(len(imap[0])):
                imap_list.append([imap[1]['%d'%i],imap[1]['%d'%i]])
                nmap_list.append([nmap, nmap])
            
        elif map_set == 'si' and params['ps_type'] == '2df':
            print 'prepare maps of simulation maps for 2df power sepctrum'
            imaps_a = functions.get_mapdict(params['sim_root'], selection='delta')
            nmaps_a = functions.get_mapdict(params['opt_root'], selection='selection')

            imaps_b = functions.get_mapdict(params['sim_root'], selection='delta')
            nmaps_b = functions.get_mapdict(params['opt_root'], selection='selection')

            for i in range(params['sim_numb']):
                imap_list.append([imaps_a[1]['%d'%i], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a, nmaps_b])
        else:
            print 'error: map_set error!'
            exit()
        return imap_list, nmap_list, tabs_list

    def save_section(self, params, ps_2d, kn_2d, step, tab):
        if step == 'ne' and params['ps_type'] == 'auto':
            ps_2d = np.sqrt(ps_2d)

        k_bin = np.logspace(np.log10(params['kbin_min']), 
                            np.log10(params['kbin_max']), 
                            num=params['kbin_num'])

        k_axes_2d = ("k_p", "k_v")
        info_2d = {'axes': k_axes_2d, 'type': 'vect'}
        info_2d['k_p_delta']  = k_bin[1]/k_bin[0]
        info_2d['k_p_centre'] = k_bin[params['kbin_num']//2]
        info_2d['k_v_delta']  = k_bin[1]/k_bin[0]
        info_2d['k_v_centre'] = k_bin[params['kbin_num']//2]

        ps_2d  = algebra.make_vect(ps_2d, axis_names=k_axes_2d)
        kn_2d  = algebra.make_vect(kn_2d, axis_names=k_axes_2d)
        ps_2d.info = info_2d
        kn_2d.info = info_2d

        if params['ps_mode'] == None:
            file_name = '%s_%s_'%(params['ps_type'], step)
        else:
            file_name = '%s_%s_%dmode_'%(params['ps_type'], step, params['ps_mode'])

        file_root = params['ps_root'] + params['ps_name'] + '/'
        if not os.path.exists(file_root):
            os.makedirs(file_root)

        algebra.save(file_root + file_name + '2dpow_%s'%tab, ps_2d)
        algebra.save(file_root + file_name + '2dkmn_%s'%tab, kn_2d)

    def save(self, params, ps_2d, kn_2d, ps_1d, kn_1d, step):

        print ps_2d.shape
        print ps_2d.flatten().max(), ps_2d.flatten().min()
        ps_2d_mean = np.mean(ps_2d, axis=0)
        ps_2d_std  = np.std(ps_2d, axis=0)
        ps_1d_mean = np.mean(ps_1d, axis=0)
        ps_1d_std  = np.std(ps_1d, axis=0)

        kn_2d_mean = np.mean(kn_2d, axis=0)
        kn_2d_std  = np.std(ps_2d, axis=0)
        kn_1d_mean = np.mean(kn_1d, axis=0)
        kn_1d_std  = np.std(ps_1d, axis=0)

        if step == 'ps' and params['ps_type'] == 'auto':
            print ps_2d.shape[0]
            ps_2d_std /= sqrt(ps_2d.shape[0])
            ps_1d_std /= sqrt(ps_1d.shape[0])
        if step == 'ne' and params['ps_type'] == 'auto':
            ps_2d_mean = np.sqrt(ps_2d_mean)
            ps_2d_std = np.sqrt(ps_2d_std)
            ps_2d = np.sqrt(ps_2d)

        k_bin = np.logspace(np.log10(params['kbin_min']), 
                            np.log10(params['kbin_max']), 
                            num=params['kbin_num'])

        k_axes_2d = ("k_p", "k_v")
        info_2d = {'axes': k_axes_2d, 'type': 'vect'}
        info_2d['k_p_delta']  = k_bin[1]/k_bin[0]
        info_2d['k_p_centre'] = k_bin[params['kbin_num']//2]
        info_2d['k_v_delta']  = k_bin[1]/k_bin[0]
        info_2d['k_v_centre'] = k_bin[params['kbin_num']//2]

        ps_2d_mean = algebra.make_vect(ps_2d_mean, axis_names=k_axes_2d)
        kn_2d_mean = algebra.make_vect(kn_2d_mean, axis_names=k_axes_2d)
        ps_2d_std  = algebra.make_vect(ps_2d_std, axis_names=k_axes_2d)
        kn_2d_std  = algebra.make_vect(kn_2d_std, axis_names=k_axes_2d)

        ps_2d_mean.info = info_2d
        kn_2d_mean.info = info_2d
        ps_2d_std.info = info_2d
        kn_2d_std.info = info_2d

        k_axes_1d = ("k",)
        info_1d = {'axes': k_axes_1d, 'type': 'vect'}
        info_1d['k_delta']  = k_bin[1]/k_bin[0]
        info_1d['k_centre'] = k_bin[params['kbin_num']//2]

        ps_1d_mean = algebra.make_vect(ps_1d_mean, axis_names=k_axes_1d)
        kn_1d_mean = algebra.make_vect(kn_1d_mean, axis_names=k_axes_1d)
        ps_1d_std  = algebra.make_vect(ps_1d_std, axis_names=k_axes_1d)
        kn_1d_std  = algebra.make_vect(kn_1d_std, axis_names=k_axes_1d)

        ps_1d_mean.info = info_1d
        kn_1d_mean.info = info_1d
        ps_1d_std.info = info_1d
        kn_1d_std.info = info_1d

        if params['ps_mode'] == None:
            file_name = '%s_%s_'%(params['ps_type'], step)
        else:
            file_name = '%s_%s_%dmode_'%(params['ps_type'], step, params['ps_mode'])

        file_root = params['ps_root'] + params['ps_name'] + '/'
        if not os.path.exists(file_root):
            os.makedirs(file_root)

        print file_root

        algebra.save(file_root + file_name + '2dpow', ps_2d_mean)
        algebra.save(file_root + file_name + '2derr', ps_2d_std )
        algebra.save(file_root + file_name + '2dkmn', kn_2d_mean)
        algebra.save(file_root + file_name + '1dpow', ps_1d_mean)
        algebra.save(file_root + file_name + '1derr', ps_1d_std )
        algebra.save(file_root + file_name + '1dkmn', kn_1d_mean)

        #if step == 'ps':
        k_axes_each_2d = ("map", "k_p", "k_v")
        info_each_2d = {'axes': k_axes_each_2d, 'type': 'vect'}
        info_each_2d['map_delta']  = 1.
        info_each_2d['map_centre'] = range(ps_2d.shape[0])[ps_2d.shape[0]//2]
        info_each_2d['k_p_delta']  = k_bin[1]/k_bin[0]
        info_each_2d['k_p_centre'] = k_bin[params['kbin_num']//2]
        info_each_2d['k_v_delta']  = k_bin[1]/k_bin[0]
        info_each_2d['k_v_centre'] = k_bin[params['kbin_num']//2]

        k_axes_each_1d = ("map", "k")
        info_each_1d = {'axes': k_axes_each_1d, 'type': 'vect'}
        info_each_1d['map_delta']  = 1.
        info_each_1d['map_centre'] = range(ps_1d.shape[0])[ps_1d.shape[0]//2]
        info_each_1d['k_delta']  = k_bin[1]/k_bin[0]
        info_each_1d['k_centre'] = k_bin[params['kbin_num']//2]

        ps_2d_each = algebra.make_vect(ps_2d, axis_names=k_axes_each_2d)
        ps_2d_each.info = info_each_2d
        algebra.save(file_root + file_name + '2draw', ps_2d_each)

        ps_1d_each = algebra.make_vect(ps_1d, axis_names=k_axes_each_1d)
        ps_1d_each.info = info_each_1d
        algebra.save(file_root + file_name + '1draw', ps_1d_each)


if __name__ == '__main__':
    import sys
