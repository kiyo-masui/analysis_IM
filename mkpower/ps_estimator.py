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
import mkpower
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

    # 'None': return error, 'auto': , 'cros': , 'wigg'
    'ps_type'   : 'None', 
    'ps_mode'   : 10,

    'est_transfer' : True,
    'est_powerspc' : True,
    'est_gausserr' : True,
    'est_powersim' : True
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
        if self.params['est_gausserr']:
            self.estimate(comm, rank, size, 'ns')
        if self.params['est_powersim']:
            self.estimate(comm, rank, size, 'si')

    def estimate(self, comm, rank, size, step):
        params = self.params

        if rank == 0:
            imap_list, nmap_list = self.prepare_maps(params, step)
        else:
            imap_list = []
            nmap_list = []

        if comm != None:
            print 'broadcast map list'
            comm.barrier()
            imap_list = comm.bcast(imap_list, root=0)
            nmap_list = comm.bcast(nmap_list, root=0)
            comm.barrier()

        pair_numb = len(imap_list)
        if rank < pair_numb:
            imap_list = imap_list[rank::size]
            nmap_list = nmap_list[rank::size]
            print "rank %03d have %d maps to process"%(rank, len(imap_list))

            if params['sim_fact'] != 0:
                if step == 'rf':   factors = [params['sim_fact'], params['sim_fact']]
                elif step == 'tr': factors = [1., params['sim_fact']]
                else:  factors = None
            else: factors = None
            ps, kn = self.process_maps(params, imap_list, nmap_list, factors)

        if comm != None:
            if rank < pair_numb:
                if rank != 0:
                    print "rank %d: sending data"%rank
                    comm.ssend(ps, dest=0, tag=11)
                    comm.ssend(kn, dest=0, tag=12)
                    print "rank %d: data sent"%rank
                else:
                    if size > pair_numb:
                        active_rank_num = pair_numb
                    else:
                        active_rank_num = size
                    for i in range(1, active_rank_num):
                        print "rank %d: receiving data from rank %d"%(rank, i)
                        ps = np.concatenate([ps, comm.recv(source=i, tag=11)], axis=0)
                        kn = np.concatenate([ps, comm.recv(source=i, tag=12)], axis=0)
                        print "rank %d: received data from rank %d"%(rank, i)
                    for i in range(active_rank_num, size):
                        message = comm.recv(source=i, tag=11)
            else:
                comm.ssend(1, dest=0, tag=11)

            comm.barrier()

        if rank == 0:
            self.save(params, ps, kn, step)

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
        
        imap1 = algebra.make_vect(algebra.load(imap_pair[0]))
        imap2 = algebra.make_vect(algebra.load(imap_pair[1]))
        if factors != None:
            imap1 *= factors[0]
            imap2 *= factors[1]
        if nmap_pair[0] != None:
            nmap1 = algebra.make_vect(algebra.load(nmap_pair[0]))
        else:
            nmap1 = np.ones_like(imap1)
        if nmap_pair[1] != None:
            nmap2 = algebra.make_vect(algebra.load(nmap_pair[1]))
        else:
            nmap2 = np.ones_like(imap2)
        
        self.apply_cut_list(nmap1)
        self.apply_cut_list(nmap2)

        k_edges_p = self.get_kbin_edges(params)
        k_edges_v = self.get_kbin_edges(params)

        ps_box = functions.BOX(params['boxshape'], imap1, imap2, nmap1, nmap2)
        ps_box.mapping_to_xyz()
        ps_box.estimate_ps_3d()
        ps_box.convert_ps_to_unitless()
        ps_box.convert_3dps_to_2dps(k_edges_p, k_edges_v)

        return ps_box.ps_2d, ps_box.kn_2d

    def process_maps(self, params, imap_list, nmap_list, factors=None):

        ps = np.zeros(shape=(len(imap_list), params['kbin_num'], params['kbin_num']))
        kn = np.zeros(shape=(len(imap_list), params['kbin_num'], params['kbin_num']))
        for i in range(len(imap_list)):
            ps[i], kn[i] = self.ps_estimate(params, imap_list[i], nmap_list[i], factors)

        return ps, kn

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

        if map_set == 'rf' or (map_set == 'si' and params['ps_type'] == 'cros'):
            print 'prepare maps for reference calculation'
            imaps_a = functions.get_mapdict(params['sim_root'], selection='raw')
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            #imaps_b = functions.get_mapdict(params['sim_root'], selection='delta')
            imaps_b = functions.get_mapdict(params['sim_root'], selection='raw')
            nmaps_b = functions.get_mapdict(params['opt_root'], selection='selection')

            for i in range(params['sim_numb']):
                imap_list.append([imaps_a[1]['%d'%i], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  nmaps_b])
        elif map_set == 'tr':
            print 'prepare maps for transfer function calculation'
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            #imaps_b = functions.get_mapdict(params['sim_root'], selection='delta')
            imaps_b = functions.get_mapdict(params['sim_root'], selection='raw')
            nmaps_b = functions.get_mapdict(params['opt_root'], selection='selection')

            for i in range(params['sim_numb']):
                imaps_a = functions.get_mapdict(params['ssm_root']%i)
                imap_list.append([imaps_a[1]['map;%dmodes'%ps_mode], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  nmaps_b])
        elif map_set == 'ps' and params['ps_type'] == 'auto':
            print 'prepare maps for power spectrum calculation'
            imaps = functions.get_mapdict(params['gbt_root'])
            section = ['A', 'B', 'C', 'D']
            for i in range(len(section)):
                for j in range(i+1, len(section)):
                    s1 = '%s_with_%s'%(section[i], section[j])
                    s2 = '%s_with_%s'%(section[j], section[i])
                    imap_list.append([imaps[1]['%s;map;%dmodes'%(s1, ps_mode)],
                                      imaps[1]['%s;map;%dmodes'%(s2, ps_mode)]])
                    nmap_list.append([imaps[1]['%s;noise_inv;%dmodes'%(s1, ps_mode)],
                                      imaps[1]['%s;noise_inv;%dmodes'%(s2, ps_mode)]])
        elif map_set == 'ns' and params['ps_type'] == 'auto':
            imaps = functions.get_mapdict(params['gbt_root'])
            section = ['A', 'B', 'C', 'D']
            for i in range(len(section)):
                for j in range(len(section)):
                    if i == j:
                        continue
                    s = '%s_with_%s'%(section[i], section[j])
                    imap_list.append([imaps[1]['%s;map;%dmodes'%(s, ps_mode)],
                                      imaps[1]['%s;map;%dmodes'%(s, ps_mode)]])
                    nmap_list.append([imaps[1]['%s;noise_inv;%dmodes'%(s, ps_mode)],
                                      imaps[1]['%s;noise_inv;%dmodes'%(s, ps_mode)]])
        elif map_set == 'si' and params['ps_type'] == 'auto':
            print 'prepare maps for simulation calculation'
            imaps_a = functions.get_mapdict(params['sim_root'], selection='raw')
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            imaps_b = functions.get_mapdict(params['sim_root'], selection='raw')
            nmaps_b = functions.get_mapdict(params['gbt_root'])

            for i in range(params['sim_numb']):
                imap_list.append([imaps_a[1]['%d'%i], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  nmaps_b[1]['weight;%dmodes'%ps_mode]])
        else:
            print 'error: map_set error!'
            exit()
        return imap_list, nmap_list

    def save(self, params, ps, kn, step):

        ps_mean = np.mean(ps, axis=0)
        ps_std  = np.std(ps, axis=0)

        kn_mean = np.mean(kn, axis=0)
        kn_std  = np.std(ps, axis=0)

        if step == 'ps' and params['ps_type'] == 'auto':
            print ps.shape[0]
            ps_std /= sqrt(ps.shape[0])

        k_bin = np.logspace(np.log10(params['kbin_min']), 
                            np.log10(params['kbin_max']), 
                            num=params['kbin_num'])

        k_axes = ("k_p", "k_v")
        info = {'axes': k_axes, 'type': 'vect'}
        info['k_p_delta']  = k_bin[1]/k_bin[0]
        info['k_p_centre'] = k_bin[params['kbin_num']//2]
        info['k_v_delta']  = k_bin[1]/k_bin[0]
        info['k_v_centre'] = k_bin[params['kbin_num']//2]

        #k_edges_p = self.get_kbin_edges(params)
        #k_edges_v = self.get_kbin_edges(params)
        #info['k_p_edges'] = k_edges_p
        #info['k_v_edges'] = k_edges_v

        ps_mean = algebra.make_vect(ps_mean, axis_names=k_axes)
        kn_mean = algebra.make_vect(kn_mean, axis_names=k_axes)
        ps_std  = algebra.make_vect(ps_std, axis_names=k_axes)
        kn_std  = algebra.make_vect(kn_std, axis_names=k_axes)

        ps_mean.info = info
        kn_mean.info = info
        ps_std.info = info
        kn_std.info = info

        file_name = '%s_%s_%dmode_'%(params['ps_type'], step,params['ps_mode'])

        file_root = params['ps_root'] + params['ps_name'] + '/'
        if not os.path.exists(file_root):
            os.makedirs(file_root)

        print file_root

        algebra.save(file_root + file_name + '2dpow', ps_mean)
        algebra.save(file_root + file_name + '2derr', ps_std )
        algebra.save(file_root + file_name + '2dkmn', kn_mean)

if __name__ == '__main__':
    import sys
    if len(sys.argv)==2 :
        PowerSpectrumMaker(str(sys.argv[1])).execute()
    elif len(sys.argv)>2 :
        print 'Maximun one argument, a parameter file name.'
    else :
        PowerSpectrumMaker().execute()
