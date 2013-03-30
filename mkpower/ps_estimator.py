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

    'cut_list'  : [],

    'boxshape'  : (128,128,128),
    'discrete'  : 3,

    'kbin_num'  : 20,
    'kbin_min'  : -1,
    'kbin_max'  : -1,

    'ps_root'   : './',
    'ps_name'   : 'ps',

    # 'None': return error, 'auto': , 'cros': , 'wigg'
    'ps_type'   : 'None', 
    'ps_mode'   : 10,
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

        execute(self, nprocesses=1, comm=None, rank=None, size=None):

    def execute(self, nprocesses=1, comm=None, rank=0, size=1):

        self.estimate(comm, rank, size, 'rf')
        self.estimate(comm, rank, size, 'tr')
        self.estimate(comm, rank, size, 'ps')
        self.estimate(comm, rank, size, 'ns')

    def estimate(self, comm, rank, size, step):
        params = self.params

        if rank == 0:
            imap_list, nmap_list = self.prepare_maps(params, step)
        if comm != None:
            comm.bcast(imap_list, root=0)
            comm.bcast(nmap_list, root=0)
            comm.barrier()

        pair_numb = len(imap_list)
        if rank < pair_numb:
            imap_list = imap_list[rank::size]
            nmap_list = nmap_list[rank::size]

            ps, kn = self.process_maps(params, imap_list, nmap_list)

        if comm != None:
            if rank != 0 and rank < pair_numb:
                print "rank %d: sending data"%rank
                comm.ssend(ps, dest=0, tag=11)
                comm.ssend(kn, dest=0, tag=11)
                print "rank %d: data sent"%rank
            elif rank == 0:
                if size > pair_numb:
                    active_rank_num = pair_numb
                else:
                    active_rank_num = size
                for i in range(1, active_rank_num):
                    print "rank %d: receiving data from rank %d"%(rank, i)
                    ps = np.append(ps, comm.recv(source=i, tag=11), axis=0)
                    kn = np.append(kn, comm.recv(source=i, tag=11), axis=0)
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
        k_e = np.logspace(params['kbin_min'], 
                          params['kbin_max'], 
                          num=params['kbin_num'])
        k_s = k_e[-1]/k_e[-2]
        k_e = np.append( k_e, k_e[-1]*k_s ) / (0.5 * k_s)
        return k_e

    def apply_cut_list(self, nmap):
        if len(self.params) != 0:
            nmap[self.params['cut_list']] = 0.
        return nmap

    def ps_estimate(self, params, imap_pair, nmap_pair):
        
        imap1 = algebra.make_vect(algebra.load(imap_pair[0]))
        imap2 = algebra.make_vect(algebra.load(imap_pair[1]))
        nmap1 = algebra.make_vect(algebra.load(nmap_pair[0]))
        nmap2 = algebra.make_vect(algebra.load(nmap_pair[1]))
        
        self.apply_cut_list(nmap1)
        self.apply_cut_list(nmap2)

        k_edges_p = self.get_kbin_edges()
        k_edges_v = self.get_kbin_edges()

        ps_box = functions.BOX(params['boxshape'], imap1, imap2, nmap1, nmap2)
        ps_box.mapping_to_xyz()
        ps_box.estimat_ps_3d()
        ps_box.convert_ps_to_unitless()
        ps_box.convert_3dps_to_2dps(k_edges_p, k_edges_v)

        return ps_box.ps_2d, ps_box.kn_2d

    def process_maps(self, params, imap_list, nmap_list):

        ps = np.zeros(shape=(len(imap_list), params['kbin_num'], params['kbin_num']))
        kn = np.zeros(shape=(len(imap_list), params['kbin_num'], params['kbin_num']))
        for i in range(len(imap_list)):
            ps[i], kn[i] = self.ps_estimate(params, imap_list[i], nmap_list[i])

        return ps, kn

    def prepare_maps(self, params, map_set):
        '''
        map_set: 'rf' : prepare maps for reference calculation, 
                 'tr' : prepare maps for transfer function calculation,
                 'ps' : prepare maps for power spectrum calculation,
        '''

        ps_mode = params['ps_mode']
        imap_list = []
        nmap_list = []

        if map_set == 'rf':
            imaps_a = functions.get_mapdict(params['sim_root'], selection='raw')
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            imaps_b = functions.get_mapdict(params['sim_root'], selection='delta')
            nmaps_b = functions.get_mapdict(params['opt_root'], selection='selection')

            for i in range(params['sim_numb']):
                imap_list.append([imaps_a[1]['%d'%i], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  nmaps_b])
        elif map_set == 'tr':
            nmaps_a = functions.get_mapdict(params['gbt_root'])

            imaps_b = functions.get_mapdict(params['sim_root'], selection='delta')
            nmaps_b = functions.get_mapdict(params['opt_root'], selection='selection')

            for i in range(params['sim_numb']):
                imaps_a = functions.get_mapdict(params['ssm_root']%i)
                imap_list.append([imaps_a[1]['map;%dmodes'%ps_mode], 
                                  imaps_b[1]['%d'%i]])
                nmap_list.append([nmaps_a[1]['weight;%dmodes'%ps_mode], 
                                  nmaps_b])
         elif map_set == 'ps' and params['ps_type'] == 'auto':
            imaps = function.get_mapdict[params['gbt_root']]
            section = ['A', 'B', 'C', 'D']
            for i in range(len(section)):
                for j in range(i+1, len(section)):
                    s1 = '%s_with_%s'%(section[i], section[j])
                    s2 = '%s_with_%s'%(section[j], section[i])
                    imaps_list.append([imaps[1]['%s;map;%dmodes'%(s1, ps_mode)],
                                       imaps[1]['%s;map;%dmodes'%(s2, ps_mode)]])
                    nmaps_list.append([imaps[1]['%s;weight;%dmodes'%(s1, ps_mode)],
                                       imaps[1]['%s;weight;%dmodes'%(s2, ps_mode)]])
         elif map_set == 'ns' and params['ps_type'] == 'auto':
            imaps = function.get_mapdict[params['gbt_root']]
            section = ['A', 'B', 'C', 'D']
            for i in range(len(section)):
                for j in range(len(section)):
                    s = '%s_with_%s'%(section[i], section[j])
                    imaps_list.append([imaps[1]['%s;map;%dmodes'%(s, ps_mode)],
                                       imaps[1]['%s;map;%dmodes'%(s, ps_mode)]])
                    nmaps_list.append([imaps[1]['%s;weight;%dmodes'%(s, ps_mode)],
                                       imaps[1]['%s;weight;%dmodes'%(s, ps_mode)]])
         else:
            print 'error: map_set error!'
            exit()
         return imap_list, nmap_list

    def save(self, params, ps, kn, step):

        ps_mean = np.mean(ps, axis=0)
        ps_std  = np.std(ps, axis=0)

        kn_mean = np.mean(kn, axis=0)
        kn_std  = np.std(ps, axis=0)

        if setp == 'power' and params['ps_type'] == 'auto':
            ps_std /= sqrt(ps.shape[0])

        k_bin = np.logspace(params['kbin_min'], 
                            params['kbin_max'], 
                            num=params['kbin_num'])

        k_axes = ("k_p", "k_v")
        info = {'axes': k_axes, 'type': 'vect'}
        info['k_p_delta']  = k_bin[1]/k_bin[0]
        info['k_p_centre'] = k_bin[params['kbin_num']//2]
        info['k_v_delta']  = k_bin[1]/k_bin[0]
        info['k_v_centre'] = k_bin[params['kbin_num']//2]

        ps_mean = algebra.make_vect(ps_mean, axis_names=k_axes)
        kn_mean = algebra.make_vect(kn_mean, axis_names=k_axes)
        ps_std  = algebra.make_vect(ps_std, axis_names=k_axes)
        kn_std  = algebra.make_vect(kn_std, axis_names=k_axes)

        ps_mean.info = info
        kn_mean.info = info
        ps_std.info = info
        kn_std.info = info

        file_name = '%s_%s_%dmode_'%(params['ps_type'], step, params['ps_mode'])

        algebra.save(file_name + '2dpow', ps_mean)
        algebra.save(file_name + '2derr', ps_std )
        algebra.save(file_name + '2dkmn', kn_mean)

if __name__ == '__main__':
    import sys
    if len(sys.argv)==2 :
        PowerSpectrumMaker(str(sys.argv[1])).execute()
    elif len(sys.argv)>2 :
        print 'Maximun one argument, a parameter file name.'
    else :
        PowerSpectrumMaker().execute()
