#! /usr/bin/env python

import scipy as sp
import numpy as np
from numpy.fft import *
import scipy.linalg as linalg
import multiprocessing as mp
import pdb
from core import algebra, hist
from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from scipy import integrate
from math import *
from sys import *
from mpi4py import MPI
import os
import h5py

import functions


pi = np.pi
deg2rad = pi/180.

params_init = {
    'processes' : 1,

    'sim_root'  : './',
    'ssm_root'  : './',
    'gbt_root'  : './',
    'opt_root'  : './',
    'trans_root' : './',
    'transsim_root' : './',

    'sim_numb'  : 50,
    'sim_fact'  : 0.,

    'cut_list'  : [],

    'colour' : None,

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
    'est_noiselev' : False,

    'do_2dpower' : True,
    'do_1dpower' : True,

    # for transfer function 
    'degrade_factor' : 1.4, 
    'telescope' : 'Parkes', 
    'beam_list' : np.arange(13) + 1,
    'modes' : 0,
    'svd_path' : '', 
    'svd_file' : ['',],
    'freq_list' : [],
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
        print 'rank', rank, 'size', size
        self.execute(nprocesses=1, comm=comm, rank=rank, size=size)

    def execute(self, nprocesses=1, comm=None, rank=0, size=1):

        params = self.params

        self.file_root = params['ps_root'] + params['ps_name'] + '/'

        if rank == 0:
            if not os.path.exists(self.file_root):
                os.makedirs(self.file_root)
            self.ps_result = h5py.File(self.file_root + 'ps_result.hd5', 'a')
            print 'ps_result =', self.ps_result
        if self.params['est_powersim']:
            self.estimate(comm, rank, size, 'si')
            holder = params['ps_type']
            self.params['ps_type'] = '2df'
            #self.estimate(comm, rank, size, 'si')
            self.params['ps_type'] = holder
        if self.params['est_transfer'] :
            self.estimate(comm, rank, size, 'tr')
            self.estimate(comm, rank, size, 'rf')
            if params['ps_type'] == 'cros':
                self.estimate(comm, rank, size, 'rb')
        if self.params['est_powerspc']:
            self.estimate(comm, rank, size, 'ps')
            holder = params['ps_type']
            self.params['ps_type'] = '2df'
       	    #self.estimate(comm, rank, size, 'ps')
            self.params['ps_type'] = holder
        if self.params['est_powershn']:
            self.estimate(comm, rank, size, 'sn')
            holder = params['ps_type']
            self.params['ps_type'] = '2df'
            #self.estimate(comm, rank, size, 'sn')
            self.params['ps_type'] = holder
        if self.params['est_gausserr']:
            self.estimate(comm, rank, size, 'ns')
        if self.params['est_noiseerr']:
            self.estimate(comm, rank, size, 'ne')
        if self.params['est_noiselev']:
            self.estimate(comm, rank, size, 'nl')

        if rank == 0:
            self.ps_result.close()

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
            ps_2d, kn_2d, ps_1d, kn_1d,\
                    ps_2d_kpp, kn_2d_kpp, ps_2d_kpn, kn_2d_kpn, \
                    ps_1d_kpp, kn_1d_kpp, ps_1d_kpn, kn_1d_kpn =\
                self.process_maps(params, imap_list, nmap_list, tabs_list, step, factors)
        if comm != None:
            if rank < pair_numb:
                if rank != 0:
                    print "rank %d: sending data"%rank
                    comm.ssend(ps_2d, dest=0, tag=11)
                    comm.ssend(kn_2d, dest=0, tag=12)
                    comm.ssend(ps_1d, dest=0, tag=13)
                    comm.ssend(kn_1d, dest=0, tag=14)
                    comm.ssend(ps_2d_kpp, dest=0, tag=21)
                    comm.ssend(kn_2d_kpp, dest=0, tag=22)
                    comm.ssend(ps_2d_kpn, dest=0, tag=23)
                    comm.ssend(kn_2d_kpn, dest=0, tag=24)
                    comm.ssend(ps_1d_kpp, dest=0, tag=31)
                    comm.ssend(kn_1d_kpp, dest=0, tag=32)
                    comm.ssend(ps_1d_kpn, dest=0, tag=33)
                    comm.ssend(kn_1d_kpn, dest=0, tag=34)
                    comm.ssend(tabs_list, dest=0, tag=35)
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
                        kn_2d = np.concatenate([kn_2d, comm.recv(source=i, 
                                                                 tag=12)], 
                                                                 axis=0)
                        ps_1d = np.concatenate([ps_1d, comm.recv(source=i, 
                                                                 tag=13)],
                                                                 axis=0)
                        kn_1d = np.concatenate([kn_1d, comm.recv(source=i, 
                                                                 tag=14)], 
                                                                 axis=0)
                        ps_2d_kpp = np.concatenate([ps_2d_kpp, comm.recv(source=i, 
                                                                     tag=21)],
                                                                     axis=0)
                        kn_2d_kpp = np.concatenate([kn_2d_kpp, comm.recv(source=i, 
                                                                     tag=22)], 
                                                                     axis=0)
                        ps_2d_kpn = np.concatenate([ps_2d_kpn, comm.recv(source=i, 
                                                                     tag=23)],
                                                                     axis=0)
                        kn_2d_kpn = np.concatenate([kn_2d_kpn, comm.recv(source=i, 
                                                                     tag=24)], 
                                                                     axis=0)
                        ps_1d_kpp = np.concatenate([ps_1d_kpp, comm.recv(source=i, 
                                                                     tag=31)],
                                                                     axis=0)
                        kn_1d_kpp = np.concatenate([kn_1d_kpp, comm.recv(source=i, 
                                                                     tag=32)], 
                                                                     axis=0)
                        ps_1d_kpn = np.concatenate([ps_1d_kpn, comm.recv(source=i, 
                                                                     tag=33)],
                                                                     axis=0)
                        kn_1d_kpn = np.concatenate([kn_1d_kpn, comm.recv(source=i, 
                                                                     tag=34)], 
                                                                     axis=0)
                        print 'up to tabs'
                        tabs_list = np.concatenate([tabs_list, comm.recv(source=i,
                                                                     tag=35)],
                                                                     axis=0)
                        print "rank %d: received data from rank %d"%(rank, i)
                    for i in range(active_rank_num, size):
                        message = comm.recv(source=i, tag=15)
            else:
                comm.ssend(1, dest=0, tag=15)

            comm.barrier()

        if rank == 0:
            if step == 'ps' and len(tabs_list) != 0:
                for i in range(len(ps_2d)):
                    self.save_section(params, ps_2d[i], kn_2d[i], step, tabs_list[i])
            self.save(params, ps_2d, kn_2d, ps_1d, kn_1d, 
                    ps_2d_kpp, kn_2d_kpp, ps_2d_kpn, kn_2d_kpn, 
                    ps_1d_kpp, kn_1d_kpp, ps_1d_kpn, kn_1d_kpn, 
                    step)

        if comm != None:
            comm.barrier()

    def get_kbin_centr(self, params):

        k_e = np.logspace(np.log10(params['kbin_min']), 
                          np.log10(params['kbin_max']), 
                          num=params['kbin_num'] + 1)
        k_s = k_e[-1]/k_e[-2]
        k_c = k_e[:-1]*np.sqrt(k_s)
        return k_c

    def get_kbin_edges(self, params):
        k_e = np.logspace(np.log10(params['kbin_min']), 
                          np.log10(params['kbin_max']), 
                          num=params['kbin_num'] + 1)
        # using the logspace as the kbin edges, instead of kbin centre, that is
        # the same treatment as Eric's
        #k_s = k_e[-1]/k_e[-2]
        #k_e = np.append( k_e, k_e[-1]*k_s ) / (np.sqrt(k_s))
        return k_e

    def apply_cut_list(self, nmap):
        if len(self.params['cut_list']) != 0:
            for i in range(len(nmap)):
                nmap[i][self.params['cut_list']] = 0.
        return nmap

    def ps_estimate(self, params, imap_pair, nmap_pair, step, factors=None):
        
        print "Estimate power using map from"
        print len(imap_pair[0])*'-' +\
              "\n [%s\n *%s]\n x[%s\n *%s]\n"%\
              ( nmap_pair[0], imap_pair[0], nmap_pair[1], imap_pair[1]) +\
              len(imap_pair[0])*'-' + '\n'

        #imap1 = algebra.make_vect(algebra.load(imap_pair[0]))
        #imap2 = algebra.make_vect(algebra.load(imap_pair[1]))
        imap1 = [algebra.make_vect(algebra.load(x)) for x in imap_pair[0]]
        imap2 = [algebra.make_vect(algebra.load(x)) for x in imap_pair[1]]

        #box_path = []
        #if not os.path.exists(self.file_root + 'box/'):
        #    os.makedirs(self.file_root + 'box/')
        #map_name = imap_pair[0].split('/')[-1].split('.')[0]
        #box_path.append(self.file_root + 'box/' + map_name)
        #map_name = imap_pair[1].split('/')[-1].split('.')[0]
        #box_path.append(self.file_root + 'box/' + map_name)

        if factors != None:
            if factors[0] == 'square':
                imap1 = [x**2 for x in imap1]
                #imap1 = imap1**2
            else:
                imap1 = [x*factors[0] for x in imap1]
            if factors[1] == 'square':
                imap2 = [x**2 for x in imap2]
                #imap2 = imap2**2
            else:
                imap2 = [x*factors[1] for x in imap2]
                #imap2 *= factors[1]

        if nmap_pair[0] != None:
            #nmap1 = [np.ones_like(x) for x in imap1]
            nmap1 = [algebra.make_vect(algebra.load(x)) for x in nmap_pair[0]]
        else:
            print "No weighting given, use 1"
            nmap1 = [np.ones_like(x) for x in imap1]
        if nmap_pair[1] != None:
            #nmap2 = [np.ones_like(x) for x in imap2]
            nmap2 = [algebra.make_vect(algebra.load(x)) for x in nmap_pair[1]]
        else:
            print "No weighting given, use 1"
            nmap2 = [np.ones_like(x) for x in imap2]
        
        self.apply_cut_list(nmap1)
        self.apply_cut_list(nmap2)

        #if (step == 'si' and params['ps_type'] == 'auto'):
        #    import transf
        #    svd_path = params['svd_path']
        #    svd_file_list = params['svd_file']
        #    for ii in range(len(imap1)):
        #        map_temp = algebra.zeros_like(imap1[ii])
        #        map_numb = np.float(len(svd_file_list))
        #        for svd_file in svd_file_list:
        #            print "subtract modes from %s"%svd_path + svd_file
        #            svd_data  = h5py.File(svd_path + svd_file, 'r')
        #            map_temp += transf.subtract_svd_modes(
        #                    params, imap1[ii], svd_data, nmap1[ii], imap_pair[0][ii])
        #            svd_data.close()
        #        imap1[ii] = map_temp / map_numb
        #        del map_temp

        #if step == 'si' and params['ps_type'] == 'auto':
        #    import transf
        #    svd_path = params['svd_path']
        #    svd_file_list = params['svd_file']
        #    for ii in range(len(imap2)):
        #        map_temp = algebra.zeros_like(imap2[ii])
        #        map_numb = np.float(len(svd_file_list))
        #        for svd_file in svd_file_list:
        #            print "subtract modes from %s"%svd_path + svd_file
        #            svd_data  = h5py.File(svd_path + svd_file, 'r')
        #            map_temp += transf.subtract_svd_modes(
        #                    params, imap2[ii], svd_data, nmap2[ii], imap_pair[1][ii])
        #            svd_data.close()
        #        imap2[ii] = map_temp / map_numb
        #        del map_temp

        k_edges_p = self.get_kbin_edges(params)
        k_edges_v = self.get_kbin_edges(params)

        ps_box = functions.BOX(imap1, imap2, nmap1, nmap2)
        #ps_box.mapping_to_xyz(box_path)
        #ps_box.mapping_to_xyz(cube_force=None)
        ps_box.mapping_to_xyz(largeangle=False)
        #ps_box.mapping_to_xyz(largeangle=True, cube_force=0.4)
        if params['ps_type'] == 'auto':        
            print 'Estimating auto'
            ps_box.auto_estimate_ps_3d()
        else:
            ps_box.estimate_ps_3d()
        ps_box.convert_ps_to_unitless()
        ps_box.convert_3dps_to_2dps(k_edges_p, k_edges_v)
        ps_box.convert_3dps_to_2dps_dipole(k_edges_p, k_edges_v)
        ps_box.convert_3dps_to_1dps(k_edges_p)

        return ps_box.ps_2d, ps_box.kn_2d, ps_box.ps_1d, ps_box.kn_1d,\
                ps_box.ps_2d_kpp, ps_box.kn_2d_kpp, ps_box.ps_2d_kpn, ps_box.kn_2d_kpn,\
                ps_box.ps_1d_kpp, ps_box.kn_1d_kpp, ps_box.ps_1d_kpn, ps_box.kn_1d_kpn

    def process_maps(self, params, imap_list, nmap_list, tabs_list, step, 
                     factors=None):

        ps_2d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        kn_2d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        ps_2d_kpp = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        kn_2d_kpp = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        ps_2d_kpn = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        kn_2d_kpn = np.zeros(shape=(len(imap_list), 
                                params['kbin_num'], 
                                params['kbin_num']))
        ps_1d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        kn_1d = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        ps_1d_kpp = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        kn_1d_kpp = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        ps_1d_kpn = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        kn_1d_kpn = np.zeros(shape=(len(imap_list), 
                                params['kbin_num']))
        for i in range(len(imap_list)):
            ps_2d[i], kn_2d[i], ps_1d[i], kn_1d[i], \
                    ps_2d_kpp[i], kn_2d_kpp[i], ps_2d_kpn[i], kn_2d_kpn[i],\
                    ps_1d_kpp[i], kn_1d_kpp[i], ps_1d_kpn[i], kn_1d_kpn[i]=\
                self.ps_estimate(params, imap_list[i], nmap_list[i], step, factors)
            #if len(tabs_list) != 0:
            #    self.ps_result = h5py.File(self.file_root + 'ps_result.hd5', 'r+', driver='mpio', comm=MPI.COMM_WORLD)
            #    self.save_section(params, ps_2d[i], kn_2d[i], step, tabs_list[i])

        return ps_2d, kn_2d, ps_1d, kn_1d,\
                ps_2d_kpp, kn_2d_kpp, ps_2d_kpn, kn_2d_kpn,\
                ps_1d_kpp, kn_1d_kpp, ps_1d_kpn, kn_1d_kpn

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
        simroot = '/scratch2/p/pen/nluciw/parkes/simulations/bandpass/alpha20/ra%s/'%os.getenv('PKS_FIELD')
        rawroot = '/scratch2/p/pen/nluciw/parkes/simulations/corr_effic/rawsim/ra%s/'%os.getenv('PKS_FIELD')
        #rawroot = '/scratch2/p/pen/nluciw/parkes/simulations/corr_effic/rawsim/k_reduced/ra%s/'%os.getenv('PKS_FIELD')
        seleroot = '/scratch2/p/pen/nluciw/parkes/maps/2df/pks_%s_2dfmock_full/'%os.getenv('OPT_FIELD')
        if self.params['colour'] =='blue':
            colorroot = '/scratch2/p/pen/nluciw/parkes/simulations/color_gals/blue/%s/'%os.getenv('OPT_FIELD')
        elif self.params['colour'] =='red':
            colorroot = '/scratch2/p/pen/nluciw/parkes/simulations/color_gals/red/%s/'%os.getenv('OPT_FIELD')

        if map_set == 'si' and params['ps_type'] == 'cros':
            for i in range(params['sim_numb']):
                print 'prepping', i
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    #imaps_a = functions.get_mapdict(simroot, selection='pks')
                    #imaps_a = functions.get_mapdict(params['sim_root'][j], selection='pks')
                    #imaps_a = functions.get_mapdict(params['sim_root'][j], selection='delta')
                    imaps_a = functions.get_mapdict(colorroot, selection='optsim')
                    nmaps_a = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    #imaps_b = functions.get_mapdict(params['sim_root'][j], selection='optsim')
                    #imaps_b = functions.get_mapdict(params['sim_root'][j], selection='pks')
                    #imaps_b = functions.get_mapdict(params['opt_root'][j], selection='real')
                    imaps_b = functions.get_mapdict(colorroot, selection='optsim')
                    #nmaps_b = functions.get_mapdict(seleroot, selection='sele')
                    nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    #nmaps_b = functions.get_mapdict(params['gbt_root'][j])

                    #imapa_temp.append(imaps_a[1]['%d'%(i)])
                    imapa_temp.append(imaps_a[1]['%d'%(i)])
                    #imapa_temp.append(imaps_a[1]['map;%dmodes;%03d'%(ps_mode,i)])
                    imapb_temp.append(imaps_b[1]['%d'%(i)])
                    #imapb_temp.append(imaps_b)
                    #nmapa_temp.append(nmaps_a[1]['weight;10modes'])
                    nmapa_temp.append(nmaps_a)
                    nmapb_temp.append(nmaps_b)
 
                imap_list.append([imapa_temp, imapb_temp])
                nmap_list.append([nmapa_temp, nmapb_temp])
                #nmap_list.append([None,None])
                #nmap_list.append([nmapa_temp, None])
 
        elif map_set == 'ps' and params['ps_type'] == 'cros':

            imapa_temp = []
            imapb_temp = []
            nmapa_temp = []
            nmapb_temp = []
            for i in range(len(params['gbt_root'])):

                #imaps_a = functions.get_mapdict(params['gbt_root'][i])
                imaps_a = functions.get_mapdict(params['opt_root'][i], selection='real')
                #nmaps_a = functions.get_mapdict(params['gbt_root'][i])
                nmaps_a = functions.get_mapdict(params['opt_root'][i], selection='sele')
                #imap = functions.get_mapdict(params['gbt_root'][i])
                imap = functions.get_mapdict(params['opt_root'][i], selection='real')
                #nmap = functions.get_mapdict(seleroot, selection='sele')
                nmap = functions.get_mapdict(params['opt_root'][i], selection='sele')

                #imapa_temp.append(imaps_a[1]['map;%dmodes'%ps_mode])
                imapa_temp.append(imaps_a)
                imapb_temp.append(imap)
                #imapb_temp.append(imap[1]['map;%dmodes'%ps_mode])
                #nmapa_temp.append(nmaps_a[1]['weight;%dmodes'%ps_mode])
                nmapa_temp.append(nmap)
                nmapb_temp.append(nmap)

            imap_list = [[imapa_temp, imapb_temp],]
            nmap_list = [[nmapa_temp, nmapb_temp],]
            #nmap_list = [[None, None],]

        elif map_set == 'nl' and params['ps_type'] == 'cros':

            imapa_temp = []
            imapb_temp = []
            nmapa_temp = []
            nmapb_temp = []
            for i in range(len(params['gbt_root'])):

                imaps_a = functions.get_mapdict(params['gbt_root'][i])
                nmaps_a = functions.get_mapdict(params['gbt_root'][i])
                imaps_b = functions.get_mapdict(params['gbt_root'][i])
                nmaps_b = functions.get_mapdict(params['gbt_root'][i])

                imapa_temp.append(imaps_a[1]['map;%dmodes'%ps_mode])
                imapb_temp.append(imaps_b[1]['map;%dmodes'%ps_mode])
                nmapa_temp.append(nmaps_a[1]['weight;%dmodes'%ps_mode])
                nmapb_temp.append(nmaps_b[1]['weight;%dmodes'%ps_mode])

            imap_list = [[imapa_temp, imapb_temp],]
            nmap_list = [[nmapa_temp, nmapb_temp],]

        elif map_set == 'sn' and params['ps_type'] == 'cros':
            print 'prepare maps of mock catalog for cros power sepctrum short noise'
            for i in range(100,200):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    imaps_a = functions.get_mapdict(params['gbt_root'][j])
                    nmaps_a = functions.get_mapdict(params['gbt_root'][j])
                    imap = functions.get_mapdict(params['opt_root'][j], selection='mock')
                    #nmap = functions.get_mapdict(seleroot, selection='sele')
                    nmap = functions.get_mapdict(params['opt_root'][j], selection='sele')

                    #imapa_temp.append(imaps_a[1]['map;%dmodes'%ps_mode])
                    imapa_temp.append(imap[1]['%d'%i])
                    imapb_temp.append(imap[1]['%d'%i])
                    nmapa_temp.append(nmaps_a[1]['weight;%dmodes'%ps_mode])
                    nmapb_temp.append(nmap)
                    #imap = functions.get_mapdict(params['opt_root'][j])
                    #nmap = functions.get_mapdict(params['opt_root'][j], selection='separable')

                imap_list.append([imapa_temp, imapb_temp])
#                nmap_list.append([nmapa_temp, nmapb_temp])
                nmap_list.append([None, None])

        elif map_set == 'rf' and params['ps_type'] == 'cros':
            print 'prepare maps for reference calculation'
            for i in range(params['sim_numb']):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    imaps_a = functions.get_mapdict(params['sim_root'][j], selection='raw')
                    nmaps_a = functions.get_mapdict(params['gbt_root'][j])
                    imaps_b = functions.get_mapdict(params['sim_root'][j], selection='delta')
                    #nmaps_b = functions.get_mapdict(seleroot, selection='sele')
                    nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    imapa_temp.append(imaps_a[1]['%d'%(i)])
                    imapb_temp.append(imaps_b[1]['%d'%(i)])
                    nmapa_temp.append(nmaps_a[1]['weight;10modes'])
                    nmapb_temp.append(nmaps_b)

                imap_list.append([imapa_temp, imapb_temp])
                nmap_list.append([nmapa_temp, nmapb_temp])
                #nmap_list.append([nmapa_temp, None])
                #nmap_list.append([None,None])

        elif map_set == 'rb' and params['ps_type'] == 'cros':
            print 'prepare maps for reference calculation'
            for i in range(params['sim_numb']):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    #imaps_a = functions.get_mapdict(params['sim_root'][j], selection='beammeansub')
                    imaps_a = functions.get_mapdict(params['sim_root'][j], selection='degradebeam')
                    nmaps_a = functions.get_mapdict(params['gbt_root'][j])
                    imaps_b = functions.get_mapdict(params['sim_root'][j], selection='delta')
                    #nmaps_b = functions.get_mapdict(seleroot, selection='sele')
                    nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    imapa_temp.append(imaps_a[1]['%d'%(i)])
                    imapb_temp.append(imaps_b[1]['%d'%(i)])
                    nmapa_temp.append(nmaps_a[1]['weight;10modes'])
                    nmapb_temp.append(nmaps_b)

                imap_list.append([imapa_temp, imapb_temp])
                #nmap_list.append([nmapa_temp, nmapb_temp])
                nmap_list.append([None,None])

        elif map_set == 'tr' and params['ps_type'] == 'cros':
            print 'prepare maps for reference calculation'
            #for i in range(params['sim_numb']):
            for i in range(100):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    #imaps_a = functions.get_mapdict(params['beammeansub_root'][j], selection='beammeansub')
                    #imaps_a = functions.get_mapdict(params['sim_root'][j], selection='pks')
                    imaps_a = functions.get_mapdict(params['trans_root'][j])
                    #nmaps_a = functions.get_mapdict(params['trans_root'][j])
                    nmaps_a = functions.get_mapdict(params['gbt_root'][j])
                    imaps_b = functions.get_mapdict(params['sim_root'][j], selection='optsim')
                    #imaps_b = functions.get_mapdict(params['opt_root'][j], selection='mock')
                    #imaps_b = functions.get_mapdict(colorroot, selection='optsim')
                    #nmaps_b = functions.get_mapdict(seleroot, selection='sele')
                    nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    imapa_temp.append(imaps_a[1]['map;%dmodes;%03d'%(ps_mode,i)])
                    #imapa_temp.append(imaps_a[1]['%d'%(i)])
                    #imapa_temp.append(imaps_a[1]['map;%dmodes'%(ps_mode)])
                    imapb_temp.append(imaps_b[1]['%d'%(i)])
                    #nmapa_temp.append(nmaps_a[1]['weight;%dmodes;%03d'%(ps_mode,i)])
                    nmapa_temp.append(nmaps_a[1]['weight;%dmodes'%(ps_mode)])
                    nmapb_temp.append(nmaps_b)

                imap_list.append([imapa_temp, imapb_temp])
                nmap_list.append([nmapa_temp, nmapb_temp])
                #nmap_list.append([nmapa_temp, None])

        #elif map_set == 'tr' and params['ps_type'] == 'cros':
        #    print 'prepare maps for transfer function calculation'
        #    for i in range(params['sim_numb']):
        #        imapa_temp = []
        #        imapb_temp = []
        #        nmapa_temp = []
        #        nmapb_temp = []
        #        for j in range(len(params['sim_root'])):
        #            imaps_a = functions.get_mapdict(params['ssm_root'][j]%i)
        #            nmaps_a = functions.get_mapdict(params['gbt_root'][j])
        #            imaps_b = functions.get_mapdict(params['sim_root'][j], selection='delta')
        #            nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')
        #            imapa_temp.append(imaps_a[1]['map;%dmodes'%ps_mode])
        #            imapb_temp.append(imaps_b[1]['%d'%i])
        #            #nmapa_temp.append(nmaps_a[1]['weight;%dmodes'%ps_mode])
        #            nmapa_temp.append(nmaps_a[1]['weight;0modes'])
        #            nmapb_temp.append(nmaps_b)

        #        imap_list.append([imapa_temp, imapb_temp])
        #        nmap_list.append([nmapa_temp, nmapb_temp])
        elif map_set == 'rf' and params['ps_type'] == 'auto':
            print 'prepare maps for reference calculation'
            for i in range(params['sim_numb']):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    imaps_a = functions.get_mapdict(rawroot,selection='raw')
                    nmaps_a = functions.get_mapdict(params['gbt_root'][j])

                    imaps_b = functions.get_mapdict(rawroot,selection='raw')
                    nmaps_b = functions.get_mapdict(params['gbt_root'][j])
                    imapa_temp.append(imaps_a[1]['%d'%i])
                    imapb_temp.append(imaps_b[1]['%d'%i])
                    nmapa_temp.append(nmaps_a[1]['beam123_with_beam456;noise_inv;0modes'])
                    nmapb_temp.append(nmaps_b[1]['beam123_with_beam789;noise_inv;0modes'])

                imap_list.append([imapa_temp, imapb_temp])
                #nmap_list.append([nmapa_temp, nmapb_temp])
                nmap_list.append([None,None])

        elif map_set == 'tr' and params['ps_type'] == 'auto':
            print 'prepare maps for reference calculation'
            #for i in range(params['sim_numb']):
            for i in range(params['sim_numb']):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    #imaps_a = functions.get_mapdict(params['beammeansub_root'][j], selection='beammeansub')
                    #imaps_a = functions.get_mapdict(params['sim_root'][j], selection='degradebeam')
                    imaps_a = functions.get_mapdict(params['trans_root'][j])
                    #nmaps_a = functions.get_mapdict(params['trans_root'][j])
                    nmaps_a = functions.get_mapdict(params['gbt_root'][j])
                    imaps_b = functions.get_mapdict(rawroot, selection='raw')
                    nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    imapa_temp.append(imaps_a[1]['map;%dmodes;%03d'%(ps_mode,i)])
                    #imapa_temp.append(imaps_a[1]['map;%dmodes'%(ps_mode)])
                    imapb_temp.append(imaps_b[1]['%d'%i])
                    #nmapa_temp.append(nmaps_a[1]['weight;%dmodes;%03d'%(ps_mode,i)])
                    nmapa_temp.append(nmaps_a[1]['weight;%dmodes'%(ps_mode)])
                    nmapb_temp.append(nmaps_b)

                imap_list.append([imapa_temp, imapb_temp])
                nmap_list.append([nmapa_temp, None])


        #elif map_set == 'tr' and params['ps_type'] == 'auto':
        #    print 'prepare maps for transfer function calculation'
        #    for i in range(params['sim_numb']):
        #        imapa_temp = []
        #        imapb_temp = []
        #        nmapa_temp = []
        #        nmapb_temp = []
        #        for j in range(len(params['sim_root'])):
        #            imaps_a = functions.get_mapdict(params['sim_root'][j],selection='pks')
        #            nmaps_a = functions.get_mapdict(params['gbt_root'][j])
        #            imaps_b = functions.get_mapdict(params['sim_root'][j],selection='raw')
        #            nmaps_b = functions.get_mapdict(params['opt_root'][j],selection='sele')
        #            imapa_temp.append(imaps_a[1]['A;B;%d'%i])
        #            imapb_temp.append(imaps_b[1]['%d'%i])
        #            nmapa_temp.append(nmaps_a[1]['123;inv;456;%dmodes'%ps_mode])
        #            nmapb_temp.append(nmaps_b)

        #        imap_list.append([imapa_temp, imapb_temp])
        #        nmap_list.append([nmapa_temp, nmapb_temp])

        elif map_set == 'ps' and params['ps_type'] == 'auto':
            print 'prepare maps for power spectrum calculation'
            #imaps = functions.get_mapdict(params['gbt_root'])
            #section = ['A', 'B', 'C', 'D', 'E']
            section = ['123', '456', '789', '10to13']
            #section = ['A', 'B', 'C', 'D']
            for i in range(len(section)):
                for j in range(i+1, len(section)):
                    s1 = '%s'%(section[i])
                    s2 = '%s'%(section[j])
                    imapa_temp = []
                    imapb_temp = []
                    nmapa_temp = []
                    nmapb_temp = []
                    for k in range(len(params['gbt_root'])):
                        imaps = functions.get_mapdict(params['gbt_root'][k])
                        imapa_temp.append(imaps[1]['beam%s_with_beam%s;map;%dmodes'%(s1, s2, ps_mode)])
                        imapb_temp.append(imaps[1]['beam%s_with_beam%s;map;%dmodes'%(s2, s1, ps_mode)])
                        nmapa_temp.append(imaps[1]['beam%s_with_beam%s;noise_inv;%dmodes'%(s1, s2, ps_mode)])
                        nmapb_temp.append(imaps[1]['beam%s_with_beam%s;noise_inv;%dmodes'%(s2, s1,ps_mode)])
                    imap_list.append([imapa_temp, imapb_temp])
                    nmap_list.append([nmapa_temp, nmapb_temp])
                    tabs_list.append('%s%sx%s%s'%(section[i], section[j], 
                                                  section[j], section[i]))

        elif map_set == 'ns' and params['ps_type'] == 'auto':
            #imaps = functions.get_mapdict(params['gbt_root'])
            section = ['A', 'B', 'C', 'D']
            #section = ['A', 'B', 'C']
            #section = ['A', 'B', 'C', 'D', 'E']
            for i in range(len(section)):
                for j in range(len(section)):
                    if i == j:
                        continue
                    s = '%'%(section[i], section[j])
                    imapa_temp = []
                    imapb_temp = []
                    nmapa_temp = []
                    nmapb_temp = []
                    for k in range(len(params['gbt_root'])):
                        imaps = functions.get_mapdict(params['gbt_root'][k])
                        imapa_temp.append(imaps[1]['%s;map;%dmodes'%(s, ps_mode)])
                        imapb_temp.append(imaps[1]['%s;map;%dmodes'%(s, ps_mode)])
                        nmapa_temp.append(imaps[1]['%s;noise_inv;%dmodes'%(s, ps_mode)])
                        nmapb_temp.append(imaps[1]['%s;noise_inv;%dmodes'%(s, ps_mode)])
                    imap_list.append([imapa_temp, imapb_temp])
                    nmap_list.append([nmapa_temp, nmapb_temp])
                    tabs_list.append('%s%s'%(section[i], section[j]))

        elif map_set == 'si' and params['ps_type'] == 'auto':
            print 'prepare maps for simulation calculation'

            sections = ['123', '456', '789', '10to13']
            section = ['A', 'B', 'C', 'D']
            #section = ['A']
            for i in range(len(section)):
                for j in range(i+1, len(section)):
                    s1 = '%s'%(sections[i])
                    s2 = '%s'%(sections[j])
                    s3 = '%s'%sections[i]
                    s4 = '%s'%sections[j]
                    for k in range(params['sim_numb']):
                        print 'prepping', k
                        imapa_temp = []
                        imapb_temp = []
                        nmapa_temp = []
                        nmapb_temp = []
                        for l in range(len(params['sim_root'])):
                            imaps_a =functions.get_mapdict(params['sim_root'][l],selection='pks')
                            nmaps_a = functions.get_mapdict(params['gbt_root'][l])
                            imaps_b = functions.get_mapdict(params['sim_root'][l],selection='pks')
                            nmaps_b = functions.get_mapdict(params['gbt_root'][l])
                            imapa_temp.append(imaps_a[1]['%s;%s;%smodes;%d'%(s1,s2,ps_mode,k)])
                            imapb_temp.append(imaps_b[1]['%s;%s;%smodes;%d'%(s2,s1,ps_mode,k)])
                            #imapa_temp.append(imaps_a[1]['%d'%i])
                            #imapb_temp.append(imaps_b[1]['%d'%i])
                            nmapa_temp.append(nmaps_a[1]['beam%s_with_beam%s;noise_inv;%dmodes'%(s3,s4,ps_mode)])
                            nmapb_temp.append(nmaps_b[1]['beam%s_with_beam%s;noise_inv;%dmodes'%(s4,s3,ps_mode)])

                        imap_list.append([imapa_temp, imapb_temp])
                        #nmap_list.append([None,None])
                        nmap_list.append([nmapa_temp, nmapb_temp])

        elif map_set == 'ne' and params['ps_type'] == 'auto':
            print 'prepare maps for noise error calculation'
            #imaps = functions.get_mapdict(params['gbt_root'])
            #section = ['A', 'B', 'C', 'D', 'E']
            section = ['A', 'B', 'C', 'D']
            #section = ['A', 'B', 'C']
            for i in range(len(section)):
                for j in range(len(section)):
                    if i == j:
                        continue
                    s = '%s_with_%s'%(section[i], section[j])
                    imapa_temp = []
                    imapb_temp = []
                    nmapa_temp = []
                    nmapb_temp = []
                    for k in range(len(params['gbt_root'])):
                        imaps = functions.get_mapdict(params['gbt_root'][k])
                        imapa_temp.append(imaps[1]['%s;%s;inv;%dmodes'%(section[i], section[j], ps_mode)])
                        imapb_temp.append(imaps[1]['%s;%s;inv;%dmodes'%(section[i], section[j], ps_mode)])
                        nmapa_temp.append(None)
                        nmapb_temp.append(None)
                    imap_list.append([imapa_temp, imapb_temp])
                    nmap_list.append([nmapa_temp, nmapb_temp])
                    tabs_list.append('%s%s'%(section[i], section[j]))

        elif map_set == 'ps' and params['ps_type'] == '2df':
            print 'prepare maps of real catalog for 2df power sepctrum'
            imap_temp = []
            nmap_temp = []
            imap2_temp = []
            nmap2_temp = []

            for i in range(len(params['gbt_root'])):

                imap = functions.get_mapdict(params['opt_root'][i], selection='real')
                nmap = functions.get_mapdict(params['opt_root'][i], selection='sele')

                imap_temp.append(imap)
                nmap_temp.append(nmap)

                #second set of maps for colour-separation
                imap2 = functions.get_mapdict(params['opt_root'][i].replace('red', 'blue'), selection='real')
                nmap2 = functions.get_mapdict(params['opt_root'][i].replace('red', 'blue'), selection='sele')

                imap2_temp.append(imap2)
                nmap2_temp.append(nmap2)

            imap_list = [[imap_temp, imap_temp],]
            nmap_list = [[nmap_temp, nmap_temp],]

        elif map_set == 'sn' and params['ps_type'] == '2df':
            print 'prepare maps of mock catalog for 2df power sepctrum short noise'
            for i in range(params['sim_numb']):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    #imaps_a = functions.get_mapdict(params['sim_root'][j], selection='optsim')
                    imaps_a = functions.get_mapdict(params['opt_root'][j], selection='mock')
                    nmaps_a = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    #imaps_b = functions.get_mapdict(params['sim_root'][j], selection='optsim')
                    imaps_b = functions.get_mapdict(params['opt_root'][j], selection='mock')
                    nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')

                    imapa_temp.append(imaps_a[1]['%d'%(i)])
                    imapb_temp.append(imaps_b[1]['%d'%(i)])
                    nmapa_temp.append(nmaps_a)
                    nmapb_temp.append(nmaps_b)

                imap_list.append([imapa_temp, imapb_temp])
                nmap_list.append([nmapa_temp, nmapb_temp])

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
            for i in range(params['sim_numb']):
                imapa_temp = []
                imapb_temp = []
                nmapa_temp = []
                nmapb_temp = []
                for j in range(len(params['sim_root'])):
                    #imaps_a = functions.get_mapdict(params['sim_root'][j], selection='optsim')
                    imaps_a = functions.get_mapdict(params['sim_root'][j], selection='delta')
                    nmaps_a = functions.get_mapdict(params['opt_root'][j], selection='sele')
                    #imaps_b = functions.get_mapdict(params['sim_root'][j], selection='optsim')
                    imaps_b = functions.get_mapdict(params['sim_root'][j], selection='delta')
                    nmaps_b = functions.get_mapdict(params['opt_root'][j], selection='sele')

                    imapa_temp.append(imaps_a[1]['%d'%(i)])
                    imapb_temp.append(imaps_b[1]['%d'%(i)])
                    nmapa_temp.append(nmaps_a)
                    nmapb_temp.append(nmaps_b)

                imap_list.append([imapa_temp, imapb_temp])
                #nmap_list.append([nmapa_temp, nmapb_temp])
                nmap_list.append([None,None])
        else:
            print 'error: map_set error!'
            exit()
        return imap_list, nmap_list, tabs_list

    def save_section(self, params, ps_2d, kn_2d, step, tab):
        if step == 'ne' and params['ps_type'] == 'auto':
            ps_2d = np.sqrt(ps_2d)

        #k_bin = np.logspace(np.log10(params['kbin_min']), 
        #                    np.log10(params['kbin_max']), 
        #                    num=params['kbin_num'])
        k_bin = self.get_kbin_centr(params)

        k_axes_2d = ("k_p", "k_v")
        info_2d = {'axes': k_axes_2d, 'type': 'vect'}
        info_2d['k_p_delta']  = k_bin[1]/k_bin[0]
        info_2d['k_p_centre'] = k_bin[params['kbin_num']//2]
        info_2d['k_v_delta']  = k_bin[1]/k_bin[0]
        info_2d['k_v_centre'] = k_bin[params['kbin_num']//2]

        ps_2d = algebra.make_vect(ps_2d, axis_names=k_axes_2d)
        kn_2d = algebra.make_vect(kn_2d, axis_names=k_axes_2d)
        ps_2d.info = info_2d
        kn_2d.info = info_2d

        if params['ps_mode'] == None:
            file_name = '%s_%s_'%(params['ps_type'], step)
        else:
            file_name = '%s_%s_%dmode_'%(params['ps_type'], step, params['ps_mode'])

        #file_root = params['ps_root'] + params['ps_name'] + '/'
        #if not os.path.exists(file_root):
        #    os.makedirs(file_root)
        file_root = self.file_root

        algebra.save_h5(self.ps_result, file_name + '2dpow_%s'%tab, ps_2d)
        algebra.save_h5(self.ps_result, file_name + '2dkmn_%s'%tab, kn_2d)

    def save(self, params, ps_2d, kn_2d, ps_1d, kn_1d, 
            ps_2d_kpp, kn_2d_kpp, ps_2d_kpn, kn_2d_kpn, 
            ps_1d_kpp, kn_1d_kpp, ps_1d_kpn, kn_1d_kpn, 
            step):

        #print ps_2d.shape
        #print ps_2d.flatten().max(), ps_2d.flatten().min()
        ps_2d_mean = np.mean(ps_2d, axis=0)
        ps_2d_std  = np.std(ps_2d, axis=0)
        ps_1d_mean = np.mean(ps_1d, axis=0)
        ps_1d_std  = np.std(ps_1d, axis=0)

        kn_2d_mean = np.mean(kn_2d, axis=0)
        kn_2d_std  = np.std(ps_2d, axis=0)
        kn_1d_mean = np.mean(kn_1d, axis=0)
        kn_1d_std  = np.std(ps_1d, axis=0)

        #print ps_2d_kpp.shape
        ps_2d_kpp_mean = np.mean(ps_2d_kpp, axis=0)
        ps_2d_kpp_std  = np.std(ps_2d_kpp, axis=0)
        kn_2d_kpp_mean = np.mean(kn_2d_kpp, axis=0)
        kn_2d_kpp_std  = np.std(ps_2d_kpp, axis=0)

        ps_2d_kpn_mean = np.mean(ps_2d_kpn, axis=0)
        ps_2d_kpn_std  = np.std(ps_2d_kpn, axis=0)
        kn_2d_kpn_mean = np.mean(kn_2d_kpn, axis=0)
        kn_2d_kpn_std  = np.std(ps_2d_kpn, axis=0)

        ps_1d_kpp_mean = np.mean(ps_1d_kpp, axis=0)
        ps_1d_kpp_std  = np.std(ps_1d_kpp, axis=0)
        kn_1d_kpp_mean = np.mean(kn_1d_kpp, axis=0)
        kn_1d_kpp_std  = np.std(ps_1d_kpp, axis=0)

        ps_1d_kpn_mean = np.mean(ps_1d_kpn, axis=0)
        ps_1d_kpn_std  = np.std(ps_1d_kpn, axis=0)
        kn_1d_kpn_mean = np.mean(kn_1d_kpn, axis=0)
        kn_1d_kpn_std  = np.std(ps_1d_kpn, axis=0)

        if step == 'ps' and params['ps_type'] == 'auto':
            print ps_2d.shape[0]
            ps_2d_std /= sqrt(ps_2d.shape[0])
            ps_1d_std /= sqrt(ps_1d.shape[0])
        if step == 'ne' and params['ps_type'] == 'auto':
            ps_2d_mean = np.sqrt(ps_2d_mean)
            ps_2d_std = np.sqrt(ps_2d_std)
            ps_2d = np.sqrt(ps_2d)

        #k_bin = np.logspace(np.log10(params['kbin_min']), 
        #                    np.log10(params['kbin_max']), 
        #                    num=params['kbin_num'])
        k_bin = self.get_kbin_centr(params)

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

        #print k_axes_2d
        #print ps_2d_kpp_mean.shape
        #print ps_2d_kpp_mean

        ps_2d_kpp_mean = algebra.make_vect(ps_2d_kpp_mean, axis_names=k_axes_2d)
        kn_2d_kpp_mean = algebra.make_vect(kn_2d_kpp_mean, axis_names=k_axes_2d)
        ps_2d_kpp_std  = algebra.make_vect(ps_2d_kpp_std, axis_names=k_axes_2d)
        kn_2d_kpp_std  = algebra.make_vect(kn_2d_kpp_std, axis_names=k_axes_2d)

        ps_2d_kpp_mean.info = info_2d
        kn_2d_kpp_mean.info = info_2d
        ps_2d_kpp_std.info = info_2d
        kn_2d_kpp_std.info = info_2d

        ps_2d_kpn_mean = algebra.make_vect(ps_2d_kpn_mean, axis_names=k_axes_2d)
        kn_2d_kpn_mean = algebra.make_vect(kn_2d_kpn_mean, axis_names=k_axes_2d)
        ps_2d_kpn_std  = algebra.make_vect(ps_2d_kpn_std, axis_names=k_axes_2d)
        kn_2d_kpn_std  = algebra.make_vect(kn_2d_kpn_std, axis_names=k_axes_2d)

        ps_2d_kpn_mean.info = info_2d
        kn_2d_kpn_mean.info = info_2d
        ps_2d_kpn_std.info = info_2d
        kn_2d_kpn_std.info = info_2d

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

        ps_1d_kpp_mean = algebra.make_vect(ps_1d_kpp_mean, axis_names=k_axes_1d)
        kn_1d_kpp_mean = algebra.make_vect(kn_1d_kpp_mean, axis_names=k_axes_1d)
        ps_1d_kpp_std  = algebra.make_vect(ps_1d_kpp_std, axis_names=k_axes_1d)
        kn_1d_kpp_std  = algebra.make_vect(kn_1d_kpp_std, axis_names=k_axes_1d)

        ps_1d_kpp_mean.info = info_1d
        kn_1d_kpp_mean.info = info_1d
        ps_1d_kpp_std.info = info_1d
        kn_1d_kpp_std.info = info_1d

        ps_1d_kpn_mean = algebra.make_vect(ps_1d_kpn_mean, axis_names=k_axes_1d)
        kn_1d_kpn_mean = algebra.make_vect(kn_1d_kpn_mean, axis_names=k_axes_1d)
        ps_1d_kpn_std  = algebra.make_vect(ps_1d_kpn_std, axis_names=k_axes_1d)
        kn_1d_kpn_std  = algebra.make_vect(kn_1d_kpn_std, axis_names=k_axes_1d)

        ps_1d_kpn_mean.info = info_1d
        kn_1d_kpn_mean.info = info_1d
        ps_1d_kpn_std.info = info_1d
        kn_1d_kpn_std.info = info_1d

        if params['ps_mode'] == None:
            file_name = '%s_%s_'%(params['ps_type'], step)
        else:
            file_name = '%s_%s_%dmode_'%(params['ps_type'], step, params['ps_mode'])

        #file_root = params['ps_root'] + params['ps_name'] + '/'
        #if not os.path.exists(file_root):
        #    os.makedirs(file_root)
        file_root = self.file_root

        print file_root
        algebra.save_h5(self.ps_result, file_name + '2dpow', ps_2d_mean)
        algebra.save_h5(self.ps_result, file_name + '2derr', ps_2d_std )
        algebra.save_h5(self.ps_result, file_name + '2dkmn', kn_2d_mean)
        algebra.save_h5(self.ps_result, file_name + '1dpow', ps_1d_mean)
        algebra.save_h5(self.ps_result, file_name + '1derr', ps_1d_std )
        algebra.save_h5(self.ps_result, file_name + '1dkmn', kn_1d_mean)

        algebra.save_h5(self.ps_result, file_name + '2dpow_kpp', ps_2d_kpp_mean)
        algebra.save_h5(self.ps_result, file_name + '2derr_kpp', ps_2d_kpp_std )
        algebra.save_h5(self.ps_result, file_name + '2dkmn_kpp', kn_2d_kpp_mean)
        algebra.save_h5(self.ps_result, file_name + '2dpow_kpn', ps_2d_kpn_mean)
        algebra.save_h5(self.ps_result, file_name + '2derr_kpn', ps_2d_kpn_std )
        algebra.save_h5(self.ps_result, file_name + '2dkmn_kpn', kn_2d_kpn_mean)

        algebra.save_h5(self.ps_result, file_name + '1dpow_kpp', ps_1d_kpp_mean)
        algebra.save_h5(self.ps_result, file_name + '1derr_kpp', ps_1d_kpp_std )
        algebra.save_h5(self.ps_result, file_name + '1dkmn_kpp', kn_1d_kpp_mean)
        algebra.save_h5(self.ps_result, file_name + '1dpow_kpn', ps_1d_kpn_mean)
        algebra.save_h5(self.ps_result, file_name + '1derr_kpn', ps_1d_kpn_std )
        algebra.save_h5(self.ps_result, file_name + '1dkmn_kpn', kn_1d_kpn_mean)

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
        algebra.save_h5(self.ps_result, file_name + '2draw', ps_2d_each)

        ps_2d_each = algebra.make_vect(kn_2d, axis_names=k_axes_each_2d)
        ps_2d_each.info = info_each_2d
        algebra.save_h5(self.ps_result, file_name + '2draw_kmn', ps_2d_each)

        ps_1d_each = algebra.make_vect(ps_1d, axis_names=k_axes_each_1d)
        ps_1d_each.info = info_each_1d
        algebra.save_h5(self.ps_result, file_name + '1draw', ps_1d_each)

        ps_2d_each = algebra.make_vect(ps_2d_kpp, axis_names=k_axes_each_2d)
        ps_2d_each.info = info_each_2d
        algebra.save_h5(self.ps_result, file_name + '2draw_kpp', ps_2d_each)

        ps_2d_each = algebra.make_vect(kn_2d_kpp, axis_names=k_axes_each_2d)
        ps_2d_each.info = info_each_2d
        algebra.save_h5(self.ps_result, file_name + '2draw_kmn_kpp', ps_2d_each)

        ps_2d_each = algebra.make_vect(ps_2d_kpn, axis_names=k_axes_each_2d)
        ps_2d_each.info = info_each_2d
        algebra.save_h5(self.ps_result, file_name + '2draw_kpn', ps_2d_each)

        ps_2d_each = algebra.make_vect(kn_2d_kpn, axis_names=k_axes_each_2d)
        ps_2d_each.info = info_each_2d
        algebra.save_h5(self.ps_result, file_name + '2draw_kmn_kpn', ps_2d_each)

if __name__ == '__main__':
    import sys
