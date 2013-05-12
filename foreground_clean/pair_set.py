r"""perform operations on sets of map pairs

Some improvements to consider:
    tag the outputs with the same 15hr_session etc as the data; right now these
    are fully specified by the directory
    pro to change: what if files move to a differently-named dir?
    con: what if tags grow sec_A_15hr_test1_4waysplit_etc.etc_etc._etc._with_B
    currently prefer more verbose directory names registered in the file db

    have the cleaning act on the weights -- this is probably only relevant when
    the full N^-1 rather than its diagonal is considered
"""
import os
import sys
import copy
import time
import scipy as sp
import numpy as np
from utils import file_tools as ft
from utils import data_paths as dp
from quadratic_products import corr_estimation as ce
from kiyopy import parse_ini
import kiyopy.utils
from core import algebra
from foreground_clean import map_pair
from foreground_clean import find_modes
from multiprocessing import Process, current_process
from utils import batch_handler
from mpi4py import MPI
# TODO: make map cleaning multiprocess; could also use previous cleaning, e.g.
# 5 modes to clean 10 modes = modes 5 to 10 (but no need to do this)
# TODO: move all magic strings to __init__ or params
# TODO: replace print with logging
# Special parameters that rarely need to get used
# no_weights: do not weight before finding nu-nu' covariance (False by default)
# SVD_root: use the SVD from another cleaning run on this run (None by default)
# regenerate_noise_inv: do not use memoize to save previous diag(N^-1) calc


params_init = {
               'SVD_root': None,
               'SVD_file': None,
               'output_root': "local_test",
               'pairlist' : None,
               'pairdict' : None,
               'map1': 'GBT_15hr_map',
               'map2': 'GBT_15hr_map',
               'noise_inv1': 'GBT_15hr_map',
               'noise_inv2': 'GBT_15hr_map',
               'calc_diagnal' : False,
               'simnum' : None,
               'simfile1': None,
               'simfile2': None,
               'sim_multiplier': 1.,
               'subtract_inputmap_from_sim': False,
               'subtract_sim_from_inputmap': False,
               'subtract_realmap_from_sim': False,
               'realmap_dir' : '',
               'freq_list1': (),
               'freq_list2': (),
               'tack_on': None,
               'convolve': True,
               'clip_weight_percent': None,
               'degrade_factor': 1.1,
               'weighted_SVD' : False,
               'factorizable_noise': True,
               'sub_weighted_mean': True,
               'regenerate_noise_inv': True,
               'modes': [10, 15],
               'no_weights': False,
               'save_section': True,
               }
prefix = 'fs_'


def wrap_find_weight(filename, regenerate=False, calc_diagnal=False):
    if not calc_diagnal:
        if regenerate:
            retval = find_weight(filename)
        else:
            retval = memoize_find_weight(filename)
    else:
        if regenerate:
            retval = find_weight_re_diagnal(filename)
        else:
            retval = memoize_find_weight_re_diagnal(filename)

    return batch_handler.repackage_kiyo(retval)


@batch_handler.memoize_persistent
def memoize_find_weight(filename):
    return find_weight(filename)


def find_weight(filename):
    r"""rather than read the full noise_inv and find its diagonal, cache the
    diagonal values.

    Note that the .info does not get shelved (class needs to be made
    serializeable). Return the info separately.
    """
    noise_inv_diag = algebra.make_vect(algebra.load(filename))

    return noise_inv_diag, noise_inv_diag.info

@batch_handler.memoize_persistent
def memoize_find_weight_re_diagnal(filename):
    return find_weight_re_diagnal(filename)


def find_weight_re_diagnal(filename):
    r"""rather than read the full noise_inv and find its diagonal, cache the
    diagonal values.

    Note that the .info does not get shelved (class needs to be made
    serializeable). Return the info separately.
    """
    print "loading noise: " + filename
    noise_inv = algebra.make_mat(algebra.open_memmap(filename, mode='r'))
    noise_inv_diag = noise_inv.mat_diag()

    return noise_inv_diag, noise_inv_diag.info

def inv_diag_noise(source_dict, target):
    if not os.path.exists(target):
        os.mkdir(target)
    for map_pair in source_dict.keys():
            fileroot = source_dict[map_pair]
            filename = fileroot.split('/')[-1]
            filename = filename.replace('diag', 'weight')

            noise_diag = algebra.make_vect(algebra.load(fileroot))
            noise_diag[noise_diag==0] = np.inf
            noise_inv_diag = 1./noise_diag
            algebra.save(target+filename, noise_inv_diag)

            source_dict[map_pair] = target + filename

    return source_dict

def diag_noise(source_dict, target):
    if not os.path.exists(target):
        os.mkdir(target)
    for map_pair in source_dict.keys():
            fileroot = source_dict[map_pair]
            filename = fileroot.split('/')[-1]
            filename = filename.replace('inv', 'weight')
            if not os.path.exists(target+filename):
                #os.symlink(fileroot, target + filename)
                noise_inv_diag, info = find_weight_re_diagnal(fileroot)
                algebra.save(target+filename, noise_inv_diag)

            source_dict[map_pair] = target + filename

    return source_dict

def mklink_for_mappair(source_dict, target):
    r"""to make a link from target direction to the source_dic direction.
    """
    if not os.path.exists(target):
        os.mkdir(target)
    for map_pair in source_dict.keys():
        for key in ['map1', 'map2']:
            fileroot = source_dict[map_pair][key]
            filename = fileroot.split('/')[-1]
            if not os.path.exists(target+filename):
                print target+filename
                os.symlink(fileroot, target + filename)
                os.symlink(fileroot+'.meta', target+filename+'.meta')

            source_dict[map_pair][key] = target + filename

        for key in ['noise_inv1', 'noise_inv2']:
            fileroot = source_dict[map_pair][key]
            filename = fileroot.split('/')[-1]
            filename = filename.replace('inv', 'weight')
            if not os.path.exists(target+filename):
                #os.symlink(fileroot, target + filename)
                noise_inv_diag, info = find_weight_re_diagnal(fileroot)
                algebra.save(target+filename, noise_inv_diag)

            source_dict[map_pair][key] = target + filename

    return source_dict

def extend_iqu_map(source_dict=None, target_dict=None, map_dict=None):
    if source_dict != None:
        imap = algebra.make_vect(algebra.load(source_dict['imap']))
        qmap = algebra.make_vect(algebra.load(source_dict['qmap']))
        umap = algebra.make_vect(algebra.load(source_dict['umap']))

        if source_dict.has_key('imap_weight'):
            imap_weight = algebra.make_vect(algebra.load(source_dict['imap_weight']))
            qmap_weight = algebra.make_vect(algebra.load(source_dict['qmap_weight']))
            umap_weight = algebra.make_vect(algebra.load(source_dict['umap_weight']))
        elif source_dict.has_key('imap_inv'):
            imap_weight, info = find_weight_re_diagnal(source_dict['imap_inv'])
            qmap_weight, info = find_weight_re_diagnal(source_dict['qmap_inv'])
            umap_weight, info = find_weight_re_diagnal(source_dict['umap_inv'])
        else:
            print 'Warning: no weight'
            imap_weight = algebra.ones_like(imap)
            qmap_weight = algebra.ones_like(imap)
            umap_weight = algebra.ones_like(imap)
    elif map_dict != None:
        imap = map_dict['imap']
        qmap = map_dict['qmap']
        umap = map_dict['umap']

        if 'imap_weight' in map_dict.keys():
            imap_weight = map_dict['imap_weight']
            qmap_weight = map_dict['qmap_weight']
            umap_weight = map_dict['umap_weight']
        else:
            print 'Warning: no weight'
            imap_weight = algebra.ones_like(imap)
            qmap_weight = algebra.ones_like(imap)
            umap_weight = algebra.ones_like(imap)
    else:
        print "Error: Can not find I Q U maps"
        exit()

    iqu = algebra.info_array(imap.tolist() + qmap.tolist() + umap.tolist())
    iqu = algebra.make_vect(iqu)
    iqu.info = imap.info
    iqu.copy_axis_info(imap)

    iqu_weight = algebra.info_array(imap_weight.tolist() + 
                                    qmap_weight.tolist() + 
                                    umap_weight.tolist())
    iqu_weight = algebra.make_vect(iqu_weight)
    iqu_weight.info = imap_weight.info
    iqu_weight.copy_axis_info(imap_weight)

    if target_dict != None:
        algebra.save(target_dict['map'], iqu)
        algebra.save(target_dict['weight'], iqu_weight)
    else:
        map_dict = {}
        map_dict['map']    = iqu
        map_dict['weight'] = iqu_weight
        return map_dict

def extend_iquv_map(source_dict=None, target_dict=None, map_dict=None):
    if source_dict != None:
        imap = algebra.make_vect(algebra.load(source_dict['imap']))
        qmap = algebra.make_vect(algebra.load(source_dict['qmap']))
        umap = algebra.make_vect(algebra.load(source_dict['umap']))
        vmap = algebra.make_vect(algebra.load(source_dict['vmap']))

        if source_dict.has_key('imap_weight'):
            imap_weight = algebra.make_vect(algebra.load(source_dict['imap_weight']))
            qmap_weight = algebra.make_vect(algebra.load(source_dict['qmap_weight']))
            umap_weight = algebra.make_vect(algebra.load(source_dict['umap_weight']))
            vmap_weight = algebra.make_vect(algebra.load(source_dict['vmap_weight']))
        elif source_dict.has_key('imap_inv'):
            imap_weight, info = find_weight_re_diagnal(source_dict['imap_inv'])
            qmap_weight, info = find_weight_re_diagnal(source_dict['qmap_inv'])
            umap_weight, info = find_weight_re_diagnal(source_dict['umap_inv'])
            vmap_weight, info = find_weight_re_diagnal(source_dict['vmap_inv'])
        else:
            print 'Warning: no weight'
            imap_weight = algebra.ones_like(imap)
            qmap_weight = algebra.ones_like(imap)
            umap_weight = algebra.ones_like(imap)
            vmap_weight = algebra.ones_like(imap)
    elif map_dict != None:
        imap = map_dict['imap']
        qmap = map_dict['qmap']
        umap = map_dict['umap']
        vmap = map_dict['vmap']

        if 'imap_weight' in map_dict.keys():
            imap_weight = map_dict['imap_weight']
            qmap_weight = map_dict['qmap_weight']
            umap_weight = map_dict['umap_weight']
            vmap_weight = map_dict['vmap_weight']
        else:
            print 'Warning: no weight'
            imap_weight = algebra.ones_like(imap)
            qmap_weight = algebra.ones_like(imap)
            umap_weight = algebra.ones_like(imap)
            vmap_weight = algebra.ones_like(imap)
    else:
        print "Error: Can not find I Q U V maps"
        exit()

    iquv = algebra.info_array(imap.tolist() + qmap.tolist() +
                              umap.tolist() + vmap.tolist())
    iquv = algebra.make_vect(iquv)
    iquv.info = imap.info
    iquv.copy_axis_info(imap)

    iquv_weight = algebra.info_array(imap_weight.tolist() + qmap_weight.tolist() + 
                                     umap_weight.tolist() + vmap_weight.tolist())
    iquv_weight = algebra.make_vect(iquv_weight)
    iquv_weight.info = imap_weight.info
    iquv_weight.copy_axis_info(imap_weight)

    if target_dict != None:
        algebra.save(target_dict['map'], iquv)
        algebra.save(target_dict['weight'], iquv_weight)
    else:
        map_dict = {}
        map_dict['map']    = iquv
        map_dict['weight'] = iquv_weight
        return map_dict

def divide_iqu_map(source_dict=None, target_dict=None, map_dict=None):
    if source_dict != None:
        iqu        = algebra.make_vect(algebra.load(source_dict['map']))
        iqu_weight = algebra.make_vect(algebra.load(source_dict['weight']))
    elif map_dict != None:
        iqu        = algebra.make_vect(map_dict['map'])
        iqu_weight = algebra.make_vect(map_dict['weight'])
    else:
        print "Error: Can not find iqu map"

    nfreq = iqu.shape[0]/3

    imap = algebra.make_vect(iqu[ 0*nfreq : 1*nfreq, ...])
    qmap = algebra.make_vect(iqu[ 1*nfreq : 2*nfreq, ...])
    umap = algebra.make_vect(iqu[ 2*nfreq : 3*nfreq, ...])

    imap.info = iqu.info
    qmap.info = iqu.info
    umap.info = iqu.info

    imap.copy_axis_info(iqu)
    qmap.copy_axis_info(iqu)
    umap.copy_axis_info(iqu)

    imap_weight = algebra.make_vect(iqu_weight[ 0*nfreq : 1*nfreq, ...])
    qmap_weight = algebra.make_vect(iqu_weight[ 1*nfreq : 2*nfreq, ...])
    umap_weight = algebra.make_vect(iqu_weight[ 2*nfreq : 3*nfreq, ...])

    imap_weight.info = iqu_weight.info
    qmap_weight.info = iqu_weight.info
    umap_weight.info = iqu_weight.info

    imap_weight.copy_axis_info(iqu_weight)
    qmap_weight.copy_axis_info(iqu_weight)
    umap_weight.copy_axis_info(iqu_weight)

    if target_dict != None:
        algebra.save(target_dict['imap'], imap)
        algebra.save(target_dict['qmap'], qmap)
        algebra.save(target_dict['umap'], umap)

        algebra.save(target_dict['imap_weight'], imap_weight)
        algebra.save(target_dict['qmap_weight'], qmap_weight)
        algebra.save(target_dict['umap_weight'], umap_weight)
    else:
        map_dict = {}
        map_dict['imap'] = imap
        map_dict['qmap'] = qmap
        map_dict['umap'] = umap
        map_dict['imap_weight'] = imap_weight
        map_dict['qmap_weight'] = qmap_weight
        map_dict['umap_weight'] = umap_weight
        return map_dict

def divide_iquv_map(source_dict=None, target_dict=None, map_dict=None):
    if source_dict != None:
        iquv        = algebra.make_vect(algebra.load(source_dict['map']))
        iquv_weight = algebra.make_vect(algebra.load(source_dict['weight']))
    elif map_dict != None:
        iquv        = algebra.make_vect(map_dict['map'])
        iquv_weight = algebra.make_vect(map_dict['weight'])
    else:
        print "Error: Can not find iquv map"

    nfreq = iquv.shape[0]/4

    imap = algebra.make_vect(iquv[ 0*nfreq : 1*nfreq, ...])
    qmap = algebra.make_vect(iquv[ 1*nfreq : 2*nfreq, ...])
    umap = algebra.make_vect(iquv[ 2*nfreq : 3*nfreq, ...])
    vmap = algebra.make_vect(iquv[ 3*nfreq : 4*nfreq, ...])

    imap.info = iquv.info
    qmap.info = iquv.info
    umap.info = iquv.info
    vmap.info = iquv.info

    imap.copy_axis_info(iquv)
    qmap.copy_axis_info(iquv)
    umap.copy_axis_info(iquv)
    vmap.copy_axis_info(iquv)

    imap_weight = algebra.make_vect(iquv_weight[ 0*nfreq : 1*nfreq, ...])
    qmap_weight = algebra.make_vect(iquv_weight[ 1*nfreq : 2*nfreq, ...])
    umap_weight = algebra.make_vect(iquv_weight[ 2*nfreq : 3*nfreq, ...])
    vmap_weight = algebra.make_vect(iquv_weight[ 3*nfreq : 4*nfreq, ...])

    imap_weight.info = iquv_weight.info
    qmap_weight.info = iquv_weight.info
    umap_weight.info = iquv_weight.info
    vmap_weight.info = iquv_weight.info

    imap_weight.copy_axis_info(iquv_weight)
    qmap_weight.copy_axis_info(iquv_weight)
    umap_weight.copy_axis_info(iquv_weight)
    vmap_weight.copy_axis_info(iquv_weight)

    if target_dict != None:
        algebra.save(target_dict['imap'], imap)
        algebra.save(target_dict['qmap'], qmap)
        algebra.save(target_dict['umap'], umap)
        algebra.save(target_dict['vmap'], vmap)

        algebra.save(target_dict['imap_weight'], imap_weight)
        algebra.save(target_dict['qmap_weight'], qmap_weight)
        algebra.save(target_dict['umap_weight'], umap_weight)
        algebra.save(target_dict['vmap_weight'], vmap_weight)
    else:
        map_dict = {}
        map_dict['imap'] = imap
        map_dict['qmap'] = qmap
        map_dict['umap'] = umap
        map_dict['vmap'] = vmap
        map_dict['imap_weight'] = imap_weight
        map_dict['qmap_weight'] = qmap_weight
        map_dict['umap_weight'] = umap_weight
        map_dict['vmap_weight'] = vmap_weight
        return map_dict

class PairSet():
    r"""Class to manage a set of map pairs
    """

    @batch_handler.log_timing
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        # recordkeeping
        self.pairs = {}
        self.pairs_parallel_track = {}
        self.pairlist = []
        self.datapath_db = dp.DataPath()

        self.params = params_dict
        if parameter_file:
            self.params = parse_ini.parse(parameter_file, params_init,
                                          prefix=prefix)

        self.freq_list1 = sp.array(self.params['freq_list1'], dtype=int)
        if len(self.params['freq_list2']) == 0:
            self.freq_list2 = self.freq_list1
        else:
            self.freq_list2 = sp.array(self.params['freq_list2'], dtype=int)


    @batch_handler.log_timing
    def mpiexecute(self, processes):
        r"""The MPI method is mainly build for cleanning map+sim"""
        
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

        output_root_tmp = params['output_root']
        simfile1_tmp = params['simfile1']
        comm.barrier()

        if rank<n_sim:
            sim_list = range(rank, n_sim, size)
            for sim in sim_list:
                self.params['output_root'] = output_root_tmp%sim
                self.params['simfile1'] = simfile1_tmp%sim
                print "RANK %d : sim%03d output %s"\
                      %(rank, sim, self.params['output_root'])
                print
                self.execute(processes)
                print "RANK %d : sim%03d job done"%(rank, sim)
        comm.barrier()

    @batch_handler.log_timing
    def execute(self, processes):
        r"""prepare direction"""
        #self.output_root = self.datapath_db.fetch(self.params['output_root'],
        #                                          tack_on=self.params['tack_on'])
        self.output_root = self.params['output_root']
        if not os.path.isdir(self.output_root):
            os.makedirs(self.output_root)


        if self.params['SVD_root']:
            if os.path.exists(self.params['SVD_root']):
                self.SVD_root = self.params['SVD_root']
            else:
                self.SVD_root = self.datapath_db.fetch(self.params['SVD_root'],
                                                   intend_write=True)
            print "WARNING: using %s to clean (intended?)" % self.SVD_root
        else:
            self.SVD_root = self.output_root

        # Write parameter file.
        parse_ini.write_params(self.params, self.output_root + 'params.ini',
                               prefix=prefix)

        r"""main call to execute the various steps in foreground removal"""
        self.load_pairs()

        self.preprocess_pairs()

        if self.params['weighted_SVD']:
            self.call_pairs("apply_map_weights")

        self.calculate_correlation()
        #self.calculate_svd()

        mode_list_stop = self.params['modes']
        mode_list_start = copy.deepcopy(self.params['modes'])
        mode_list_start[1:] = mode_list_start[:-1]

        #self.uncleaned_pairs = copy.deepcopy(self.pairs)
        for (n_modes_start, n_modes_stop) in zip(mode_list_start,
                                             mode_list_stop):
            self.subtract_foregrounds(n_modes_start, n_modes_stop)

            if self.params['weighted_SVD']:
                self.call_pairs("apply_map_weights")

            self.save_data(n_modes_stop)

            if self.params['weighted_SVD']:
                self.call_pairs("apply_map_weights")

            # NOTE: if you use this you also need to copy the parallel pairs!
            #self.pairs = copy.deepcopy(self.uncleaned_pairs)

    @batch_handler.log_timing
    def load_pairs(self):
        r"""load the set of map/noise pairs specified by keys handed to the
        database. This sets up operations on the quadratic product
            Q = map1^T noise_inv1 B noise_inv2 map2
        """
        par = self.params
        if not par['pairlist']:
            if par['calc_diagnal']:
                noise_inv_suffix = ";noise_inv"
            else:
                noise_inv_suffix = ";noise_weight"
            (self.pairlist, pairdict) = dp.cross_maps(par['map1'], par['map2'],
                                                 par['noise_inv1'],
                                                 par['noise_inv2'],
                                                 noise_inv_suffix=noise_inv_suffix,
                                                 verbose=False,
                                                 db_to_use=self.datapath_db)
        else:
            self.pairlist = par['pairlist']
            pairdict = par['pairdict']

        for pairitem in self.pairlist:
            pdict = pairdict[pairitem]
            print "-" * 80
            dp.print_dictionary(pdict, sys.stdout,
                                key_list=['map1', 'noise_inv1',
                                          'map2', 'noise_inv2'])

            # map1 & noise_inv1
            map1 = algebra.make_vect(algebra.load(pdict['map1']))
            if par['simfile1'] is not None:
                print "adding %s with multiplier %s" % (par['simfile1'],
                                                        par['sim_multiplier'])

                sim1 = algebra.make_vect(algebra.load(par['simfile1']))
                sim1 *= par['sim_multiplier']
            else:
                sim1 = algebra.zeros_like(map1)
            if not par['no_weights']:
                noise_inv1 = wrap_find_weight(pdict['noise_inv1'],
                                regenerate=par['regenerate_noise_inv'],
                                calc_diagnal = par['calc_diagnal'])
            else:
                noise_inv1 = algebra.ones_like(map1)

            # map2 & noise_inv2
            #if pairitem == 'I_with_E':
            if len(self.freq_list2) == 4*len(self.freq_list1):
                '''For IQUV case'''
                print 'Construct E map using I Q U V'
                iquvdict = {}
                iquvdict['imap'] = pdict['map2'].replace('_E', '_I')
                iquvdict['qmap'] = pdict['map2'].replace('_E', '_Q')
                iquvdict['umap'] = pdict['map2'].replace('_E', '_U')
                iquvdict['vmap'] = pdict['map2'].replace('_E', '_V')
                iquvdict['imap_weight'] = pdict['noise_inv2'].replace('_E', '_I')
                iquvdict['qmap_weight'] = pdict['noise_inv2'].replace('_E', '_Q')
                iquvdict['umap_weight'] = pdict['noise_inv2'].replace('_E', '_U')
                iquvdict['vmap_weight'] = pdict['noise_inv2'].replace('_E', '_V')
                map_dict = extend_iquv_map(source_dict=iquvdict)
                map2 = map_dict['map']
                noise_inv2 = map_dict['weight']

                sim2 = copy.deepcopy(sim1)
                map_dict = {}
                map_dict['imap'] = sim2
                map_dict['qmap'] = algebra.zeros_like(sim1)
                map_dict['umap'] = algebra.zeros_like(sim1)
                map_dict['vmap'] = algebra.zeros_like(sim1)
                map_dict = extend_iquv_map(map_dict=map_dict)
                sim2 = map_dict['map']
            elif len(self.freq_list2) == 3*len(self.freq_list1):
                '''For IQU case'''
                print 'Construct E map using I Q U'
                iquvdict = {}
                iquvdict['imap'] = pdict['map2'].replace('_E', '_I')
                iquvdict['qmap'] = pdict['map2'].replace('_E', '_Q')
                iquvdict['umap'] = pdict['map2'].replace('_E', '_U')
                iquvdict['imap_weight'] = pdict['noise_inv2'].replace('_E', '_I')
                iquvdict['qmap_weight'] = pdict['noise_inv2'].replace('_E', '_Q')
                iquvdict['umap_weight'] = pdict['noise_inv2'].replace('_E', '_U')
                map_dict = extend_iqu_map(source_dict=iquvdict)
                map2 = map_dict['map']
                noise_inv2 = map_dict['weight']

                sim2 = copy.deepcopy(sim1)
                map_dict = {}
                map_dict['imap'] = sim2
                map_dict['qmap'] = algebra.zeros_like(sim1)
                map_dict['umap'] = algebra.zeros_like(sim1)
                map_dict = extend_iqu_map(map_dict=map_dict)
                sim2 = map_dict['map']
            else:
                '''For common case'''
                map2 = algebra.make_vect(algebra.load(pdict['map2']))
                if par['simfile2'] is not None:
                    print "adding %s with multiplier %s" % (par['simfile2'],
                                                            par['sim_multiplier'])
                    sim2 = algebra.make_vect(algebra.load(par['simfile2']))
                    sim2 *= par['sim_multiplier']
                else:
                    sim2 = algebra.zeros_like(map2)
                if not par['no_weights']:
                    noise_inv2 = wrap_find_weight(pdict['noise_inv2'],
                                    regenerate=par['regenerate_noise_inv'],
                                    calc_diagnal = par['calc_diagnal'])
                else:
                    noise_inv2 = algebra.ones_like(map2)

            #if self.params['clip_weight_percent'] is not None:
            #    print "Note: your are clipping the weight maps"
            #    mask1 = self.define_weightmask(noise_inv1, 
            #                percentile=self.params['clip_weight_percent'])
            #    mask2 = self.define_weightmask(noise_inv2, 
            #                percentile=self.params['clip_weight_percent'])
            #    noise_inv1 = self.saturate_weight(noise_inv1, mask1)
            #    noise_inv2 = self.saturate_weight(noise_inv2, mask2)

            pair = map_pair.MapPair(map1 + sim1, map2 + sim2,
                                    noise_inv1, noise_inv2,
                                    self.freq_list1, self.freq_list2)
            pair.set_names(pdict['tag1'], pdict['tag2'])

            pair.params = self.params
            self.pairs[pairitem] = pair

            if par['subtract_inputmap_from_sim'] or \
               par['subtract_sim_from_inputmap']:
                if par['subtract_inputmap_from_sim']:
                    pair_parallel_track = map_pair.MapPair(map1, map2,
                                                  noise_inv1, noise_inv2,
                                                  self.freq_list1, self.freq_list2)

                if par['subtract_sim_from_inputmap']:
                    pair_parallel_track = map_pair.MapPair(sim1, sim2,
                                                  noise_inv1, noise_inv2,
                                                  self.freq_list1, self.freq_list2)

                pair_parallel_track.set_names(pdict['tag1'], pdict['tag2'])
                pair_parallel_track.params = self.params
                self.pairs_parallel_track[pairitem] = pair_parallel_track


    @batch_handler.log_timing
    def preprocess_pairs(self):
        r"""perform several preparation tasks on the data
        1. convolve down to a common beam
        2. make the noise factorizable
        3. subtract the weighted mean from each freq. slice
        """
        if self.params["convolve"]:
            self.call_pairs("degrade_resolution")

        if self.params["factorizable_noise"]:
            self.call_pairs("make_noise_factorizable")

        if self.params["sub_weighted_mean"]:
            self.call_pairs("subtract_weighted_mean")

    @batch_handler.log_timing
    def define_weightmask(self, input_weight, percentile=50.):
        # flatten to Ra/Dec
        input_weight = np.ma.array(input_weight)
        input_weight[input_weight==np.inf] = np.ma.masked
        input_weight[np.isnan(input_weight)] = np.ma.masked
        weight_2d = np.mean(input_weight, axis=0)
        print "define_weightmask: shape ", weight_2d.shape
        w_at_percentile = np.percentile(weight_2d, percentile)
    
        # return a boolean mask (here True = masked)
        return (weight_2d <= w_at_percentile)
    
    @batch_handler.log_timing
    def saturate_weight(self, input_weight, weightmask):
        r"""replace all weights above a given mask with the average in the
        mask. Note that this clobbers the input array"""
    
        for freq_index in range(input_weight.shape[0]):
            mweight = np.ma.array(input_weight[freq_index, ...], mask=weightmask)
            meanslice = np.ma.mean(mweight)
            input_weight[freq_index, np.logical_not(weightmask)] = meanslice
    
        return input_weight


    @batch_handler.log_timing
    def calculate_correlation(self):
        r"""Note that multiprocessing's map() is more elegant than Process,
        but fails for handing in complex map_pair objects
        """
        process_list = []
        for pairitem in self.pairlist:
            filename = self.output_root
            filename_svd = filename + "SVD_pair_%s.pkl" % pairitem

            map1 = copy.deepcopy(np.array(self.pairs[pairitem].map1))
            map2 = copy.deepcopy(np.array(self.pairs[pairitem].map2))
            weight1 = copy.deepcopy(np.array(self.pairs[pairitem].noise_inv1))
            weight2 = copy.deepcopy(np.array(self.pairs[pairitem].noise_inv2))
            freqs1 = copy.deepcopy(self.pairs[pairitem].freq1)
            freqs2 = copy.deepcopy(self.pairs[pairitem].freq2)

            if self.params['clip_weight_percent'] is not None:
                print "Note: your are clipping the weight maps"
                percentile = self.params['clip_weight_percent']
                mask1   = self.define_weightmask(weight1, percentile=percentile)
                mask2   = self.define_weightmask(weight2[:weight1.shape[0], ...], 
                                                 percentile=percentile)
                weight1 = self.saturate_weight(weight1, mask1)
                weight2 = self.saturate_weight(weight2, mask1)

            (freq_cov, counts) = find_modes.freq_covariance(map1, map2,
                                        weight1, weight2,
                                        freqs1, freqs2,
                                        no_weight=self.params['weighted_SVD'])

            n_modes  = min(len(freqs1),len(freqs2))
            svd_info = find_modes.get_freq_svd_modes(freq_cov, n_modes)
            ft.save_pickle(svd_info, filename_svd)

        #    multi = Process(target=wrap_corr, args=([self.pairs[pairitem],
        #                    filename]), name=pairitem)

        #    process_list.append(multi)

        #    multi.start()

        #for process in process_list:
        #    process.join()

    @batch_handler.log_timing
    def subtract_foregrounds(self, n_modes_start, n_modes_stop):
        for pairitem in self.pairlist:
            if self.params['SVD_file'] != None:
                filename_svd = "%s/%s" % (self.SVD_root, self.params['SVD_file'])
            else:
                filename_svd = "%s/SVD_pair_%s.pkl" % (self.SVD_root, pairitem)
            print "subtracting %d to %d modes from %s using %s" % (n_modes_start, \
                                                                n_modes_stop, \
                                                                pairitem, \
                                                                filename_svd)

            # svd_info: 0 is vals, 1 is modes1 (left), 2 is modes2 (right)
            svd_info = ft.load_pickle(filename_svd)

            self.pairs[pairitem].subtract_frequency_modes(
                                    svd_info[1][n_modes_start:n_modes_stop],
                                    svd_info[2][n_modes_start:n_modes_stop])

            if self.params['subtract_inputmap_from_sim'] or \
               self.params['subtract_sim_from_inputmap']:
                self.pairs_parallel_track[pairitem].subtract_frequency_modes(
                                        svd_info[1][n_modes_start:n_modes_stop],
                                        svd_info[2][n_modes_start:n_modes_stop])

    @batch_handler.log_timing
    def save_data(self, n_modes):
        prodmap_list = []
        weight_list = []

        n_modes = "%dmodes" % n_modes
        for pairitem in self.pairlist:
            pair = self.pairs[pairitem]
            (tag1, tag2) = (pair.map1_name, pair.map2_name)
            clnoise = "cleaned_noise_inv"
            map1_file = "%s/sec_%s_cleaned_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, tag2, n_modes)
            noise_inv1_file = "%s/sec_%s_%s_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, clnoise, tag2, n_modes)
            modes1_file = "%s/sec_%s_modes_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag1, tag2, n_modes)

            #if pair.map1.shape == pair.map2.shape:
            map2_file = "%s/sec_%s_cleaned_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, tag1, n_modes)
            noise_inv2_file = "%s/sec_%s_%s_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, clnoise, tag1, n_modes)
            modes2_file = "%s/sec_%s_modes_clean_map_I_with_%s_%s.npy" % \
                            (self.output_root, tag2, tag1, n_modes)

            if self.params['subtract_inputmap_from_sim'] or \
               self.params['subtract_sim_from_inputmap']:
                map1 = pair.map1 - self.pairs_parallel_track[pairitem].map1
                map2 = pair.map2 - self.pairs_parallel_track[pairitem].map2
            elif self.params['subtract_realmap_from_sim']:
                if not os.path.exists(self.params['realmap_dir']):
                    print "Error: Real map directory does not exists"
                    exit()
                else:
                    realmap_file = "%s/sec_%s_cleaned_clean_map_I_with_%s_%s.npy"%\
                                   (self.params['realmap_dir'], tag1, tag2, n_modes)
                    realmap = algebra.make_vect(algebra.load(realmap_file))
                    print "Subtract realmap from result"
                    map1 = copy.deepcopy(pair.map1) - realmap
                    map2 = copy.deepcopy(pair.map2)
                    if map2.shape == map1.shape:
                        map2 -= realmap
            else:
                map1 = copy.deepcopy(pair.map1)
                map2 = copy.deepcopy(pair.map2)

            prodmap_list.append(map1 * pair.noise_inv1)
            prodmap_list.append(map2 * pair.noise_inv2)
            weight_list.append(pair.noise_inv1)
            weight_list.append(pair.noise_inv2)

            if self.params['save_section']:
                algebra.save(map1_file, map1)
                algebra.save(noise_inv1_file, pair.noise_inv1)
                algebra.save(modes1_file, pair.left_modes)

                if pair.map1.shape == pair.map2.shape:
                    algebra.save(map2_file, map2)
                    algebra.save(noise_inv2_file, pair.noise_inv2)
                    algebra.save(modes2_file, pair.right_modes)

            #if map2.shape[0] == 3*map1.shape[0]:
            #    #source_dict = {}
            #    #source_dict['map'] = map2_file
            #    #source_dict['weight'] = noise_inv2_file
            #    map_dict = {}
            #    map_dict['map'] = map2
            #    map_dict['weight'] = pair.noise_inv2
            #    target_dict = {}
            #    target_dict['imap'] = map2_file.replace('_'+tag2, '_'+tag2+'_I')
            #    target_dict['qmap'] = map2_file.replace('_'+tag2, '_'+tag2+'_Q')
            #    target_dict['umap'] = map2_file.replace('_'+tag2, '_'+tag2+'_U')
            #    target_dict['imap_weight'] =\
            #                    noise_inv2_file.replace('_'+tag2, '_'+tag2+'_I')
            #    target_dict['qmap_weight'] =\
            #                    noise_inv2_file.replace('_'+tag2, '_'+tag2+'_Q')
            #    target_dict['umap_weight'] =\
            #                    noise_inv2_file.replace('_'+tag2, '_'+tag2+'_U')
            #    divide_iqu_map(map_dict=map_dict, target_dict=target_dict)

        if map1.shape != map2.shape:
            print "Shape of map1 and map2 are different, can not get combined map."
        else:
            cumulative_product = algebra.zeros_like(prodmap_list[0])
            cumulative_weight = algebra.zeros_like(prodmap_list[0])
            for mapind in range(0, len(prodmap_list)):
                cumulative_product += prodmap_list[mapind]
                cumulative_weight += weight_list[mapind]

            algebra.compressed_array_summary(cumulative_weight, "weight map")
            algebra.compressed_array_summary(cumulative_product, "product map")

            cumulative_weight[cumulative_weight < 1.e-20] = 0.
            cumulative_product[cumulative_weight < 1.e-20] = 0.

            newmap = cumulative_product / cumulative_weight

            # if the new map is nan or inf, set it and the wieghts to zero
            nan_array = np.isnan(newmap)
            newmap[nan_array] = 0.
            cumulative_product[nan_array] = 0.
            cumulative_weight[nan_array] = 0.
            inf_array = np.isinf(newmap)
            newmap[inf_array] = 0.
            cumulative_product[inf_array] = 0.
            cumulative_weight[inf_array] = 0.
            algebra.compressed_array_summary(newmap, "new map")
            algebra.compressed_array_summary(cumulative_product,"final map * weight")
            algebra.compressed_array_summary(cumulative_weight, "final weight map")

            combined = "combined_clean"
            combined_map_file = "%s/%s_map_%s.npy" % \
                                (self.output_root, combined, n_modes)
            combined_weight_file = "%s/%s_weight_%s.npy" % \
                                (self.output_root, combined, n_modes)
            combined_product_file = "%s/%s_product_%s.npy" % \
                                (self.output_root, combined, n_modes)
            combined_ones_file = "%s/%s_ones_%s.npy" % \
                                (self.output_root, combined, n_modes)

            algebra.save(combined_map_file, newmap)
            algebra.save(combined_product_file, cumulative_product)
            algebra.save(combined_weight_file, cumulative_weight)
            algebra.save(combined_ones_file, algebra.ones_like(newmap))

    # Service functions ------------------------------------------------------
    def call_pairs(self, call):
        r"""call some operation on the map pairs"""
        for pairitem in self.pairlist:
            pair = self.pairs[pairitem]
            try:
                method_to_call = getattr(pair, call)
            except AttributeError:
                print "ERROR: %s missing call %s" % (pairitem, call)
                sys.exit()

            print "calling %s() on pair %s" % (call, pairitem)
            method_to_call()

        if self.params['subtract_inputmap_from_sim'] or \
           self.params['subtract_sim_from_inputmap']:
            for pairitem in self.pairlist:
                pair_parallel_track = self.pairs_parallel_track[pairitem]
                try:
                    method_to_call_parallel_track = getattr(pair_parallel_track, call)
                except AttributeError:
                    print "ERROR: %s missing call %s" % (pairitem, call)
                    sys.exit()

                print "calling %s() on pair %s" % (call, pairitem)
                method_to_call_parallel_track()


def wrap_corr(pair, filename):
    r"""Do the correlation for a map_pair `pair`.
    Correlations in the `pair` (map_pair type) are saved to `filename`
    """
    name = current_process().name
    (freq_cov, counts) = pair.freq_covariance()
    ft.save_pickle((freq_cov, counts), filename)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        PairSet(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'
