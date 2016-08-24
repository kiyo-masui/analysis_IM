#! /usr/bin/env python

import os 
import numpy as np
from core import algebra
from kiyopy import parse_ini
from mpi4py import MPI

params_init = {
              'sourceroot' : '',
              'sourcefile' : '',
              'targetroot' : '',
              'targetfile' : '',
              'mode_list' : [],
              'simnum' : 0,
              'sim_file': '',
              }
prefix = 'cb_'

class MakeCombined():
    def __init__(self, parameter_file=None, params_dict=None, feedback=0):
        self.params = parse_ini.parse(parameter_file, params_init, prefix=prefix)
    def mpiexecute(self, processes):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        params = self.params

        if params['simnum'] !=0:
            n_sim = params['simnum']
        else: 
            print 'MPI need total simulation number. '
            n_sim = 1
            exit()

        comm.barrier()
        if rank<n_sim:
            sim_list = range(rank, n_sim, size)
            for sim in sim_list:
                self.params['sourcefile'] = self.params['sim_file']%sim
                self.params['targetfile'] = self.params['sim_file']%sim
                print "RANK %d : combine file in [%s] to [%s]"\
                    %(rank, self.params['sourcefile'], self.params['targetfile'])
                self.execute(processes)
        comm.barrier()


    def execute(self, processes):
        params = self.params
        for mode in params['mode_list']:
            source_mapname = 'sec_%s_cleaned_clean_map_I_with_%s'+'_%dmodes.npy'%mode
            source_maproot = params['sourceroot'] + params['sourcefile']

            source_maplist = self.find_map_list(source_maproot, source_mapname)

            target_mapname = 'combined_clean_map_%dmodes.npy'%mode
            target_maproot = params['targetroot'] + params['targetfile']

            self.make_combined(source_maplist, target_maproot, target_mapname)

    def find_map_list(self, sourceroot, file):
        map_list = []
        for tag1 in ['A', 'B', 'C', 'D']:
            for tag2 in ['A', 'B', 'C', 'D']:
                if tag2 == tag1:
                    continue
                source_file_dir = sourceroot%(tag1) + file%(tag1, tag2)
                map_list.append(source_file_dir)
                print source_file_dir
        return map_list

    def make_combined(self, map_list, targetroot, file):
        imap_temp = algebra.make_vect(algebra.load(map_list[0]))
        cumulative_product = algebra.zeros_like(imap_temp)
        cumulative_weight  = algebra.zeros_like(imap_temp)
        for dir in map_list:
            imap = algebra.make_vect(algebra.load(dir))
            wmap = algebra.make_vect(algebra.load(dir.replace('clean_map', 'noise_inv')))
            cumulative_product += imap*wmap
            cumulative_weight  += wmap
        algebra.compressed_array_summary(cumulative_weight, "weight map")
        algebra.compressed_array_summary(cumulative_product, "product map")
    
        cumulative_weight[cumulative_weight < 1.e-20] = 0.
        cumulative_product[cumulative_weight < 1.e-20] = 0.
    
        cumulative_weight[cumulative_weight == 0] = np.inf
        cumulative_weight[np.isnan(cumulative_weight)] = np.inf
        newmap = cumulative_product / cumulative_weight
        cumulative_weight[np.isinf(cumulative_weight)] = 0.
    
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
    
        combined_map_file = targetroot + file
        combined_weight_file = targetroot + file.replace('map', 'weight')
        print combined_map_file
        print combined_weight_file
        algebra.save(combined_map_file, newmap)
        algebra.save(combined_weight_file, cumulative_weight)
        
    
