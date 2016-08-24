#! /usr/bin/env python 

#import tables
import h5py
import numpy as np
from core import algebra as alg
import sys
import os

temp_root = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps/'
#temp_name = 'parkes_parallel_thread_ra33decn30_p08_568by106_beam1_clean_map_XX_1316.npy'
#temp_name = 'parkes_parallel_thread_ran18decn30_p08_488by106_beam1_clean_map_XX_1316.npy'
#temp_name = 'parkes_parallel_thread_ra165dec0_p08_440by136_beam1_clean_map_XX_1316.npy'
temp_name = 'parkes_parallel_thread_ra216dec0_p08_440by136_beam1_clean_map_XX_1316.npy'
temp = alg.make_vect(alg.load(temp_root + temp_name))


#noise1 = 'parkes_parallel_thread_ran18decn30_p08_488by106_beam' + os.getenv('BEAM') + '_noise_inv_YY_1316.hdf5'
#noise2 = 'parkes_parallel_thread_ran18decn30_p08_488by106_beam' + os.getenv('BEAM') + '_noise_inv_XX_1316.hdf5'

noise1 = 'parkes_parallel_thread_ra216dec0_p08_440by136_beam' + os.getenv('BEAM') + '_noise_inv_YY_1316.hdf5'
noise2 = 'parkes_parallel_thread_ra216dec0_p08_440by136_beam' + os.getenv('BEAM') + '_noise_inv_XX_1316.hdf5'

file_list = [noise1, noise2]

#file_list = [noise2]

#file_list = [noise1]

#file_list = [
#  'parkes_parallel_thread_ra33decn30_p08_568by106_beam9_noise_inv_YY_1316.hdf5'
#    ]


for file_name in file_list:
    file_root = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps/'
    print file_name
    file = h5py.File(file_root + file_name, mode='r')
    
    inv = file['inv_cov']
    
    shape = inv.shape
    
    for i in range(shape[0]):
        print '%03d .. '%i,
        temp[i,...] = np.diag(
            inv[i,...].reshape(shape[1]*shape[2],shape[3]*shape[4])).reshape(shape[1],shape[2])
        sys.stdout.flush()
    
    alg.save( file_root + '%s'%file_name.replace('inv', 'inv_diag').replace('hdf5', 'npy'), temp)
    print 
    
    file.close()
