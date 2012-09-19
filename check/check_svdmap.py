#! /usr/bin/env python 

import os
#from plotting import plot_cube
import numpy as np
from core import algebra
import matplotlib.pyplot as plt
from mkpower import functions


#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_5goodmode_weightedjump/"
#mapname = "sec_A_modes_clean_map_I_with_B_100modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gwj'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_unchanged/"
#mapname = "sec_A_modes_clean_map_I_with_B_25modes.npy"
#tag = mapname.split('.')[0]

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_5gwj/"
#mapname = "sec_A_modes_clean_map_I_with_B_25modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gwj'
#
#
##movieroot = "/mnt/raid-project/gmrt/ycli/movies/"
#movieroot = "/cita/d/www/home/ycli/movies/"
#cube_frame_dir = "/mnt/raid-project/gmrt/ycli/cube_frames/"
#plot_cube.make_cube_movie(maproot + mapname, "Temperature (mK)", cube_frame_dir,
#                          sigmarange=-1., outputdir=movieroot, multiplier=1000.,
#                          transverse=False, tag=tag, convolve=False, logscale=True,
#                          physical=False)


#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IxI5svd/"
#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"
#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_themselves/"
#mapfile = "sec_U_modes_clean_map_I_with_I_5modes.npy"

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_1gwj/IxI5svd/"
maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/Parkes_sept12-14_west_legendre_modes_0gwj/mapmode_map/"
mapfile = "sec_I_modes_clean_map_I_with_I_5modes.npy"

map = algebra.load(maproot + mapfile)
map = algebra.make_vect(map)
print map.shape

for i in range(map.shape[0]):
    #plt.figure(figsize=(8, 4))
    plt.figure(figsize=(34, 8))
    plt.imshow(map[i].swapaxes(0,1), interpolation='nearest', origin='lower')
    #plt.imshow(map[mapidex].swapaxes(0,1))
    plt.colorbar()
    
    pngname = mapfile.replace(mapfile.split('_')[-1], '%dmodes')
    pngname = './png/' + pngname + '.png'
    plt.savefig(pngname % (i), format='png')
    
