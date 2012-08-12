#! /usr/bin/env python 

import numpy as np
import core.algebra as algebra
import matplotlib.pyplot as plt


#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQV_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"
#mapdict = {
#          'imap' : maproot + "sec_I_modes_clean_map_I_with_Q_1modes.npy",
#          'qmap' : maproot + "sec_Q_modes_clean_map_I_with_I_1modes.npy",
#          'umap' : maproot + "sec_U_modes_clean_map_I_with_I_1modes.npy",
#          'name' : 'iqu_1modes'
#          }

class plotIQU(object):
    def __init__(self, mapdict):
        self.imap = algebra.load(mapdict['imap'])
        self.qmap = algebra.load(mapdict['qmap'])
        self.umap = algebra.load(mapdict['umap'])

        self.name = mapdict['name']

        self.imap = algebra.make_vect(self.imap)
        #ra = self.imap.get_axis('ra')
        #print ra
   
    def plotiqu(self):
        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')
        extent = (ra.min(), ra.max(), dec.min(), dec.max())
        for i in range(self.imap.shape[0]):
            plt.figure(figsize=(8,4))
            plt.imshow(self.imap[i].swapaxes(0,1), 
                       origin='lower', 
                       extent=extent,
                       cmap='gist_rainbow_r')
            plt.colorbar()
            plt.quiver(ra, dec, 
                       self.qmap[i].swapaxes(0,1), 
                       self.umap[i].swapaxes(0,1), 
                       color='w')

            plt.savefig('./png/%s_%d.png'%(self.name, i), format='png')
        

if __name__=='__main__':
    maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQV_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"
    mapdict = {}

    mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_1modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_1modes.npy"
    mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_1modes.npy"
    mapdict['name'] = 'iqu_1modes'
    a = plotIQU(mapdict)
    a.plotiqu()

    mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_2modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_2modes.npy"
    mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_2modes.npy"
    mapdict['name'] = 'iqu_2modes'
    a = plotIQU(mapdict)
    a.plotiqu()

    mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_3modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_3modes.npy"
    mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_3modes.npy"
    mapdict['name'] = 'iqu_3modes'
    a = plotIQU(mapdict)
    a.plotiqu()

    maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQV_legendre_modes_0gwj/IQUmap_clean_themselves/"

    mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_5modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_5modes.npy"
    mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_5modes.npy"
    mapdict['name'] = 'iqu_5modes'
    a = plotIQU(mapdict)
    a.plotiqu()

    mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_10modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_10modes.npy"
    mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_10modes.npy"
    mapdict['name'] = 'iqu_10modes'
    a = plotIQU(mapdict)
    a.plotiqu()




