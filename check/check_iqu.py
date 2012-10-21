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
    def __init__(self, mapdict, imap=None, qmap=None, umap=None):
        if imap is None:
            self.imap = algebra.load(mapdict['imap'])
            self.qmap = algebra.load(mapdict['qmap'])
            self.umap = algebra.load(mapdict['umap'])
        else:
            self.imap = imap
            self.qmap = qmap
            self.umap = umap

        self.angle = np.arctan(self.umap/self.qmap)*180/np.pi
        self.angle[self.qmap<0] += 180
        self.angle *= 0.5

        self.name = mapdict['name']

        self.imap = algebra.make_vect(self.imap)
        #ra = self.imap.get_axis('ra')
        #print ra
   
    def plotiqu(self, num=None):
        ra = self.imap.get_axis('ra')
        dec= self.imap.get_axis('dec')
        dra = ra[1]-ra[0]
        ddec = dec[1]-dec[0]
        ra = np.resize(ra-0.5*dra,  ra.shape[0]+1)
        dec= np.resize(dec-0.5*ddec,dec.shape[0]+1)
        ra[-1] = ra[-2]+dra
        dec[-1]= dec[-2]+ddec
        extent = (ra.min(), ra.max(), dec.min(), dec.max())
        if num==None:
            for i in range(self.imap.shape[0]):
                plt.figure(figsize=(16,8))
                #plt.imshow(self.imap[i].swapaxes(0,1), 
                #           origin='lower', 
                #           extent=extent,
                #           cmap='gist_rainbow_r')
                map = self.imap[i].swapaxes(0,1)
                max = map.flatten().max()
                min = map.flatten().min()
                if np.fabs(min) > np.fabs(max):
                    max = np.fabs(min)
                plt.pcolor(ra, dec, map, vmax=max, vmin=-max)
                plt.colorbar()
                plt.quiver(ra[:-1]+0.5*dra, dec[:-1]+0.5*ddec, 
                           self.qmap[i].swapaxes(0,1), 
                           self.umap[i].swapaxes(0,1), 
                           angles = self.angle[i].swapaxes(0,1),
                           pivot = 'middle',
                           headwidth = 1,
                           color='k')

                plt.xlim(ra.min(), ra.max())
                plt.ylim(dec.min(), dec.max())
                plt.xlabel('RA [deg]')
                plt.ylabel('Dec[deg]')
                plt.savefig('./png/%s_%d.png'%(self.name, i), format='png')
        else:
            i = num
            plt.figure(figsize=(16,8))
            #plt.imshow(self.imap[i].swapaxes(0,1), 
            #           origin='lower', 
            #           extent=extent,)
            map = self.imap[i].swapaxes(0,1)
            #max = map.flatten().max()
            #min = map.flatten().min()
            max = 0.03
            min = -0.03
            if np.fabs(min) > np.fabs(max):
                max = np.fabs(min)
            else:
                min = -np.fabs(max)
            plt.pcolor(ra, dec, map, vmax=max, vmin=-max)
            plt.colorbar()
            plt.quiver(ra[:-1]+0.5*dra, dec[:-1]+0.5*ddec,
                       self.qmap[i].swapaxes(0,1), 
                       self.umap[i].swapaxes(0,1), 
                       angles = self.angle[i].swapaxes(0,1),
                       pivot = 'middle',
                       headwidth = 1,
                       color='k')

            plt.xlim(ra.min(), ra.max())
            plt.ylim(dec.min(), dec.max())
            plt.xlabel('RA [deg]')
            plt.ylabel('Dec[deg]')
            plt.savefig('./png/%s_%d.png'%(self.name, i), format='png')

        

if __name__=='__main__':
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_themselves/"
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"
    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_1gwj/IQUmap_clean_withIxIsvd/"
    maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IQU_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"
    mapdict = {}

    #mapdict['imap'] = maproot + "fir_1hr_41-18_avg_fdgp_clean_map_I_800.npy"
    #mapdict['qmap'] = maproot + "fir_1hr_41-18_avg_fdgp_clean_map_Q_800.npy"
    #mapdict['umap'] = maproot + "fir_1hr_41-18_avg_fdgp_clean_map_U_800.npy"
    #mapdict['name'] = '1hr_noise_iqu_cleaned_clean_0'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_0modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_0modes.npy"
    mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_0modes.npy"
    mapdict['name'] = '1hr_iqu_cleaned_clean_0'
    a = plotIQU(mapdict)
    a.plotiqu(150)

    mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_1modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_1modes.npy"
    mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_1modes.npy"
    mapdict['name'] = '1hr_iqu_cleaned_clean_1'
    a = plotIQU(mapdict)
    a.plotiqu(150)

    mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_2modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_2modes.npy"
    mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_2modes.npy"
    mapdict['name'] = '1hr_iqu_cleaned_clean_2'
    a = plotIQU(mapdict)
    a.plotiqu(150)

    maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_IQU_legendre_modes_0gwj/IQUmap_clean_themselves/"

    mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_I_5modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_5modes.npy"
    mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_5modes.npy"
    mapdict['name'] = '1hr_iqu_cleaned_clean_5'
    a = plotIQU(mapdict)
    a.plotiqu(150)

    mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_I_10modes.npy"
    mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_10modes.npy"
    mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_10modes.npy"
    mapdict['name'] = '1hr_iqu_cleaned_clean_10'
    a = plotIQU(mapdict)
    a.plotiqu(150)

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_0modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_0modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_0modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_0'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_1modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_1modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_1modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_1'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_2modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_2modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_2modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_2'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_3modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_3modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_3modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_3'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    #mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_1modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_1modes.npy"
    #mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_1modes.npy"
    #mapdict['name'] = 'iqu_1modes'
    #a = plotIQU(mapdict)
    #a.plotiqu()

    #mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_2modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_2modes.npy"
    #mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_2modes.npy"
    #mapdict['name'] = 'iqu_2modes'
    #a = plotIQU(mapdict)
    #a.plotiqu()

    #mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_3modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_3modes.npy"
    #mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_3modes.npy"
    #mapdict['name'] = 'iqu_3modes'
    #a = plotIQU(mapdict)
    #a.plotiqu()

    #----------------------

    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_themselves/"

    #mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_5modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_5modes.npy"
    #mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_5modes.npy"
    #mapdict['name'] = 'iqu_5modes'
    #a = plotIQU(mapdict)
    #a.plotiqu()

    #mapdict['imap'] = maproot + "sec_I_modes_clean_map_I_with_Q_10modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_modes_clean_map_I_with_I_10modes.npy"
    #mapdict['umap'] = maproot + "sec_U_modes_clean_map_I_with_I_10modes.npy"
    #mapdict['name'] = 'iqu_10modes'
    #a = plotIQU(mapdict)
    #a.plotiqu()

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_5modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_5modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_5modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_8'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_20modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_20modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_20modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_23'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_I_20modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_20modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_20modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_23_II'
    #a = plotIQU(mapdict)
    #a.plotiqu(150)

    #----------------------

    #maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQU_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"

    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_1modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_1modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_1modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_2'
    #a = plotIQU(mapdict)
    #imap1 = a.imap
    #qmap1 = a.qmap
    #umap1 = a.umap
    #mapdict['imap'] = maproot + "sec_I_cleaned_clean_map_I_with_Q_2modes.npy"
    #mapdict['qmap'] = maproot + "sec_Q_cleaned_clean_map_I_with_I_2modes.npy"
    #mapdict['umap'] = maproot + "sec_U_cleaned_clean_map_I_with_I_2modes.npy"
    #mapdict['name'] = 'iqu_cleaned_clean_2'
    #b = plotIQU(mapdict)
    #imap2 = b.imap
    #qmap2 = b.qmap
    #umap2 = b.umap


    #imap = imap1 - imap2
    #qmap = qmap1 - qmap2
    #umap = umap1 - umap2

    #mapdict['name'] = 'iqu_cleaned_clean_2-1'
    #c = plotIQU(mapdict, imap, qmap, umap)

    #c.plotiqu(150)


