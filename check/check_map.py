#! /usr/bin/env python
import sys
import numpy as np
from core import algebra
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from mkpower import functions
import math

maproot = "/Users/ycli/DATA/2df/map_2929.5/"
mapname = "real_map_2df"
boxroot = "/Users/ycli/DATA/2df/map_2929.5/"
#boxname = "fftbox_real_map_2df"
boxname = "kiyo_fftbox_real_map_2df"

#maproot = "/Users/ycli/DATA/maps/"
#mapname = "secA_15hr_41-80_avg_fdgp_new_clean_map_I_800"
#boxroot = "/Users/ycli/DATA/maps/"
#boxname = "fftbox_secA_15hr_41-80_avg_fdgp_new_clean_map_I_800"

mapidex = 50

map = algebra.make_vect(algebra.load(maproot + mapname + '.npy'))

if mapidex>=map.shape[0]:
    print 'index out of the maps!!'
    exit()

if len(sys.argv) == 1:
    freq = map.get_axis('freq') - 0.5*map.info['freq_delta']
    ra = map.get_axis('ra') - 0.5*map.info['ra_delta']
    dec= map.get_axis('dec') - 0.5*map.info['dec_delta']

    freq = np.resize(freq, freq.shape[0]+1)
    ra = np.resize(ra, ra.shape[0]+1)
    dec = np.resize(dec, dec.shape[0]+1)

    freq[-1] = freq[-2] + map.info['freq_delta']
    ra[-1] = ra[-2] + map.info['ra_delta']
    dec[-1] = dec[-2] + map.info['dec_delta']

    badfreq = []
    #badfreq = [0, 9, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23,]
    badfreq = range(10)
    if len(badfreq) != 0:
        mask = np.zeros(map.shape)
        mask[badfreq,:,:] = 1.
        map = np.ma.masked_array(map, mask=mask)

    #vmax = map[mapidex].flatten().max()
    #vmin = map[mapidex].flatten().min()
    #vmax = map.flatten().max()
    #vmin = map.flatten().min()
    #vmax = 0.5
    #vmin = -0.5

    vmax = 2
    vmin = -2
    #vmax = 0.001
    #vmin = -0.001
    #vmax = 0.0003
    #vmin = -0.0003

    #vmax = 50
    #vmin = -50

    #if vmax<np.fabs(vmin):
    #    vmax = np.fabs(vmin)
    #else:
    #    vmin = -np.fabs(vmax)

    #map = np.ma.masked_equal(map, 0)
    map = np.ma.array(map)
    map[map==0] = np.ma.masked
    #vmax = 1.
    #vmin = -1.

    #plt.figure(figsize=(34, 8))
    #plt.figure(figsize=(9,8))
    #f = plt.figure(figsize=(16,9))
    #f = plt.figure(figsize=(10,16))
    f = plt.figure(figsize=(8,7))
    ax = ImageGrid(f, 111,
                   nrows_ncols = (1, 1),
                   direction = "row",
                   axes_pad = 0.05,
                   add_all = True,
                   label_mode = "L",
                   share_all = True,
                   cbar_location = "right",
                   cbar_mode = "single",
                   cbar_size = 0.2,
                   cbar_pad = 0.05,
                   )
    im = ax[0].pcolormesh(ra, dec, np.ma.mean(map, 0).swapaxes(0,1))
    #im = ax[0].pcolormesh(ra, dec, map[mapidex].swapaxes(0,1))
    im.set_clim(vmin, vmax)
    ax[0].set_xlim(ra.min(), ra.max())
    ax[0].set_ylim(dec.min(), dec.max())
    ax[0].set_xlabel('RA [deg]')
    ax[0].set_ylabel('DEC[deg]')
    ax[0].set_title(mapname)

    ax[0].cax.colorbar(im)

    plt.tick_params(length=6, width=1.)
    plt.tick_params(which='minor', length=3, width=1.)

    #plt.pcolor(ra, dec, map[mapidex].swapaxes(0,1), vmax=vmax, vmin=vmin)
    #plt.pcolor(ra, dec, np.mean(map, 0).swapaxes(0,1), vmax=vmax, vmin=vmin)

elif len(sys.argv) == 2 and sys.argv[1]=='freq':
    mapname += '_freq'

    freq = map.get_axis('freq') - 0.5*map.info['freq_delta']
    ra = map.get_axis('ra') - 0.5*map.info['ra_delta']
    dec= map.get_axis('dec') - 0.5*map.info['dec_delta']

    freq = np.resize(freq, freq.shape[0]+1)
    ra = np.resize(ra, ra.shape[0]+1)
    dec = np.resize(dec, dec.shape[0]+1)

    freq[-1] = freq[-2] + map.info['freq_delta']
    ra[-1] = ra[-2] + map.info['ra_delta']
    dec[-1] = dec[-2] + map.info['dec_delta']

    badfreq = []
    #badfreq = [0, 9, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23,]
    badfreq = range(10)
    if len(badfreq) != 0:
        mask = np.zeros(map.shape)
        mask[badfreq,:,:] = 1.
        map = np.ma.masked_array(map, mask=mask)

    #vmax = map.flatten().max()
    #vmin = map.flatten().min()
    vmax = 0.5
    vmin = -0.5

    #map = np.ma.masked_equal(map, 0)

    freq /= 1.e9
    print freq

    f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10,10))
    f.subplots_adjust(hspace=0)

    #pc1 = ax1.pcolor(freq, ra,  map.swapaxes(0,1).mean(axis=2))
    #ax1.set_ylim(ra.min(), ra.max())
    #ax1.set_ylabel('RA [deg]')
    #pc2 = ax2.pcolor(freq, dec, map.swapaxes(0,2).mean(axis=1))
    #ax2.set_ylim(dec.min(), dec.max())
    #ax2.set_ylabel('DEC[deg]')
    for i in range(ra.shape[0]-1):
        for j in range(dec.shape[0]-1):
            ax3.plot(freq[:-1]+0.5*map.info['freq_delta']/1.e9, map[:,i,j])
    
    pc1 = ax1.pcolor(freq, ra,  map.swapaxes(0,1)[:,:,30], vmin=vmin, vmax=vmax)
    ax1.set_ylim(ra.min(), ra.max())
    ax1.set_ylabel('RA [deg]')
    pc2 = ax2.pcolor(freq, dec, map.swapaxes(0,2)[:,30,:], vmin=vmin, vmax=vmax)
    ax2.set_ylim(dec.min(), dec.max())
    ax2.set_ylabel('DEC[deg]')
    #for i in range(ra.shape[0]-1):
    #    ax3.plot(freq[:-1]+0.5*map.info['freq_delta']/1.e9, map[:,i,30])
    #for i in range(dec.shape[0]-1):
    #    ax3.plot(freq[:-1]+0.5*map.info['freq_delta']/1.e9, map[:,30,i])
        
    ax3.set_ylabel('Flux [Jy]')
    ax3.set_xlabel('Freqency [GHz]')
    ax3.set_xlim(freq.min(), freq.max())

    #f.colorbar(pc1, ax1, orientation='horizontal')


elif len(sys.argv) == 2:
    freq = map.get_axis('freq')
    z = 1.42e9/freq - 1.
    r = functions.fq2r(freq)
    boxidex = int((r[mapidex]-83.7065)/0.597)
    #boxidex = int((r[mapidex]-1507.03)/1.768)
    
    box = algebra.make_vect(algebra.load(boxroot + boxname + '.npy'))
    r_box = box.get_axis('freq')
    print r_box
    print r[mapidex]
    boxidex = np.digitize([r[mapidex],], r_box)[0]
    #box = np.load(boxroot + boxname + '.npy')
    
    if sys.argv[1]=='1':
        plt.figure(figsize=(8, 7))
        plt.subplot(211)
        plt.imshow(map[mapidex].swapaxes(0,1), interpolation='nearest')
        plt.colorbar()
        plt.subplot(212)
        plt.imshow(box[boxidex].swapaxes(0,1), interpolation='nearest')
        plt.colorbar()
        print map[mapidex].flatten().max()
        print map[mapidex].flatten().min()
        print box[boxidex].flatten().max()
        print box[boxidex].flatten().min()
    
    elif sys.argv[1]=='2':
        plt.figure(figsize=(8, 10))
        plt.subplot(211)
        plt.imshow(map.swapaxes(0,1)[64].swapaxes(0,1), interpolation='nearest')
        plt.colorbar()
        plt.subplot(212)
        plt.imshow(box.swapaxes(0,1)[64].swapaxes(0,1), interpolation='nearest')
        plt.colorbar()
    elif sys.argv[1]=='3':
        #maphist = np.histogram(map.flatten(), bins=50)
        #boxhist = np.histogram(box.flatten(), bins=50)
        plt.figure(figsize=(8,8))
        #plt.hist(map.flatten(), range=(1.e-5, 0.003), bins=100, color='b',
        #         normed=True, alpha=0.5, label='map')
        #plt.hist(map.flatten(), range=(-0.003, -1.e-5), bins=100, color='b',
        #         normed=True, alpha=0.5)
        #plt.hist(box.flatten(), range=(1.e-5, 0.003), bins=100, color='g',
        #         normed=True, alpha=0.5, label='box')
        #plt.hist(box.flatten(), range=(-0.003, -1.e-5), bins=100, color='g',
        #         normed=True, alpha=0.5)
        #plt.ylim(ymax=5000)
    
        plt.hist(map.flatten(), bins=100, color='b', 
                 normed=True, alpha=0.5, label='map')
        plt.hist(box.flatten(), bins=100, color='g', 
                 normed=True, alpha=0.5, label='box')
        plt.ylim(ymax=1000)
        #plt.ylim(ymax=100000)
        plt.legend()



pngname = './png/'+mapname+'%03d.png'%mapidex
plt.savefig(pngname, format='png')
plt.show()



