#! /usr/bin/env python
import sys
import numpy as np
from core import algebra
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from mkpower import functions
import math

# my first map
#maproot = "/mnt/raid-project/gmrt/ycli/map_result/maps/may04.2012/"
#mapname = "fir_15hr_41-90_clean_map_I_762"
#mapidex = 20

#maproot = "/mnt/raid-project/gmrt/ycli/map_result/maps/may04.2012/"
#mapname = "secA_15hr_41-90_clean_map_I_762"
#mapidex = 20

#maproot = "/mnt/raid-project/gmrt/ycli/map_result/maps/may10.2012/"
#mapname = "fir_1hr_41-90_clean_map_I_762"
#mapidex = 10

#maproot = "/mnt/raid-project/gmrt/ycli/map_result/maps/may10.2012/"
#mapname = "secD_1hr_41-90_dirty_map_I_762"
#mapidex = 10

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_15hr_map_oldcal_allsvd_mapmode_map/"
#mapname = "combined_clean_map_15modes"
#mapidex = 10

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_15hr_map_oldcal_legendre_modes_5gwj_mapmode_map/"
#mapname = "combined_clean_map_15modes"
#mapidex = 10

#maproot = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/"
#mapname = "sec_A_15hr_41-90_clean_map_I"
#mapidex = 20

#maproot = "/mnt/raid-project/gmrt/tcv/maps/"
#mapname = "15hr_41-90_fdgp_RM_clean_map_V"
#mapidex = 20

#maproot = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/"
#mapname = "sec_A_15hr_41-90_noise_weight_I"
#mapidex = 20

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_1hr_map_oldcal_legendre_modes_0gwj/mapmode_map/"
#mapname = "combined_clean_map_30modes"
#mapidex = 10

#maproot = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal/"
#mapname = "secA_1hr_41-18_clean_map_I_800"
#mapidex = 10

#maproot = "/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_1hr_map_oldcal/"
#mapname = "combined_clean_map_50modes"
##mapname = "sec_A_cleaned_clean_map_I_with_B_50modes"
#mapidex = 10

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_1hr_map_oldcal_legendre_modes_0gwj/mapmode_map/"
#mapname = "combined_clean_map_50modes"
#mapname = "sec_A_cleaned_clean_map_I_with_B_50modes"
#mapidex = 10

# maps and fftboxes
#maproot = "/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_oldmap_ideal/"
#mapname = "sim_temperature_000"
#mapidex = 200
#boxroot = "/mnt/raid-project/gmrt/ycli/ps_result/simulation_auto_sim_15hr_oldmap_ideal_15/fftbox/"
#boxname = "fftbox_sim_temperature_000"

#maproot = "/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_oldmap_str/"
#mapname = "sim_temperature_000"
#mapidex = 250
#boxroot = "/mnt/raid-project/gmrt/ycli/ps_result/simulation_auto_sim_15hr_oldmap_str_15/fftbox/"
#boxname = "fftbox_sim_temperature_000"

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQV_legendre_modes_0gwj/IxI5svd/"
#mapname = "sec_I_cleaned_clean_map_I_with_I_3modes"
#mapidex = 150

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/IQV_legendre_modes_0gwj/IQUmap_clean_withIxIsvd/"
#mapname = "sec_I_cleaned_clean_map_I_with_Q_1modes"
##mapname = "sec_Q_cleaned_clean_map_I_with_I_0modes"
##mapname = "sec_U_cleaned_clean_map_I_with_I_1modes"
#mapidex = 150

#maproot = "/mnt/data-pen3/ycli/map_result/maps/parkes_linear_highres/"
#mapname = "fir_dirty_map_I_1315"
#mapname = "fir_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12_west_dirty_map_I_1315"
#mapname = "fir_parkes_2008_09_12_west_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_13_west_dirty_map_I_1315"
#mapname = "fir_parkes_2008_09_13_west_clean_map_I_1315"
#mapidex = 150

#maproot = "/mnt/data-pen3/ycli/map_result/maps/parkes_12-13-14/"
maproot = "/mnt/raid-project/gmrt/ycli/map_result/maps/parkes/"
#maproot = "/mnt/data-pen3/ycli/map_result/maps/parkes_linear_highres_combine/"
#mapname = "fir_parkes_2008_09_12-13_west_dirty_map_I_1315"
#mapname = "fir_parkes_2008_09_12-13_west_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12-13-14_west_dirty_map_I_1315"
#mapname = "fir_parkes_2008_09_11-12-13-14_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w-13-14_beam0testclean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w-13-14_beam0testdirty_map_I_1315"

#mapname = "fir_parkes_2008_09_12w_13_14dirty_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_13_14_beam_0_dirty_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_13_14_beam_1_dirty_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_13_14_beam_2_dirty_map_I_1315"

#mapname = "fir_parkes_2008_09_12w_13_14_beam_0_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_13_14_beam_1_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_13_14_beam_2_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_13_14_beam_3_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_13_14_clean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_beam0testclean_map_I_1315"
#mapname = "fir_parkes_2008_09_12w_beam1testdirty_map_I_1315"

#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_dirty_map_I_1315"
#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_beam_0_dirty_map_I_1315"
#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_beam_1_dirty_map_I_1315"
#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_beam_2_dirty_map_I_1315"
#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_clean_map_I_1315"

#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_clean_map_I_1315"
#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_beam_0_clean_map_I_1315"
#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_beam_1_clean_map_I_1315"
#mapname = "fir_RA+10_parkes_2008_09_12w_13_14_beam_2_clean_map_I_1315"

mapname = "fir_RA+05_parkes_2008_09_12w_13_14_clean_map_I_1315"

#maproot = "/mnt/data-pen3/tcv/oldmaps/1hr_41-16_fdg/"
#mapname = "secA_1hr_41-18_clean_map_I_800"
#maproot = "/mnt/raid-project/gmrt/tcv/maps/"
#mapname = "fir_1hr_41-18_avg_fdgp_clean_map_I_800"
#mapname = "fir_1hr_41-18_avg_fdgp_clean_map_Q_800"
#mapname = "fir_1hr_41-18_avg_fdgp_clean_map_U_800"
#mapname = "fir_1hr_41-18_avg_fdgp_clean_map_V_800"
#maproot = "/mnt/data-pen3/ycli/map_gbt/15hr_IQUV_extend_legendre_modes_0gwj_conv/"
#mapname = "15hr_41-80_avg_fdgp_noise_weight_I"
#mapname = "15hr_41-80_avg_fdgp_noise_weight_Q"
#mapname = "15hr_41-80_avg_fdgp_noise_weight_V"
#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_AQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/"
#mapname = "sec_A_cleaned_noise_inv_I_with_B_0modes"
#mapname = "sec_B_cleaned_noise_inv_I_with_A_0modes"
#mapname = "sec_A_cleaned_clean_map_I_with_B_0modes"
#mapname = "sec_B_cleaned_clean_map_I_with_A_0modes"
#mapname = "fir_1hr_41-18_avg_fdgp_noise_diag_I_800"
#maproot = "/mnt/data-pen3/ycli/map_gbt/1hr_AQUV_extend_legendre_modes_0gwj_conv/"
#mapname = "secA_1hr_41-18_noise_weight_I_800"
#mapname = "secA_1hr_41-18_noise_weight_Q_800"
#mapname = "secB_1hr_41-18_noise_weight_U_800"
#mapname = "secB_1hr_41-18_noise_weight_V_800"
#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/1hr_AQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_20modes"
#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/AQUV_extend_legendre_modes_0gwj_conv/Emap_clean_themselves/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_20modes"

mapidex = 32

def getedge(map):
    deg2rad = 3.1415926/180.
    ra   = map.get_axis('ra')*deg2rad
    ra   = ra - ra[round(ra.shape[0]/2.)]
    dec  = map.get_axis('dec')*deg2rad
    freq = map.get_axis('freq')
    r    = functions.fq2r(freq)

    print 'freq range: %f ~ %f'%(freq[0], freq[-1])
    print 'r    range: %f ~ %f'%(r[0], r[-1])
    print r.max() - r.min()

    r.shape   = r.shape + (1, 1)
    ra.shape  = (1,) + ra.shape + (1,)
    dec.shape = (1, 1) + dec.shape
    radec_map = np.zeros((3,)+map.shape)
    radec_map[0,...] = r
    radec_map[1,...] = ra
    radec_map[2,...] = dec
    xyz_map = np.zeros((3,)+map.shape)
    xyz_map[0,...]=radec_map[0,...]*np.cos(radec_map[2,...])*np.cos(radec_map[1,...])
    xyz_map[1,...]=radec_map[0,...]*np.cos(radec_map[2,...])*np.sin(radec_map[1,...])
    xyz_map[2,...]=radec_map[0,...]*np.sin(radec_map[2,...])

    x = xyz_map[0].flatten()
    y = xyz_map[1].flatten()
    z = xyz_map[2].flatten()

    xrange = (math.floor(x.min()), math.floor(x.max())+1.)
    yrange = (math.floor(y.min()), math.floor(y.max())+1.)
    zrange = (math.floor(z.min()), math.floor(z.max())+1.)

    plt.figure(figsize=(10,10))
    plt.subplot(311)
    plt.scatter(x,y, marker='.', alpha=0.1)
    plt.vlines(xrange[0], yrange[0], yrange[1])
    plt.vlines(xrange[1], yrange[0], yrange[1])
    plt.hlines(yrange[0], xrange[0], xrange[1])
    plt.hlines(yrange[1], xrange[0], xrange[1])
    plt.subplot(312)
    plt.scatter(y,z, marker='.', alpha=0.1)
    plt.vlines(yrange[0], zrange[0], zrange[1])
    plt.vlines(yrange[1], zrange[0], zrange[1])
    plt.hlines(zrange[0], yrange[0], yrange[1])
    plt.hlines(zrange[1], yrange[0], yrange[1])
    plt.subplot(313)
    plt.scatter(x,z, marker='.', alpha=0.1)
    plt.vlines(xrange[0], zrange[0], zrange[1])
    plt.vlines(xrange[1], zrange[0], zrange[1])
    plt.hlines(zrange[0], xrange[0], xrange[1])
    plt.hlines(zrange[1], xrange[0], xrange[1])

    plt.savefig('./png/mapsize.png', formate='png')

map = algebra.load(maproot + mapname + '.npy')
map = algebra.make_vect(map)
#getedge(map)
#exit()
print map.shape
#print map.get_axis('freq')
#exit()

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
    badfreq = [0, 9, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23,]
    if len(badfreq) != 0:
        mask = np.zeros(map.shape)
        mask[badfreq,:,:] = 1.
        map = np.ma.masked_array(map, mask=mask)

    #vmax = map[mapidex].flatten().max()
    #vmin = map[mapidex].flatten().min()
    #vmax = map.flatten().max()
    #vmin = map.flatten().min()
    #vmax = 0.05
    #vmin = -0.05

    vmax = 0.04
    vmin = -0.04

    #vmax = 10
    #vmin = -10
    #if vmax<np.fabs(vmin):
    #    vmax = np.fabs(vmin)
    #else:
    #    vmin = -np.fabs(vmax)

    map = np.ma.masked_equal(map, 0)
    #vmax = 1.
    #vmin = -1.

    #plt.figure(figsize=(34, 8))
    #plt.figure(figsize=(9,8))
    #f = plt.figure(figsize=(16,9))
    f = plt.figure(figsize=(10,16))
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
    badfreq = [0, 9, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23,]
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
    boxidex = int((r[mapidex]-1400.)/2)
    
    box = np.load(boxroot + boxname + '.npy')
    
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
#plt.show()



