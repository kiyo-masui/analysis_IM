#! /usr/bin/env python 

import numpy as np
import numpy.ma as ma
import scipy as sp
import matplotlib.pyplot as plt
from core import algebra

def make_factorizable(noise):
    r"""factorize the noise"""

    weight_prior = 2
    noise[noise < weight_prior] = 1.e-30
    noise = 1. / noise
    noise = ma.array(noise)
    # Get the freqency averaged noise per pixel.  Propagate mask in any
    # frequency to all frequencies.
    for noise_index in range(ma.shape(noise)[0]):
        if sp.all(noise[noise_index, ...] > 1.e20):
            noise[noise_index, ...] = ma.masked
    noise_fmean = ma.mean(noise, 0)
    noise_fmean[noise_fmean > 1.e20] = ma.masked
    # Get the pixel averaged noise in each frequency.
    noise[noise > 1.e20] = ma.masked
    noise /= noise_fmean
    noise_pmean = ma.mean(ma.mean(noise, 1), 1)
    # Combine.
    noise = noise_pmean[:, None, None] * noise_fmean[None, :, :]
    noise = (1. / noise).filled(0)

    return noise

def do_corr(maproot, mapname, wetname, mapname2=None, wetname2=None):

    def preparemap(maproot, mapname, wetname):
        map = algebra.load(maproot + mapname + '.npy')
        map = algebra.make_vect(map)
        
        wet = algebra.load(maproot + wetname + '.npy')
        wet = algebra.make_vect(wet)
        
        wet = make_factorizable(wet)
        #wet = algebra.as_alg_like(wet, map)
        
        mean = sp.sum(sp.sum(wet*map, -1), -1)
        mean /= sp.sum(sp.sum(wet, -1), -1)
        mean.shape += (1, 1)
        
        map -= mean
        
        map = map*wet
        
        map = ma.array(map)
        
        map = map.reshape(map.shape[0], -1)
        wet = wet.reshape(wet.shape[0], -1)
        
        return map, wet

    map, wet = preparemap(maproot, mapname, wetname)

    if mapname2 is None:
        map2 = map
        wet2 = wet
    else:
        map2, wet2 = preparemap(maproot, mapname2, wetname2)
    
    corr = np.dot(map.T, map2)
    wcorr = np.dot(wet.T, wet2)
    
    corr = corr/wcorr
    corr = np.ma.masked_where(np.isnan(corr), corr)
    #corr = wcorr
    #print corr[:10,:10]
    
    print corr.shape

    return corr


maproot = "/mnt/data-pen3/ycli/map_gbt/15hr_RM/"
mapname = "15hr_41-90_fdgp_RM_clean_map_Q"
wetname = "15hr_41-90_fdgp_RM_noise_weight_Q"
mapname2= "15hr_41-90_fdgp_RM_clean_map_U"
wetname2= "15hr_41-90_fdgp_RM_noise_weight_U"

corr = do_corr(maproot, mapname, wetname, mapname2, wetname2)

corr -= corr.T

plt.figure(figsize=(10,8))
plt.pcolor(corr[:200, :200])
plt.colorbar()
plt.savefig('./png/corr_' + mapname + 'x' +  mapname2 + '.png', formate='png')

mapname = "15hr_41-90_fdgp_RM_clean_map_U"
wetname = "15hr_41-90_fdgp_RM_noise_weight_U"
mapname2= "15hr_41-90_fdgp_RM_clean_map_Q"
wetname2= "15hr_41-90_fdgp_RM_noise_weight_Q"

corr = do_corr(maproot, mapname, wetname, mapname2, wetname2)

corr -= corr.T

plt.figure(figsize=(10,8))
plt.pcolor(corr[:200, :200])
plt.colorbar()
plt.savefig('./png/corr_' + mapname + 'x' +  mapname2 + '.png', formate='png')


