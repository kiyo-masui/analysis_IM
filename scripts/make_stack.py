# Simple script for stacking on galaxy locations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import core.algebra as al
from random import randint
import pylab

def gal_stack(gal_map, hi_map, hi_weight, pix_depth = 7):
    f_chans = gal_map.shape[0]
    stack = np.zeros(pix_depth,)
    null_stack = np.zeros(pix_depth,)
    weight = np.zeros(pix_depth,)
    null_weight = np.zeros(pix_depth,)
    not0 = np.argwhere(gal_map)
    ex = pix_depth - 1
    bool = np.logical_and(not0[:,0]>= ex/2, not0[:,0] <= f_chans - 1 - ex/2)
    not0_cut2 = not0[bool]
    w_map = hi_map * hi_weight
    #print not0_cut2.shape
    #not0_cut1 = not0[not0[:,0] != 0]
    #not0_cut2 = not0_cut1[not0_cut1[:,0] != 63]
    for i in range(not0_cut2.shape[0]):
        ind = not0_cut2[i,:]
        stack += w_map[ind[0]-ex/2:ind[0]+ex/2 + 1 ,ind[1], ind[2]]
        weight += hi_weight[ind[0]-ex/2:ind[0]+ex/2 + 1 ,ind[1], ind[2]]
        #print ind[0]
        map_slice = gal_map[ind[0],:,:]
        map_bool = map_slice == 0
        no_gals = np.argwhere(map_bool)
        rand_ind = randint(0, no_gals.shape[0] - 1)
        null_stack += w_map[ind[0]-ex/2:ind[0]+ex/2 + 1, no_gals[rand_ind][0], no_gals[rand_ind][1]]
        #print no_gals[0][rand_ind]
        #print rand_ind
        null_weight += hi_weight[ind[0]-ex/2:ind[0]+ex/2 + 1, no_gals[rand_ind][0], no_gals[rand_ind][1]]
    #stack /= not0_cut2.shape[0]
    null_stack /= null_weight
    stack /= weight
    return stack, weight, not0_cut2.shape[0], null_stack

def plot_stack(stack, output):
    len = np.shape(stack)[0]
    freqs = np.arange(len)- len/2 
    pylab.plot(freqs, stack)
    pylab.ylabel('S/N Estimate')
    pylab.xlabel('Offset in MHz')
    pylab.savefig(output + '.png')
    pylab.clf()

def plot_stacks(stacks, null_stack, output):
    len = np.shape(stack)[0]
    freqs = np.arange(len)- len/2
    pylab.plot(freqs, stack, 'b', label='2dF stack')
    pylab.plot(freqs, null_stack, 'r', label='Random no galaxy stack')
    pylab.ylabel('S/N Estimate')
    pylab.xlabel('Offset in MHz')
    pylab.legend(loc='upper right',prop={'size':8})
    pylab.savefig(output + '.png')
    pylab.clf()


if __name__ == '__main__':
    gal_dir = '/scratch2/p/pen/andersoc/2df_data/catalog_maps/'
    hi_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/hitconv_sync27_mbcal/I/'
    stack_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/plots/stacks/'
    hi_ras = ['ra165/','ra199/','ra33/','ran18/']
    gal_ras = ['pks_p1650_2dfmock_full/', 'pks_p1990_2dfmock_full/', 'pks_p3300_2dfmock_full/', 'pks_n1800_2dfmock_full/']
    for i in range(4):
        gal_map = al.load(gal_dir + gal_ras[i] + 'real_map_2df.npy')
        hi_map = al.load(hi_dir + hi_ras[i] + 'combined_clean_map_0modes.npy')
        hi_weight = al.load(hi_dir + hi_ras[i] + 'combined_clean_weight_0modes.npy')
        #Cut edges of the band
        #gal_map = gal_map[10:58]
        #hi_map = hi_map[10:58]
        #hi_weight = hi_weight[10:58]
        [stack, inv_noise, n_sources, null_stack] = gal_stack(gal_map, hi_map, hi_weight)
        #n_in_stack = gal_stack(gal_map, hi_map)[1]
        print hi_ras[i]
        #print np.mean(hi_map)
        #print np.std(hi_map)
        #print n_sources
        std = inv_noise**-0.5
        #print std
        #print stack
        stack /= (np.std(hi_map)*(n_sources**-0.5))
        null_stack /= (np.std(hi_map)*(n_sources**-0.5))
        #plot_stack(stack, stack_dir + 'clean10modes/' + '2df_stack_' + hi_ras[i][0:-2])
        #plot_stack(null_stack, stack_dir + 'clean10modes/' + 'nogal_stack_' + hi_ras[i][0:-2])
        plot_stacks(stack, null_stack, stack_dir + 'clean0modes/' + 'stack_' + hi_ras[i][0:-1])
        #print stack/std
