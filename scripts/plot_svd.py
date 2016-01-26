import core.algebra as al
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
from itertools import combinations

def make_svd_matrix(pairs, h5_path, mode, file_name, num=0):
    results = h5py.File(h5_path + 'SVD.hd5','r')
    fig, ax = plt.subplots(nrows=13, ncols=13, figsize=(12,18))
    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, top=0.97, wspace=0.01, hspace=0.01)
    i=0
    for el in pairs:
        current = str(el[0]+1) + '_with_' + str(el[1]+1)
        vals = results['svd_vals'][current]
        vals1 = results['svd_modes1'][current]
        vals2 = results['svd_modes2'][current]
        #weights = results['cov_counts'][current]
        if mode == 'svd_vals':
            #ax[el[1]][el[0]].plot(np.log10(total_map_svd[i][0]/total_map_svd[i][0][0]), 'b-')
            ax[el[1]][el[0]].plot(np.log10(vals[0:6]/vals[0]), 'b-')
        #if mode == 'cov_counts': 
        if mode == 'vect':
            print vals1[num]
            print vals2[num]
            ax[el[1]][el[0]].plot(vals1[num], 'b-')
            ax[el[0]][el[1]].plot(vals2[num], 'b-')
            if el[0] != 0:
                ax[el[1]][el[0]].yaxis.set_ticklabels([])
            ax[el[1]][el[0]].xaxis.set_ticklabels([])
            if el[1] != 13:
                ax[el[0]][el[1]].xaxis.set_ticklabels([])
            ax[el[0]][el[1]].yaxis.set_ticklabels([])
            i += 1
    plt.savefig(file_name)
    plt.clf()

def make_map_matrix(pairs, base_dir, type, modes, map_op, file_name):
    fig, ax = plt.subplots(nrows=13, ncols=13, figsize=(12,18))
    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.05, top=0.97, wspace=0.01, hspace=0.01)
    for el in pairs:
        current = str(el[0]+1) + '_with_' + str(el[1]+1)
        map1 = al.load(base_dir+'sec_beam' + str(el[0]+1) + '_' + type + '_I_with_beam' + str(el[1]+1) + '_' + str(modes) + 'modes.npy')
        map2 = al.load(base_dir+'sec_beam' + str(el[1]+1) + '_' + type + '_I_with_beam' + str(el[0]+1) + '_' + str(modes) + 'modes.npy')
        noise1 = al.load(base_dir+'sec_beam' + str(el[0]+1) + '_' + 'cleaned_noise_inv' + '_I_with_beam' + str(el[1]+1) + '_' + str(modes) + 'modes.npy')
        noise2 = al.load(base_dir+'sec_beam' + str(el[1]+1) + '_' + 'cleaned_noise_inv' + '_I_with_beam' + str(el[0]+1) + '_' + str(modes) + 'modes.npy')
        results = h5py.File(base_dir + 'SVD.hd5','r')
        vals1 = results['svd_modes1'][current]
        vals2 = results['svd_modes2'][current]
        modes1 = al.load(base_dir+'sec_beam' + str(el[0]+1) + '_' + 'modes_clean_map' + '_I_with_beam' + str(el[1]+1) + '_' + '1' + 'modes.npy')
        modes2 = al.load(base_dir+'sec_beam' + str(el[1]+1) + '_' + 'modes_clean_map' + '_I_with_beam' + str(el[0]+1) + '_' + '1' + 'modes.npy')
        if divide == True:
            map1 = np.transpose(map1)
            map2 = np.transpose(map2)
            if abs(np.max(modes1))<abs(np.min(modes1)):
                normal1 = -vals1[0]
            if abs(np.max(modes2))<abs(np.min(modes2)):
                normal2 = -vals2[0]
            map1 /= normal1
            map2 /= normal2
            map1 = np.transpose(map1)
            map2 = np.transpose(map2)
        #map1 = map_oper(map1, noise1, map_op)
        #map2 = map_oper(map2, noise2, map_op)
        map1 = mixed_map_oper(map1, noise1, map2, noise2, map_op)
        map2 = map1
        #print map1.shape
        #print map2.shape
        ax[el[1]][el[0]].plot(map1, 'b-')
        ax[el[0]][el[1]].plot(map2, 'b-')
    plt.savefig(file_name)
    plt.clf()

def map_oper(map, noise, map_op):
    if map_op == 'freq_fluc_size':
        noise = noise.reshape((map.shape[0],map.shape[1]*map.shape[2]))
        map = map.reshape((map.shape[0],map.shape[1]*map.shape[2]))
        n_av = np.divide(np.sum(np.multiply(noise, map), 1), np.sum(noise,1))
        map = np.transpose(map)
        noise = np.transpose(noise)
        map -= n_av
        power = np.sum(np.multiply(np.square(noise),np.square(map)), 0)
        power = np.sqrt(power)/np.sum(noise,0)
        #abs = np.std(map,1)
        return power

def mixed_map_oper(map1, noise1, map2, noise2, map_op):
    if map_op == 'freq_fluc_size':
        noise1 = noise1.reshape((map1.shape[0],map1.shape[1]*map1.shape[2]))
        map1 = map1.reshape((map1.shape[0],map1.shape[1]*map1.shape[2]))
        noise2 = noise2.reshape((map2.shape[0],map2.shape[1]*map2.shape[2]))
        map2 = map2.reshape((map2.shape[0],map2.shape[1]*map2.shape[2]))
        n_av1 = np.divide(np.sum(np.multiply(noise1, map1), 1), np.sum(noise1,1))
        n_av2 = np.divide(np.sum(np.multiply(noise2, map2), 1), np.sum(noise2,1))
        map1 = np.transpose(map1)
        noise1 = np.transpose(noise1)
        map2 = np.transpose(map2)
        noise2 = np.transpose(noise2)
        map1 -= n_av1
        map2 -= n_av2
        power = np.sum(np.multiply(np.multiply(noise1,noise2),np.multiply(map1,map2)), 0)
        power = power/np.sum(np.multiply(noise1,noise2),0)
        #abs = np.std(map,1)
        return power
    
        

comb_num = [comb for comb in combinations(range(13),2)]
base_dir_clip = '/scratch/p/pen/andersoc/combined_maps/I_cleaned/freq_flag_clip_med/'
#base_dir = '/scratch/p/pen/andersoc/combined_maps/I_cleaned/'
#base_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_to_share/XX/'

base_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/XX/'

ra = 'ra50/'
output_dir = '/scratch2/p/pen/andersoc/combined_maps/svd_plots/'
#make_svd_matrix(pairs=comb_num, h5_path = base_dir_clip + ra, mode='vect', file_name = output_dir + 'evects_0_clipped_ra50', num=0)

pol = '_XX'

#ra_list = ['ra10', 'ra20', 'ra30', 'ra40', 'ra50']
ra_list = ['ra33']
print comb_num
divide = False
for val in ra_list:
    #make_svd_matrix(pairs=comb_num, h5_path = base_dir + val + '/', mode='vect', file_name = output_dir + 'evects_1_' + val + pol, num=1)
    #make_svd_matrix(pairs=comb_num, h5_path = base_dir + val + '/', mode='svd_vals', file_name = output_dir + 'bp_div_svd_vals_first6_' + val + pol, num=0)
    make_map_matrix(pairs=comb_num, base_dir = base_dir + val + '/', type = 'cleaned_clean_map', modes = 15, map_op = 'freq_fluc_size', file_name = output_dir + 'noise_weighted_freq_fluc_cross_power_15modes_rem_bp_div_thenclean')
