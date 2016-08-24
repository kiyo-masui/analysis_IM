import core.algebra as al
import numpy as np
import h5py
from itertools import combinations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from foreground_clean import map_pair as mp

input_map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_to_share/'

x_beams = [1,2,3,4,5,6,7,8,9,11,12,13]
y_beams = [1,2,3,4,5,6,7,8,9,10,12]

comm_weight = al.load(input_map_dir + 'parkes_parallel_thread_ra33decn30_p08_568by106_beam1_noise_inv_diag_XX_1316.npy')
print comm_weight.info
for num in x_beams:
    curr_weight = al.load(input_map_dir + 'parkes_parallel_thread_ra33decn30_p08_568by106_beam%d_noise_inv_diag_XX_1316.npy' % num)
    min = np.min(curr_weight)
    hit = (curr_weight > min)
#    print np.sum(curr_weight)
#    print np.sum(curr_weight[hit])
    comm_weight[...] = np.minimum(comm_weight,curr_weight)

print comm_weight.info

min = np.min(comm_weight)

no_hit = (comm_weight == min)

comm_weight[no_hit] = 0


file_start = 'parkes_parallel_thread_ra33decn30_p08_568by106_beam'
noise_end = '_noise_inv_diag_XX_1316.npy'
map_end = '_clean_map_XX_1316.npy'
#Now, make list of map_pairs, using the common weight.
#Doing pairs to use existing beam convolution, noise factoring,
#and mean subtraction routines
conv = 1.4
pair_list = []
init_beam = 1
freq = range(64)
#Filling list
for num in x_beams:
    if num != init_beam:
        map1 = al.load(input_map_dir + file_start + str(init_beam) + map_end)
        map2 = al.load(input_map_dir + file_start + str(num) + map_end)
        pair = mp.MapPair(map1, map2, comm_weight, comm_weight, freq, conv_factor=conv)
        pair.degrade_resolution()
        pair.make_noise_factorizable()
        pair.subtract_weighted_mean()
        pair_list.append(pair)
