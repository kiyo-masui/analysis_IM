#Simple script to divide Parkes beam maps by their estimated bandpass.
#For first case, will construct noise weighted average of 1st SVD modes.

import core.algebra as al
import numpy as np
import h5py
from itertools import combinations
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

params_init = {'input_map_dir': '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_to_share/',
               'input_cleaned_root': '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_to_share/hitmap_fix/',
               'output_dir': '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/hitconv_sync27/'}

map_ra = os.getenv('MAP_RA')
map_dec = os.getenv('MAP_DEC')
map_size = os.getenv('MAP_SIZE')

plot_root = '/scratch2/p/pen/andersoc/second_parkes_pipe/bp_estimates/svd_estimates/'

bp_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/rebinned/bandpasses/'

#input_map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_to_share/'

#xx_cleaned_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_to_share/only_good_beams/XX/ra216/'

#yy_cleaned_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_to_share/only_good_beams/YY/ra216/'

#output_map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/beams_removed/'

effic = {'1': 1.36,
         '8': 1.45,
         '13': 1.72}

def calculate_bandpass(beam, cleaned_dir, cleaned_with, modes):
    #beam : integer, the beam to calculate the bandpass for
    #cleaned_dir : directory with SVD output.  Assuming hdf5 file called SVD.h5
    #cleaned_with : list of beams to use for bandpass estimation
    svd_data = h5py.File(cleaned_dir + 'SVD.hd5', 'r')
    weight_bp = np.zeros((64))
    weight_sum = 0
    for el in cleaned_with:
        if el != beam:
            el_smaller = el<beam
            noise_beam = al.load(cleaned_dir+'sec_beam' + str(beam) + '_' + 'cleaned_noise_inv' + '_I_with_beam' + str(el) + '_' + str(modes) + 'modes.npy')
            noise_el = al.load(cleaned_dir+'sec_beam' + str(el) + '_' + 'cleaned_noise_inv' + '_I_with_beam' + str(beam) + '_' + str(modes) + 'modes.npy')
            modes_beam = al.load(cleaned_dir+'sec_beam' + str(beam) + '_' + 'modes_clean_map' + '_I_with_beam' + str(el) + '_' + '1' + 'modes.npy').flatten()
            # Calculating sign of SVD  mode 1.
            # Since the main contribution to mode 1 should be strong point
            # sources, the sign of the bandpass will be the sign
            # of the pixel in the mode map with the largest amplitude.
            abs = np.absolute(modes_beam)
            sign = np.sign(modes_beam[np.argmax(abs)])
            weight = np.sqrt(np.sum(np.multiply(noise_beam,noise_el)))
            weight_sum += weight
            if el_smaller:
                weight_bp += sign*weight*svd_data['svd_modes2'][str(el)+ '_with_' +str(beam)][0]
            else:
                weight_bp += sign*weight*svd_data['svd_modes1'][str(beam)+ '_with_' +str(el)][0]
    weight_bp /= weight_sum
    print beam
    print weight_bp
    scale = 1.0/np.mean(weight_bp)
    print scale*weight_bp
    return scale*weight_bp

#beam_effic = [1.36, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.72, 1.72, 1.72, 1.72, 1.72, 1.72]

def sync_spectrum(spec_index = 2.7):
    sync = np.arange(1315.5-32,1315.5+32)/1315.5
    sync = sync**-spec_index
    return sync

def normalize_map(bp, map, beam, spec_index=2.7):
    #Divides out bandpass from map.  
    #Readjusts scale so that the average map deviation stays the same.
    #Then, divides by efficiency factor to convert from Jansky to Kelvin.
    #Efficiency Factors from http://www.atnf.csiro.au/research/multibeam/instrument/description.html
    map = np.transpose(map)
    #map_dev = map - np.mean(map)
    #scale_orig = np.sum(np.absolute(map_dev))
    bp /= sync_spectrum(spec_index)
    bp /= np.mean(bp)
    map /= bp
    #sync = np.arange(1315.5-32,1315.5+32)/1315.5
    #sync = sync**-spec_index
    #map *= sync
    #map_dev = map - np.mean(map)
    #factor = scale_orig/np.sum(np.absolute(map_dev))
    #map *= factor
    map = np.transpose(map)
    #Now, divide by efficiency factor
    #map = map/beam_effic[beam]
    return map

def normalize_noise(bp, noise, beam, spec_index=2.7):
    noise = np.transpose(noise)
    bp /= sync_spectrum(spec_index)
    bp /= np.mean(bp)
    noise *= bp**2
    #sync = np.arange(1315.5-32,1315.5+32)/1315.5
    #sync = sync**-spec_index
    #noise /= sync**2
    noise = np.transpose(noise)
    return noise

def plot_bp_estimates(bp_list, file_name):
    a = len(bp_list)
    print a
    fig, ax = plt.subplots(nrows=a, ncols=1, figsize=(12,18))
    for num in range(a):
        print num
        ax[num].plot(bp_list[num], 'b-')
    plt.savefig(file_name)
    plt.clf()
        
def save_bps(bp_array, beam_list, ra, dec, pol, dir):
    fname = 'bandpasses_' + str(beam_list) + '_' + ra + '_' + dec + '_' + pol
    fname = dir + fname
    np.save(fname, bp_array)

def apply_effic_factor(data, beam, type = 'bp'):
    #Data must be divided by efficiency factor to go from Jy to K.
    if beam == 1:
        fact = effic['1']
    elif beam <= 8:
        fact = effic['8']
    else:
        fact = effic['13']
    if type == 'map':
        data /= fact
    elif type == 'bp':
        data *= fact
    else:
        #data must be inverse noise weight.
        data *= fact**2
    return data

map_list = []
noise_list = []

x_beams = [1,2,3,4,5,6,7,8,9,11,12,13]
y_beams = [1,2,3,4,5,6,7,8,9,10,12]

#for num in range(1,13):

if __name__ == '__main__':
    bps = np.zeros((len(y_beams),64))
    sync = sync_spectrum()
    for num in y_beams:
        map_list.append('parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_clean_map_YY_1316.npy' % num)
        noise_list.append('parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_noise_inv_diag_YY_1316.npy' % num)
    i = 0
    for maps in zip(map_list, noise_list):
        map = al.load(params_init['input_map_dir'] + maps[0])
        noise = al.load(params_init['input_map_dir'] + maps[1])
        bp = calculate_bandpass(y_beams[i], params_init['input_cleaned_root'] + 'YY/' + map_ra + '/', y_beams, 2)
        bp /= sync
        bp /= np.mean(bp)
        bp = apply_effic_factor(bp, y_beams[i])
        bps[i,:] = bp
        #map = normalize_map(bp, map, i)
        #noise = normalize_noise(bp, noise, i) 
        #al.save(params_init['output_dir'] + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_clean_map_bp_div_YY_1316.npy' % y_beams[i], map)
        #al.save(params_init['output_dir'] + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_noise_inv_diag_bp_div_YY_1316.npy' % y_beams[i], noise)
        i +=1
    #plot_bp_estimates(bps, plot_root + map_ra + 'YY')
    save_bps(bps, y_beams, map_ra, map_dec, 'YY', bp_dir)

    map_list = []
    noise_list = []
    bps = np.zeros((len(x_beams),64))
    for num in x_beams:
        map_list.append('parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_clean_map_XX_1316.npy' % num)
        noise_list.append('parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_noise_inv_diag_XX_1316.npy' % num)

    i = 0
    for maps in zip(map_list, noise_list):
        map = al.load(params_init['input_map_dir'] + maps[0])
        noise = al.load(params_init['input_map_dir'] + maps[1])
        bp = calculate_bandpass(x_beams[i], params_init['input_cleaned_root'] + 'XX/' + map_ra + '/', x_beams, 2)
        bp /= sync
        bp /= np.mean(bp)
        bp = apply_effic_factor(bp, x_beams[i])
        bps[i,:] = bp
        #map = normalize_map(bp, map, i)
        #noise = normalize_noise(bp, noise, i)
        #al.save(params_init['output_dir'] + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_clean_map_bp_div_XX_1316.npy' % x_beams[i], map)
        #al.save(params_init['output_dir'] + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%d_noise_inv_diag_bp_div_XX_1316.npy' % x_beams[i], noise)
        i +=1
    #plot_bp_estimates(bps, plot_root + map_ra + 'XX')
    save_bps(bps, x_beams, map_ra, map_dec, 'XX', bp_dir)
