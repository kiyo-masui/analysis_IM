#Simple script to divide Parkes beam maps by their estimated bandpass.
#For first case, will construct noise weighted average of 1st SVD modes.

import core.algebra as al
import numpy as np
import h5py
from itertools import combinations

input_map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_to_share/'

xx_cleaned_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_to_share/only_good_beams/XX/ra216/'

yy_cleaned_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_to_share/only_good_beams/YY/ra216/'

output_map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/beams_removed/'

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
    return weight_bp

beam_effic = [1.36, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.72, 1.72, 1.72, 1.72, 1.72, 1.72]

def normalize_map(bp, map, beam):
    #Divides out bandpass from map.  
    #Readjusts scale so that the average map deviation stays the same.
    #Then, divides by efficiency factor to convert from Jansky to Kelvin.
    #Efficiency Factors from http://www.atnf.csiro.au/research/multibeam/instrument/description.html
    map = np.transpose(map)
    map_dev = map - np.mean(map)
    scale_orig = np.sum(np.absolute(map_dev))
    map /= bp
    map_dev = map - np.mean(map)
    factor = scale_orig/np.sum(np.absolute(map_dev))
    map *= factor
    map = np.transpose(map)
    #Now, divide by efficiency factor
    #map = map/beam_effic[beam]
    return map

map_list = []

x_beams = [1,2,3,4,5,6,7,8,9,11,12,13]
y_beams = [1,2,3,4,5,6,7,8,9,10,12]

#for num in range(1,13):
for num in y_beams:
    map_list.append('parkes_parallel_thread_ra216dec0_p08_440by136_beam%d_clean_map_YY_1316.npy' % num)

i = 0
for map in map_list:
    map = al.load(input_map_dir + map)
    map = normalize_map(calculate_bandpass(y_beams[i], yy_cleaned_dir, y_beams, 15), map, i) 
    al.save(output_map_dir + 'parkes_parallel_thread_ra216dec0_p08_440by136_beam%d_clean_map_bp_div_YY_1316.npy' % y_beams[i], map)
    i +=1
