#Script for making noise weighted averages from different beam maps.
#Will combine the 12 or 11 beam maps in each polarization into 4 maps for each polariation.

import numpy as np
import core.algebra as al
import os
import copy

xx_beams = [[1,2,3],[4,5,6],[7,8,9],[11,12,13]]
yy_beams = [[1,2,3],[4,5,6],[7,8,9],[10,12]]

map_ra = os.getenv('MAP_RA')
map_dec = os.getenv('MAP_DEC')
map_size = os.getenv('MAP_SIZE')

base_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/renorm_hitconv_sync07/'
out_dir = base_dir + 'combined/'

def weighted_av_map(maps, weights):
    #Maps is a list of algebra vector maps.  Weights is a list of their weights.  Must be same size, in correct order.
    #Returns algebra vector map that is weighted average of maps, weighted by weights.
    #Maps must be same size and cover same region of sky, same freq coverage.
    map = np.multiply(maps[0], weights[0])
    weight = weights[0]
    hitweight = copy.deepcopy(weights[0])
    hitweight[map == 0.0] = 0.0
    for maps in zip(maps[1:],weights[1:]):
        map += np.multiply(maps[0],maps[1])
        weight += maps[1]
        hitweight[maps[0] != 0.0] += maps[1][maps[0] != 0]
    #In case there are zero weights, avoid dividing by zero.  
    bool = hitweight != 0.0
    map[bool] /= hitweight[bool]
    #Hitweight will only be zero where there were no hits on any beams.
    #Put the weight prior back at these points.
    bool = hitweight == 0.0
    hitweight[bool] = weight[bool]/len(weights)
    return [map, hitweight]

def get_maps(beams, dir, pol):
    maps = []
    for beam in beams:
        maps.append(al.load(dir + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam' + str(beam) + '_clean_map_bp_div_' + pol + '_1316.npy'))
    return maps

def get_weights(beams, dir, pol):
    maps = []
    for beam in beams:
        maps.append(al.load(dir + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam' + str(beam) + '_noise_inv_diag_bp_div_' + pol + '_1316.npy'))
    return maps

def save_final(data, beams, dir, pol):
    #data[0] is the map, data[1] is the weight.
    map_name = dir + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam' + ''.join(str(x) for x in beams) + '_clean_map_bp_div_' + pol + '_1316.npy'
    weight_name = dir + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam' + ''.join(str(x) for x in beams) + '_noise_inv_diag_bp_div_' + pol + '_1316.npy'
    al.save(map_name, data[0])
    al.save(weight_name, data[1])

if __name__ == '__main__':
   for beams in xx_beams:
       maps = get_maps(beams, base_dir, 'XX')
       weights = get_weights(beams, base_dir, 'XX')
       output = weighted_av_map(maps, weights)
       save_final(output, beams, out_dir, 'XX')
   for beams in yy_beams:
       maps = get_maps(beams, base_dir, 'YY')
       weights = get_weights(beams, base_dir, 'YY')
       output = weighted_av_map(maps, weights)
       save_final(output, beams, out_dir, 'YY')
