#Quick script to check relative amplitudes of two normalizations

import numpy as np
import core.algebra as al

bp_folder = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/beams_removed/'

map_prefix = 'parkes_parallel_thread_'

radec = 'ra216dec0_p08_440by136'

suffix = '_clean_map_bp_div_'

pol = 'XX'

end = '_1316.npy'

x_beams = [1,2,3,4,5,6,7,8,9,11,12,13]
y_beams = [1,2,3,4,5,6,7,8,9,10,12]

def unnorm_map(beam, pol, region):
    string = bp_folder + map_prefix + region + '_beam' + str(beam) + suffix + pol + end
    return al.load(string)

def norm_map(beam, pol, region):
    string = bp_folder + 'renormalized/' + map_prefix + region + '_beam' + str(beam) + suffix + pol + end
    return al.load(string)  

def calc_factor(beam, pol, region):
    norm = norm_map(beam, pol, region)
    unnorm = unnorm_map(beam, pol, region)
    div = unnorm/norm
    vals = div[~np.isnan(div)]
    return [vals.mean(), vals.min(), vals.max()]


radec = 'ran18decn30_p08_488by106'
radecs = ['ran18decn30_p08_488by106', 'ra33decn30_p08_568by106', 'ra165dec0_p08_440by136', 'ra182dec0_p08_440by136', 'ra199dec0_p08_440by136', 'ra216dec0_p08_440by136']
for reg in radecs:
    print reg
    pol = 'XX'
    print '\t' + pol
    for beam in x_beams:
        print '\t\t' + 'beam' + str(beam) + ' ' + str(calc_factor(beam, pol, reg))
    pol = 'YY'
    print '\t' + pol
    for beam in y_beams:
        print '\t\t' + 'beam' + str(beam) + ' ' + str(calc_factor(beam, pol, reg))

