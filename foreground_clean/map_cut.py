import os
import numpy as np
from core import algebra
import matplotlib.pyplot as plt

def map_cut(input_file, output_file, new_centre, new_shape):

    map = algebra.load(input_file)
    map = algebra.make_vect(map)
    
    ra = map.get_axis('ra')
    dec = map.get_axis('dec')
    
    centre_indx = [((map.info['ra_centre']  - ra )**2).argmin(),
                   ((map.info['dec_centre'] - dec)**2).argmin()]
    
    print ">>cuting map : \n  %s"%input_file
    print "  old map:",
    print "centre at RA=%5.4f DEC=%5.4f"%(map.info['ra_centre'], 
                                          map.info['dec_centre'])
    print "           shape ",
    print map.shape
    
    '''>>> adjust ra <<<'''
    if new_centre[0] != None:
        centre_indx[0] = ((new_centre[0] - ra)**2).argmin()
        ra_centre_new  = ra[centre_indx[0]]
        map.info['ra_centre']  = ra_centre_new
    if new_shape[0] != map.shape[1]:
        new_start = centre_indx[0] - new_shape[0]/2
        new_stop = centre_indx[0] + new_shape[0]/2
        map = map[:, new_start:new_stop, :]
        print "  select RA range: %d - %d"%(new_start, new_stop)
    
    
    '''>>> adjust dec <<<'''
    if new_centre[1] != None:
        centre_indx[1] = ((new_centre[1] - dec)**2).argmin()
        dec_centre_new  = dec[centre_indx[1]]
        map.info['dec_centre'] = dec_centre_new
    if new_shape[1] != map.shape[2]:
        new_start = centre_indx[1] - new_shape[1]/2
        new_stop = centre_indx[1] + new_shape[1]/2
        map = map[:, :, new_start:new_stop]
        print "  select DEC range: %d - %d"%(new_start, new_stop)
    
    map = algebra.make_vect(map)
    dec = map.get_axis('dec')
    ra = map.get_axis('ra')
    print "  old map:",
    print "centre at RA=%5.4f DEC=%5.4f"%(map.info['ra_centre'], 
                                          map.info['dec_centre'])
    print "           shape ",
    print map.shape
    print 

    output_root = output_file.replace(output_file.split('/')[-1], '')
    if not os.path.exists(output_root):
        os.makedirs(output_root)
    algebra.save(output_file, map)

if  __name__=="__main__":
    new_centre = (None, 2.0)
    new_shape = (128, 64)

    input_root = '/mnt/raid-project/gmrt/tcv/maps/1hr_41-18_avg_fdgp/'
    name_list = []
    for f in os.listdir(input_root):
        if f.split('.')[-1] != 'npy':
            continue
        f_split = f.split('_')
        if f_split[0][:-1] == 'sec' and \
           (f_split[5] == 'clean' or 
            f_split[5:8] == ['noise','inv','diag']):
            name_list.append(f)
            print f

    for input_name in name_list:
        #input_name = 'secA_1hr_41-18_avg_fdgp_clean_map_I_800.npy'
        input_file = input_root + input_name
        
        output_root = '/mnt/scratch-gl/ycli/maps/1hr_41-18_avg_fdgp/'
        output_name = input_name
        output_file = output_root + output_name

        map_cut(input_file, output_file, new_centre, new_shape)


