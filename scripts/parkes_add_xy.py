import numpy as np
import core.algebra as al
import os

#map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/renorm_sync07/'

#fields = ['ra165/','ra182/','ra199/','ra216/','ran18/','ra33/']

#modes = [0,1,2,5,10,15,20,25]

map_ra = os.getenv('MAP_RA')
map_dec = os.getenv('MAP_DEC')
map_size = os.getenv('MAP_SIZE')

thebeams = {'XX' :['123','456','789','111213'], 'YY' : ['123','456','789','1012']}
finalbeams = ['123','456','789','10to13']

basemap = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/renorm_hitconv_sync07/combined/'
imap = basemap + 'I/'

def add_maps(root_dir, field, modes):
    xx_fname = root_dir + 'XX/' + field + 'combined_clean_map_' + str(modes) + 'modes.npy'
    yy_fname = root_dir + 'YY/' + field + 'combined_clean_map_' + str(modes) + 'modes.npy'
    I_map = 0.5*(al.load(xx_fname)+al.load(yy_fname))
    I_fname = root_dir + 'I/' + field + 'combined_clean_map_' + str(modes) + 'modes.npy'
    al.save(I_fname, I_map)

def add_weights(root_dir, field, modes):
    xx_fname = root_dir + 'XX/' + field + 'combined_clean_weight_' + str(modes) + 'modes.npy'
    yy_fname = root_dir + 'YY/' + field + 'combined_clean_weight_' + str(modes) + 'modes.npy'
    xx_weight = al.load(xx_fname)
    yy_weight = al.load(yy_fname)
    xx_weight[xx_weight<1.e-30] = 1.e-30
    yy_weight[yy_weight<1.e-30] = 1.e-30
    I_weight = 1./(0.25*(1./xx_weight+1./yy_weight))
    I_weight[I_weight<1.e-20] = 0
    I_fname = root_dir + 'I/' + field + 'combined_clean_weight_' + str(modes) + 'modes.npy'
    al.save(I_fname, I_weight)

if __name__ == '__main__':
#    for field in fields:
#        for mode in modes:
#            add_maps(map_dir, field, mode)
#            add_weights(map_dir, field, mode)
    for beams in zip(thebeams['XX'],thebeams['YY'],finalbeams):
        xx_map = al.load(basemap + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%(b)s_clean_map_bp_div_' % {"b" : beams[0]} + 'XX' +'_1316.npy')
        yy_map = al.load(basemap + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%(b)s_clean_map_bp_div_' % {"b" : beams[1]} + 'YY' +'_1316.npy')
        mapname = imap + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%(b)s_clean_map_bp_div_' % {"b" : beams[2]} + 'I' +'_1316.npy'
        map_out = 0.5*(xx_map+yy_map)
        #Zero out places where xx or yy don't have hits
        bool = np.logical_or(xx_map ==0.0, yy_map==0.0)
        map_out[bool] == 0.0
        al.save(mapname, map_out)
        xx_weight = al.load(basemap + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%(b)s_noise_inv_diag_bp_div_' % {"b" : beams[0]} + 'XX' +'_1316.npy')
        yy_weight = al.load(basemap + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%(b)s_noise_inv_diag_bp_div_' % {"b" : beams[1]} + 'YY' +'_1316.npy')
        weightname = imap + 'parkes_parallel_thread_' + map_ra + map_dec + '_p08_' + map_size + '_beam%(b)s_noise_inv_diag_bp_div_' % {"b" : beams[2]} + 'I' +'_1316.npy'
        xx_weight[xx_weight<1.e-30] = 1.e-30
        yy_weight[yy_weight<1.e-30] = 1.e-30
        I_weight = 1./(0.25*(1./xx_weight+1./yy_weight))
        I_weight[I_weight<1.e-20] = 0
        #For places where xx or yy map is zero (bool), set weight to min of xx and yy weights
        I_weight[bool] = np.minimum(xx_weight, yy_weight)[bool]
        al.save(weightname, I_weight) 
