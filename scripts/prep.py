from mkpower import add_beams2sim as a2s
import os
import sys
import numpy as np

sim_dir = '/scratch2/p/pen/nluciw/parkes/simulations/'
map_dir = '/scratch2/p/pen/andersoc/second_parkes_pipe/maps_bp_divide/correct_beam_effic/sync27/combined/I/'
trans_dir = '/scratch2/p/pen/nluciw/parkes/maps/'

dirty_dir = '/scratch2/p/pen/nluciw/parkes/simulations/dirty/'
clean_dir = '/scratch2/p/pen/nluciw/parkes/simulations/clean/'

fields = ['ra199', 'ra33', 'ra165', 'ran18']
beams = ['123', '456', '789', '10to13']

if str(sys.argv[1]) == 'pre':

    print 'Prepping for map clean...'
    print 'Using', map_dir

    for field in fields:
    
        if not os.path.exists(dirty_dir+field+'/'):
            os.makedirs(dirty_dir+field+'/')
        a2s.combine_for_clean(sim_dir+'/corr_effic/velo/', map_dir, dirty_dir, field, 
                beams, 100) 

if str(sys.argv[1]) == 'post':

    print 'Post cleaning processing...'

    for field in fields:
        for ii in range(100):
            clean_sim = np.load(clean_dir+'sim%03d/%s/combined_clean_map_10modes.npy'%(ii,field))
            np.save(sim_dir+'corr_effic/velo/'+field+'/sim_pks_%03d'%ii, clean_sim)

            map = np.load(trans_dir + \
                    'correct_effic/I/%s/combined_clean_map_10modes.npy'%field)
            np.save(trans_dir+'transfer/correct_effic/%s/combined_clean_map_10modes_%03d'%(field,ii),\
                    clenan_sim - map)

