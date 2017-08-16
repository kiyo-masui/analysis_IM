import numpy as np
from numpy import random
import core.algebra as al
from simulations import corr21cm
from utils import units

dir = '/scratch2/p/pen/nluciw/parkes/maps/correct_effic/I/ra%s/'
sim = '/scratch2/p/pen/nluciw/parkes/simulations/corr_effic/'
trans = '/scratch2/p/pen/nluciw/parkes/maps/transfer/'
out = '/scratch2/p/pen/nluciw/parkes/simulations/shuffle/ra%s/'
mapo = '/scratch2/p/pen/nluciw/parkes/maps/correct_effic/shuffle/ra%s/'


#types = ['delta', 'optsim', 'pks', 'raw']
types = ['optsim']
fields = ['165','199','33', 'n18', '225']
simnum = 100

c21 = corr21cm.Corr21cm(z=0.08, nu_upper=1400., nu_lower=1200.)

#freqmap = al.load(out%199 + 'sim_degradebeam_000.npy', 
#                  metafile=out%199 + 'sim_degradebeam_000.npy.meta')
#freqmap = al.make_vect(freqmap)

#z_axis = units.nu21 / (freqmap.get_axis('freq') / 1.e6) - 1.0
#T_hip = c21.T_b(z_axis) * 1.e-3
#T_alf = c21.T_b(z_axis, om_HI=4.3e-4) * 1.e-3

def shufflesim(types, simnum, dir):
    for i in range(simnum):
        optsim = dir + 'sim_%s_%03d.npy'%(types[0],i)
        map = al.load(optsim)
        np.random.shuffle(map)
        al.save(dir + 'sim_optsim_%03d.npy'%i, map)

def shufflemap(dir, sim_dir, tran, mapo, out, field):

    mapfile = dir%field + 'combined_clean_map_10modes.npy'
    weightfile = dir%field + 'combined_clean_weight_10modes.npy'
    map = al.load(mapfile)
    weight = al.load(weightfile)
    shuff_map = np.empty(map.shape, dtype=map.dtype)
    shuff_wt = np.empty(weight.shape, dtype=weight.dtype)
    ind = np.random.permutation(len(map))

    for old, new in enumerate(ind):
        shuff_map[new] = map[old]
        shuff_wt[new] = weight[old]

    shuff_map = al.make_vect(shuff_map)
    shuff_map.info = map.info
    shuff_wt = al.make_vect(shuff_wt)
    shuff_wt.info = weight.info
    al.save(mapo%field + 'combined_clean_map_10modes.npy', shuff_map)
    al.save(mapo%field + 'combined_clean_weight_10modes.npy', shuff_wt)

    for ii in range(100):

        sim = al.load(sim_dir + 'velo/ra%s/'%field + 'sim_pks_%03d.npy'%ii)
        trans = al.load(tran + 'correct_effic/ra%s/combined_clean_map_10modes_%03d.npy'%(field,ii))
        shuff_sim = np.empty(sim.shape, dtype=map.dtype)
        shuff_trans = np.empty(trans.shape, dtype=map.dtype)

        testit = 0
        for old, new in enumerate(ind):
            if not testit:
                print old, 'goes to', new
                testit += 1
            shuff_sim[new] = sim[old]
            shuff_trans[new] = trans[old]

        shuff_sim = al.make_vect(shuff_sim)
        shuff_sim.info = map.info
        al.save(out%field + 'sim_pks_%03d.npy'%ii, shuff_sim)
        shuff_trans = al.make_vect(shuff_trans)
        shuff_trans.info = map.info
        al.save(tran + 'shuffle/ra%s/combined_clean_map_10modes_%03d.npy'%(field, ii), shuff_trans)

def regmove(fields, dir, out):
    for i in range(100):
        for field in fields:
            map_dir = dir%(field) + 'sim_degradebeam_%03d.npy'%i
            out_dir = out%field + 'sim_degradebeam_%03d.npy'%i
            map = al.load(map_dir) / c21.T_b(0.08, om_HI=4.3e-4) * c21.T_b(0.08)
            map *= T_alf[:, np.newaxis, np.newaxis] / T_hip[:, np.newaxis, np.newaxis]
            print 'Saving %03d...'%i,
            al.save(out_dir, map)
            print 'Done'

for field in fields:
        shufflemap(dir, sim, trans, mapo, out, field)
#regmove(fields, out, out)
