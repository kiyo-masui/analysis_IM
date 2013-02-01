#! /usr/bin/env python 

import os
import numpy as np
import data_paths
import scipy as sp


datadb = data_paths.DataPath()
sim_t = datadb.fetch('sim_15hr_oldmap_str_temperature')
sim_d = datadb.fetch('sim_15hr_oldmap_str_delta')

for i in range(10):
    maproot_t = sim_t[1]['%d'%i]
    maproot_d = sim_d[1]['%d'%i]

    map_t = np.load(maproot_t)
    map_d = np.load(maproot_d)

    tb =  map_t/map_d

    print tb[50]
