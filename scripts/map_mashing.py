import numpy as np
import core.algebra as algebra
from utils import data_paths

def add_sim_to_data(simkey, datakey, replace=False):
    datapath_db = data_paths.DataPath()

    mapA_file = datapath_db.fetch(datakey + ":A;clean_map", intend_read=True)
    mapB_file = datapath_db.fetch(datakey + ":B;clean_map", intend_read=True)
    mapC_file = datapath_db.fetch(datakey + ":C;clean_map", intend_read=True)
    mapD_file = datapath_db.fetch(datakey + ":D;clean_map", intend_read=True)
    simfile = datapath_db.fetch(simkey + ":1", intend_read=True)

    simmap = algebra.make_vect(algebra.load(simfile))

    mapset = [mapA_file, mapB_file, mapC_file, mapD_file]
    for mapfile in mapset:
        print mapfile, simfile
        origmap = algebra.make_vect(algebra.load(mapfile))
        if replace:
            algebra.save(mapfile, simmap)
        else:
            algebra.save(mapfile, origmap + simmap)

#add_sim_to_data("sim_15hr_oldmap_str_beam", "GBT_15hr_map_fdgcal_plussim")
#add_sim_to_data("sim_15hr_oldmap_str_beam", "GBT_15hr_map_oldcal_plussim")
add_sim_to_data("sim_15hr_oldmap_str_beam", "GBT_15hr_map_signal_only", replace=True)
