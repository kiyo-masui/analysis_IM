import galaxy_calc as gc
import hdf5_to_algebra as ha
from core import algebra

for ind in [0,1,2,3,4,5,6,7]:
    print "Reading", '0.800halo_catalog000%d.hdf5'%ind
    gc.galaxy_calc('/mnt/raid-project/gmrt/mufma/0.800halo_catalog000%d.hdf5'%ind)



print "Starting to bin maps"
binned_map, info = ha.bin_val_map('/mnt/raid-project/gmrt/mufma/0.800halo_catalog0000.hdf5', 513, 'num_galaxies', 'num_galaxies')

for num in [1,2,3]:
        binned_map = binned_map + ha.bin_val_map('/mnt/raid-project/gmrt/mufma/0.800halo_catalog000%d.hdf5'%num,
513,'num_galaxies','num_galaxies')[0]
# save data for pipeline in manager.py
print "Starting to make algebra vector"
map = algebra.make_vect(binned_map, axis_names=('ra', 'dec', 'freq'))
print "Reading info"
map.info = info
print "Opened file for saving"
save_file = open('gal_val_512.npy', "w")
print "Saving"
algebra.save(save_file, map)
save_file.close()
