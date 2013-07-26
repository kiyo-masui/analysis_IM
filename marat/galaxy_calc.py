import h5py
import numpy as np
from scipy.special import erf
import random

# program takes mass of a halo and assign number of galaxies it contains.
# mapping is not 1 to 1 since we use random generator.

# list of parameters
# real one from BOSS
log_M_min = 13.09
#log_M_min = 11.5
M_1 = 10**(14.)
M_0 = 10**(13.077)
sigma_log_M = 0.596
alpha = 1.0127

def galaxy_calc(h5py_file):
    # open halo catalog
    print "Check?"
    file = h5py.File(h5py_file,'r+')
    print "Done!"
    halo_mass = file['Halo_Masses']['Halo_Masses'][:]
    # Notation: A_cen - average number of centre galaxies A_sat - average
    # number of satelite galaxies
    A_cen = 0.5*(1 + erf((np.log10(halo_mass) - log_M_min)/sigma_log_M))
    A_sat = A_cen*((halo_mass - M_0)/M_1)**alpha
    A_sat[np.isnan(A_sat)] = 0 #it works without this line, but is safer that way

    print log_M_min
    print A_cen.min(), A_cen.max()
    print A_sat.min(), A_sat.max()
    # for simulation we get N_sat from the poisson distribution with mean value
    # equal to the A_sat
    N_sat = np.random.poisson(A_sat)

    # for central galaxy we use the probability of its occurence
    check = np.random.uniform(size=len(A_sat))
    bool = (check < A_cen)
    print "length of the array to convert:", len(bool)
    N_cen = np.array(bool, dtype=float)
    galaxies = N_cen + N_sat

    # saving a new column with number of galaxies
    print "reached num_gal"
    if 'num_galaxies' in file.keys():
        file['num_galaxies']['num_galaxies'][:] = galaxies
    else:
        file.create_group("num_galaxies")
        file["num_galaxies"].create_dataset("num_galaxies", data = galaxies)
    file.close()
    
if __name__=='__main__':
    path_cat = '/mnt/raid-project/gmrt/goblot/JD_catalogs/'
    log_M_min = 10
    for i in range(8):
        print 'Step',i
        galaxy_calc(path_cat+'0.800halo_catalog000%s.hdf5'%i)


