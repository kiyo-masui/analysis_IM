import h5py
import numpy as np
from scipy.special import erf
import random

# program takes mass of a halo and assign number of galaxies it contains.
# mapping is not 1 to 1 since we use random generator.

# list of parameters
log_M_min = 13.09
M_1 = 10**(14.)
M_0 = 10**(13.077)
sigma_log_M = 0.596
alpha = 1.0127

def ones(arr):
    for ind in range(len(arr)):
        if arr[ind] == True:
            arr[ind] = 1
        if arr[ind] == False:
            arr[ind] = 0
    return arr

def galaxy_calc(h5py_file):
    # open halo catalog
    file = h5py.File(h5py_file,'r+')
    halo_mass = file['Halo_Masses']['Halo_Masses'][:]

    # Notation: A_cen - average number of centre galaxies A_sat - average
    # number of satelite galaxies
    A_cen = 0.5*(1 + erf((np.log(halo_mass) - log_M_min)/sigma_log_M))
    A_sat = A_cen*((halo_mass - M_0)/M_1)**alpha

    # for simulation we get N_sat from the poisson distribution with mean value
    # equal to the A_sat
    N_sat = np.random.poisson(A_sat)

    # for central galaxy we use the probability of its occurence
    check = np.random.random(len(A_sat))
    bool = (check < A_sat)
    print "length of the array to convert:", len(bool)
    N_cen = ones(bool)
    galaxies = N_cen + N_sat

    # saving a new column with number of galaxies
    file.create_group("Galaxies_number")
    file["Galaxies_number"].create_dataset("Galaxies_number", data = galaxies)
    file.close()


if __name__ == '__main__':
    """
    redshift  = [0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500,
                 1.600, 1.700, 1.800, 1.900, 2.000, 2.100, 2.200]
    index = [0,1,2,3,4,5,6,7]
    """
    redshift = [0.800]
    index = [1,2,3,4,5,6,7]
    path = "/mnt/raid-project/gmrt/mufma/JD_h5py_catalogs/"
    for red in redshift:
        for ind in index:
            print "finished processing index", ind
            galaxy_calc(path + '%.3fhalo_catalog000%d.hdf5'%(red,ind))
        print "finished processing redshift", red

