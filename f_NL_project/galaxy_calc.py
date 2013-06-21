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
    N_cen = np.array(bool, dtype=float)
    galaxies = N_cen + N_sat

    # saving a new column with number of galaxies
    file.create_group("Galaxies_number")
    file["Galaxies_number"].create_dataset("Galaxies_number", data = galaxies)
    file.close()


