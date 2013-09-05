import h5py
import numpy as np
from scipy.special import erf
import random

# Program takes mass of a halo and assign number of galaxies it contains.
# Mapping is not 1 to 1 since we use random generator.

# List of parameters from the BOSS experiment.
log_M_min = 13.09
M_1 = 10**(14.)
M_0 = 10**(13.077)
sigma_log_M = 0.596
alpha = 1.0127


def galaxy_calc(h5py_file):
    # open halo catalog
    print "Opening h5py catalog with halos"
    file = h5py.File(h5py_file,'r+')
    halo_mass = file['Halo_Masses']['Halo_Masses'][:]
    # Notation: A_cen - average number of centre galaxies A_sat - average
    # number of satelite galaxies
    A_cen = 0.5*(1 + erf((np.log10(halo_mass) - log_M_min)/sigma_log_M))
    A_sat = A_cen*((halo_mass - M_0)/M_1)**alpha

    # for simulation we get N_sat from the poisson distribution with mean value
    # equal to the A_sat
    N_sat = np.random.poisson(A_sat)

    # for central galaxy we use the probability of its occurence
    check = np.random.random(len(A_cen))
    bool = (check < A_cen)
    N_cen = np.array(bool, dtype=float)
    galaxies = N_cen + N_sat

    # saving a new column with number of galaxies
    print "Creating group for galaxies and saving data"
    try:
        del file['num_galaxies']
    except:
        pass
    file.create_group("num_galaxies")
    file["num_galaxies"].create_dataset("num_galaxies", data = galaxies)
    file.close()


