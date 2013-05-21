import h5py
import numpy as np
import math

# list of parameters
log_M_min = 13.09
M_1 = 10**(14.)
M_0 = 10**(13.077)
sigma_log_M = 0.596
alpha = 1.0127

# open halo catalog
file = h5py.File(h5py_file)
Halo_mass = file['Halo_masses']['Halo_masses'][:]



# Notation: A_cen - average number of centre galaxies A_sat - average number of
# satelite galaxies
A_cen = 0.5*(1+math.erf((log(M) - log_M_min)/sigma_log_M))
A_sat = A_cen*((M - M_0)/M_1)**alpha

# for simulation we get N_sat from the poisson distribution with mean value
# equal to the A_sat
N_sat = np.random.poisson(A_sat)

# for central galaxy we use the probability of its occurence ??

# saving a new column with number of galaxies
file.create_group("Galaxies_number")
file["Galaxies_number"].create_dataset("Galaxies_number", data = )



