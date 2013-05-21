from core import algebra
import numpy as np
from scipy import zeros
from math import sqrt

# opening files
weight_filename = open("/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/sec_A_15hr_41-90_noise_weight_I.npy", "r")
map_filename = open("/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_map_oldcal/sec_A_cleaned_clean_map_I_with_B_20modes.npy", "r")

# making a matrix
map = algebra.make_vect(algebra.load(map_filename
))
weight = algebra.make_vect(algebra.load(weight_filename))

# creating an array representing frequencies of the slices
freq_axis = weight.get_axis("freq")

#calculation of frequency binning
array = np.roll(freq_axis, 1) - freq_axis
delta_freq = array[4]

# creating a matrix for noise connecting to it an info package(row,col,...)
Tnoise = np.zeros(858624)
Tnoise.shape = (256, 78, 43)
Tnoise = algebra.make_mat(Tnoise,
row_axes=(0,), col_axes=(1, 2), axis_names=('freq', 'ra', 'dec'))

# calculating the total time in strange units
total_units = 0
for i in range(0, 78):
    for j in range(0, 43):
        if (weight[0][i][j]!=0):
            total_units = total_units + 1.0/weight[0][i][j]

# filling in the noise cube
for k in range(0, 256):
    for i in range(0, 78):
        for j in range(0, 43):
            if (weight[k][i][j]!=0):
                Tnoise[k][i][j] = 30.0*sqrt(1.0/(delta_freq*(360000.0/total_units)
                * (weight[k][i][j])))

# saving cube to the directory
file = open("/tmp/mufma/data/Tnoise.npy", "w")
algebra.save(file, Tnoise)
file.close()
