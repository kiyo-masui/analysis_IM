from core import algebra
import numpy as np
from scipy import zeros
from math import sqrt
import copy
import random

# opening files
weight_filename = open("/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/sec_A_15hr_41-90_noise_weight_I.npy","r")


# making a matrix
weight = algebra.make_vect(algebra.load(weight_filename)) + 0.0001

# creating an array representing frequencies of the slices
freq_axis = weight.get_axis("freq")


#calculation of frequency binning
array = np.roll(freq_axis, 1) - freq_axis
delta_freq = array[4]

# calculating the total time in strange units
total_units = np.sum(weight[10, :, :])

# creating a matrix for noise connecting to it an info package(row,col,...)
coef = 30.0/np.sqrt(delta_freq*(90000.0/total_units))
Tnoise = np.random.standard_normal(weight.shape)
Tnoise = Tnoise * coef * 1.0 / np.sqrt(weight)
Tnoise = algebra.make_vect(Tnoise, axis_names=('freq', 'ra', 'dec'))
Tnoise.copy_axis_info(weight)

# saving cube to the directory
file = open("/tmp/mufma/data/Instrumental_noise.npy","w")
algebra.save(file,Tnoise)
file.close()
