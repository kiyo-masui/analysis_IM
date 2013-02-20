"""Make simulated signal realizations in the GBT survey volume"""
import numpy as np
import scipy as sp
from simulations import corr21cm
import core.algebra as algebra
import map.beam as beam
from core import constants as cc

thetax = 5.0
thetay = 3.0
nx = 64
ny = 32

freq_near = 850.
freq_far = 650.
num_freq = 256.

z_near = cc.freq_21cm_MHz / freq_near
z_far = cc.freq_21cm_MHz / freq_far

# TODO: interpolate const redshift onto const freq

a_list = []
ab_list = []
for i in range(0, 100):
    cr = corr21cm.Corr21cm()
    ## Slices produced by Corr21cm.realisation are equally spaced in
    ## redshift (ascending). Obviously this does not map onto regular
    ## frequency slices. For the moment, just crudely assume it does (and
    ## reverse the ordering to produce ascending in freq).
    rf = cr.realisation(z_near, z_far, thetax, thetay, num_freq, nx, ny)[::-1, ...]
    # the axis info for all stuff to be done on the 15hr field.
    info = {'ra_delta': -0.075045715822411624,
            'dec_delta': 0.074999999999999997,
            'dec_centre': 2.0,
            'axes': ('freq', 'ra', 'dec'),
            'ra_centre': 217.90000000000001,
            'freq_centre': 799609375.0,
            'freq_delta': -781250.0,
            'type': 'vect'}
    a = algebra.make_vect(rf, axis_names=('freq', 'ra', 'dec'))
    a = a[:, :64, :32]
    #a.set_axis_info('freq', (freq_near+freq_far)/2.0, (freq_near-freq_far)/float(num_freq))
    #a.set_axis_info('ra', 0.0, thetax / nx)
    #a.set_axis_info('dec', 0.0, thetay / ny)
    a.info = info
    a_list.append(a)
    # Real data beam stuff.
    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    freq_data *= 1.0e6
    b = beam.GaussianBeam(beam_data, freq_data)
    ab = b.apply(a)
    ab_list.append(ab)

for i in range(0, len(a_list)):
    root = '/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/test100/'
    save_name = root + "simulated_signal_map_" + str(i + 1) + ".npy"
    algebra.save(save_name, a_list[i])  # outmap)


for i in range(0, len(ab_list)):
    root = '/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/test100/'
    save_name = root + "simulated_signal_map_" + str(i + 1) + "_with_beam.npy"
    algebra.save(save_name, ab_list[i])  # outmap)

## THE MAPS IN calinliv/simulations/maps/ HAVE simulated_signal_map_1
## ADDED IN THEM AT STRENGTH 0.1 (hence the naming of the output file.

root = '/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/'
name = root + "simulated_signal_map_" + str(1)
loadnames = (name + ".npy", name + "_with_beam.npy")

exMap1 = algebra.make_vect(algebra.load(loadnames[0]))
exMap2 = algebra.make_vect(algebra.load(loadnames[1]))

root = '/mnt/raid-project/gmrt/kiyo/wiggleZ/maps/'
filename = root + "sec_A_15hr_41-73_clean_map_I.npy"
exMap = algebra.make_vect(algebra.load(filename))

exMap += 0.001 * exMap2

root = '/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/maps/'
algebra.save(root + 'sec_A_15hr_41-73_clean_map_I_with_point001_sim1.npy', exMap)

cr = corr21cm.Corr21cm()  # make sime and apply beam

## Slices produced by Corr21cm.realisation are equally spaced in
## redshift (ascending). Obviously this does not map onto regular
## frequency slices. For the moment, just crudely assume it does (and
## reverse the ordering to produce ascending in freq).
rf = cr.realisation(z_near, z_far, thetax, thetay, num_freq, nx, ny)[::-1, ...]

# the axis info for all stuff to be done on the 15hr field.
info = {'ra_delta': -0.075045715822411624,
        'dec_delta': 0.074999999999999997,
        'dec_centre': 2.0,
        'axes': ('freq', 'ra', 'dec'),
        'ra_centre': 217.90000000000001,
        'freq_centre': 799609375.0,
        'freq_delta': -781250.0,
        'type': 'vect'}

a = algebra.make_vect(rf, axis_names=('freq', 'ra', 'dec'))
#a.set_axis_info('freq', (freq_near+freq_far)/2.0, (freq_near-freq_far)/float(num_freq))
#a.set_axis_info('ra', 0.0, thetax / nx)
#a.set_axis_info('dec', 0.0, thetay / ny)
a.info = info


# Real data beam stuff.
beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
freq_data *= 1.0e6
b = beam.GaussianBeam(beam_data, freq_data)

#b = beam.GaussianBeam(width = [0.25, 0.25*freq_far/freq_near], freq = [freq_far, freq_near])
ab = b.apply(a)

#abc = copy.deepcopy(ab)

######### End beam sim ##########################################

for i in range(0, 10):
    plt.figure()
    plt.imshow(exMap1[i, :, :])
    plt.colorbar()
    plt.figure()
    plt.imshow(exMap2[i, :, :])
    plt.colorbar()


######### Get maps and add in sim or whatever you want ############

# After mode subtraction, take out a map and noise from pair.
Map1 = copy.deepcopy(F.Pairs[0].Map1)
Noise_inv1 = copy.deepcopy(F.Pairs[0].Noise_inv1)
# Make a noise for the beamsim same size all 1.
beam_noise = sp.zeros_like(Noise_inv1) + 1.0
# Make a new pair with the map and beamsim.
beam_map = copy.deepcopy(ab)
beam_map = algebra.make_vect(beam_map)
beam_map.info = copy.deepcopy(ab.info)
beam_map *= 0.1
# Map1 += 0.1*beam_map
p = fs.MapPair(Map1, 0.1 * beam_map, Noise_inv1, beam_noise, range(50, 65))
# Run correlation.
lags = F.params['lags']
corr, counts = p.correlate(lags)
