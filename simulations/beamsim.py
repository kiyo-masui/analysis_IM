import numpy as np

import corr21cm

import core.algebra as algebra
import map.beam as beam

import units

thetax = 5.0
thetay = 3.0
nx = 128
ny = 64

f1 = 850
f2 = 650
nf = 256

z1 = units.nu21 / f1
z2 = units.nu21 / f2


cr = corr21cm.Corr21cm()

## Slices produced by Corr21cm.realisation are equally spaced in
## redshift (ascending). Obviously this does not map onto regular
## frequency slices. For the moment, just crudely assume it does (and
## reverse the ordering to produce ascending in freq).
rf = cr.realisation(z1, z2, thetax, thetay, nf, nx, ny)[::-1,...]


a = algebra.make_vect(rf, axis_names = ('freq', 'ra', 'dec'))
a.set_axis_info('freq', (f1+f2)/2.0, (f1-f2)/nf)
a.set_axis_info('ra', 0.0, thetax / nx)
a.set_axis_info('dec', 0.0, thetay / ny)

b = beam.GaussianBeam(width = [0.2, 0.2*f2/f1], freq = [f2, f1])
#ab = b.apply(a)













