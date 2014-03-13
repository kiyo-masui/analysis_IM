import numpy as np

from core import algebra
from map import beam

from utils import units

from simulations import corr21cm



thetax = 5.0
thetay = 3.0
nx = 128
ny = 64

f1 = 850
f2 = 650
nf = 256

cr = corr21cm.Corr21cm()

cr.x_width = thetax
cr.y_width = thetay
cr.x_num = nx
cr.y_num = ny

cr.nu_lower = f2
cr.nu_upper = f1
cr.nu_num = nf



rf = cr.getfield()


a = algebra.make_vect(rf, axis_names = ('freq', 'ra', 'dec'))
a.set_axis_info('freq', (f1+f2)/2.0, (f1-f2)/nf)
a.set_axis_info('ra', 0.0, thetax / nx)
a.set_axis_info('dec', 0.0, thetay / ny)

b = beam.GaussianBeam(width = [0.2, 0.2*f2/f1], freq = [f2, f1])
#ab = b.apply(a)













