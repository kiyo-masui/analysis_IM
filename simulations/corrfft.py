from simulations import corr21cm

import numpy as np

import scipy


cr = corr21cm.Corr21cm()

la = np.logspace(2.0, 4.3, 100)[:,np.newaxis]

z1 = np.array([[1.0, 1.001]])
z2 = np.array([[1.0]])

aps1 = cr.angular_powerspectrum_fft(la, z1, z2)
aps2 = cr.angular_powerspectrum_flat(la, z1, z2)
