
import matplotlib
#matplotlib.use('PDF')
from matplotlib import pyplot as plt
#from mpl_toolkits.axes_grid1 import ImageGrid

#from matplotlib import rc
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

import numpy as np

from simulations import corr21cm, foregroundsck, lofar
from utils import units

import scipy.linalg as la

nf = 128
nul = 500.0
nuh = 700.0

freq = np.linspace(nul, nuh, nf)

z = units.nu21 / freq
z1, z2 = np.meshgrid(z, z)

cr = corr21cm.Corr21cm.from_file_matterps(fname = "data/corr0.dat", redshift = 0.0)

cs = cr.correlation_angular(0.0, z1, z2)

noise_power = 1e-7

fsyn = foregroundsck.Synchrotron()

cf = fsyn.angular_correlation(0.0) * fsyn.frequency_covariance(freq[:,np.newaxis], freq[np.newaxis,:])

cn = cf + np.identity(nf) * noise_power

evals, evecs = la.eigh(cs, cn)

cr2 = corr21cm.Corr21cm()

thetax = 5.0
thetay = 5.0
nx = 128
ny = 128

vs = cr2.realisation(z[-1], z[0], thetax, thetay, nf, nx, ny)[::-1,...]

ls = lofar.LofarGDSE()

ls.nu_lower, ls.nu_upper, ls.nu_num = (nul, nuh, nf)

vf = ls.getfield() * 1e3

x0 = vs + vf

def projectf(img, fvec):

    p = fvec[:,np.newaxis,np.newaxis] * np.tensordot(fvec, img, axes=(0, 0))[np.newaxis,:,:] / np.dot(fvec, fvec)
    return img - p

def projectn(img, vecs):
    ni = np.linalg.inv(np.tensordot(vecs, vecs, axes = (0,0)))

    proj = np.tensordot(vecs, np.tensordot(ni, vecs, axes = (1,1)), axes = (1,0))

    return img - np.tensordot(proj, img, axes = (1, 0))

def projectinv(img, evecs, n):
    y = np.tensordot(evecs, img, axes = (0,0))
    y[:n,:,:] = 0.0
    return np.tensordot(np.linalg.inv(evecs), y, axes = (0,0))
    
    
projx = { ('z%i' % i) : projectinv(x0, evecs, i) for i in range(10) }

np.savez('projarray.npy', **projx)

np.savez('foreground.npy', fg=vf)
np.savez('signal.npy', sg=vs)
np.savez('evecs.npy', evecs=evecs)
