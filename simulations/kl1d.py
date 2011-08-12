
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import numpy as np

from corr import RedshiftCorrelation

import foreground

import scipy.linalg as la

#import cubicspline as cs
from cosmology import Cosmology
from units import nu21


def tb(z):
    return 23.8 * ((1.0 + z) / 10.0)**0.5

def xh(z):
    return 0.02

#ps = cs.LogInterpolater.fromfile( join(dirname(__file__),"data/ps.dat")).value

nf = 100
nul = 650.0
nuh = 750.0

freq = np.linspace(nul, nuh, nf)

z = nu21 / freq

zm = z.mean()

cr = RedshiftCorrelation.from_file_matterps(fname  = "data/corr0.dat")

z1, z2 = np.meshgrid(z, z)

cs = xh(zm)**2 * tb(zm)**2 * cr.correlation_angular(0.0, z1, z2) / (1+zm)**2

noise_power = 1e-7

cn = foreground.cf_all2(freq, freq)

cn = cn + np.identity(nf) * noise_power

evals, evecs = la.eigh(cs, cn)



foreground.fgs = [ ('syn', 1000.0, 2.90, 2.5, 4.0),
                   ('ps',   57.0, 2.07, 1.1, 1.0),
                   ('eff', 0.014, 2.10, 1.0, 35.0),
                   ('gff', 0.088, 2.15, 3.0, 35.0),
                   ('rnd', 250.0, 2.90, 2.5, 0.2)]


cn2 = foreground.cf_all2(freq, freq)

cn2 = cn2 + np.identity(nf) * noise_power

smodes = evals
nmodes = np.dot(evecs.T, np.dot(cn, evecs)).diagonal()

n2modes = np.dot(evecs.T, np.dot(cn2, evecs)).diagonal()

f = plt.figure(1, figsize = (14,12))
f.subplots_adjust(left=0.08, right = 0.9, top=0.9, bottom=0.07, wspace=0.2, hspace=0.2)

ax1 = f.add_subplot(221)

sh = ax1.semilogy(smodes)
nh = ax1.semilogy(nmodes)
nh2 = ax1.semilogy(n2modes)

ax1.set_xlabel("Mode number")
ax1.set_ylabel("Signal/Noise ratio")
ax1.set_title("SN values and Power leakage")
ax2 = f.add_subplot(222)


for i in range(1,8):
    ax2.plot(freq, evecs[:,-i], "b")

ax2.set_xlim(nul, nuh)
ax2.set_xlabel("Frequency / Mhz")
ax2.set_ylabel("Weight")
ax2.set_title("Highest signal modes")

nvals, nvecs = la.eigh(cn)
n2vals, n2vecs = la.eigh(cn2)

ax3 = f.add_subplot(223)
ax4 = f.add_subplot(224)


for i in range(1,8):
    ax3.plot(freq, nvecs[:,-i], "g")
    #ax3.plot(freq, evecs[:,i-1], "r")
    ax4.plot(freq, n2vecs[:,-i], "r")


ax3.set_xlim(nul, nuh)
ax3.set_xlabel("Frequency / Mhz")
ax3.set_ylabel("Weight")

ax3.set_title("Top Foreground Eigenmodes")

ax4.set_xlim(nul, nuh)
ax4.set_xlabel("Frequency / Mhz")
ax4.set_ylabel("Weight")

ax4.set_title("Top Small-scale Foreground Eigenmodes")



lg = plt.figlegend((sh, nh, nh2), ("Signal", "Assumed Foregrounds", "Foregrounds (small-scale)"), "upper right", prop={'size': 'small'})
lg.draw_frame(False)


f.savefig("kl1d.pdf")
