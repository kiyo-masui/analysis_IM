
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
#from mpl_toolkits.axes_grid1 import ImageGrid

from matplotlib import gridspec

import numpy as np

def projectinv(img, evecs, n):
    y = np.tensordot(evecs, img, axes = (0,0))
    y[:n,:,:] = 0.0
    return np.tensordot(np.linalg.inv(evecs), y, axes = (0,0))


projx = np.load('projarray.npy.npz')

vf = np.load('foreground.npy.npz')['fg']
vs = np.load('signal.npy.npz')['sg']

evecs = np.load('evecs.npy.npz')['evecs']


#gs = gridspec.GridSpec(1, 2)

f = plt.figure()

ax1 = f.add_subplot(111)
#ax2 = plt.subplot(gs[0,1])
ax1.set_xticks([])
ax1.set_yticks([])
#ax2.set_xticks([])
#ax2.set_yticks([])

zfm = [projx['z%i'%i][0] - projx['z%i'%i][0].mean() for i in range(10)]
zff = [projx['z%i'%i][:,0,:] for i in range(10)]

vsm = vs[0] - vs[0].mean()
vfm = vf[0] - vf[0].mean()

for i in range(10):
    ax1.imshow(zfm[i])
    f.savefig('klplots/projx%i.pdf' % i)
    #f.clf()

ax1.imshow(vs[0])
f.savefig('klplots/signal.pdf')
    
projs = [ projectinv(vs, evecs, i) for i in range(10) ]
projf = [ projectinv(vf, evecs, i) for i in range(10) ]


#for i in range(10):
#    ax1.imshow(projs[i][0] - projs[i][0].mean())
#    f.savefig('klplots/s%i.pdf' % i)
#    #f.clf()

    
#for i in range(10):
#    ax1.imshow(projf[i][0] - projf[i][0].mean(), vmin = -5e2, vmax = 5e2)
#    f.savefig('klplots/f%i.pdf' % i)
#    #f.clf()


f.clf()
ax1 = f.add_subplot(121)
ax2 = f.add_subplot(122)

ax1.set_xticks([])
ax1.set_yticks([])
ax2.set_xticks([])
ax2.set_yticks([])


for i in range(4):
    ax1.imshow(zfm[i])
    ax2.imshow(zff[i])
    f.savefig('klplots/two%i.pdf' % i)
    #f.clf()

ax1.imshow(zfm[3])
ax2.imshow(vs[0])
f.savefig('klplots/comp_noproj.pdf')
    
for i in range(4,8):
    ax1.imshow(zfm[i])
    ax2.imshow(projs[i][0] - projs[i][0].mean() )
    f.savefig('klplots/comp%i.pdf' % i)
    #f.clf()
