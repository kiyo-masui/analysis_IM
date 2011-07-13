import numpy as np

import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from simulations import lofar, pointsource

cm = matplotlib.cm.hot

ls = lofar.LofarGDSE()

tb = ls.getfield()
fs = ls.nu_pixels / ls.nu_0
tb = tb *fs[:,np.newaxis,np.newaxis]**2.55

mint = 19.0
maxt = 21.0

s_aa = tb[-1,:,:]
s_af = tb[:,0,:]

f = plt.figure(1, figsize = (12,5))
f.subplots_adjust(left=0.03, right = 0.94, top=0.98, bottom=0.15)

#grid = ImageGrid(f, 111, # similar to subplot(111)
#                 nrows_ncols = (1, 2),
#                 aspect = False,
#                 add_all=True,
#                 label_mode = "L",
#                 cbar_mode='single'
#                 )

ax1 = f.add_subplot(122)

ax1.imshow(s_aa.T, extent = (0, 5.0, 0, 5.0), aspect= 'equal', vmin = mint, vmax=maxt, cmap = cm)
ax1.set_xlabel("x / degrees")
ax1.set_ylabel("y / degrees")

ax2 = f.add_subplot(121)

im2 = ax2.imshow(s_af.T, extent = (120.0, 325.0, 0, 5.0), aspect= 'auto', vmin = mint, vmax=maxt, cmap = cm)
ax2.set_ylabel("y / degrees")
ax2.set_xlabel("Frequency / MHz")

#grid.cbar_axes[0].colorbar(im2)
cb = f.colorbar(im2, ax=ax1)
cb.set_label("$T_b$ / K")

f.savefig("foreground_syn.pdf")
f.clf()

cm = matplotlib.cm.gray

ps = pointsource.DiMatteo()
ps.x_num, ps.y_num = (256, 256)
ps.nu_lower, ps.nu_upper, ps.nu_num = (120.0, 325.0, 64)
#ps.flux_max = 50.0
psm = ps.getfield()


s_aa = psm[0,:,:]
s_af = psm.sum(axis = 1)


f = plt.figure(1, figsize = (12,5))

f.subplots_adjust(left=0.03, right = 0.94, top=0.98, bottom=0.15)
ax1 = f.add_subplot(122)

ax1.imshow(s_aa.T, extent = (0, 5.0, 0, 5.0), aspect= 'equal', cmap = cm, vmin = 0, vmax = 10)
ax1.set_ylabel("y / degrees")
ax1.set_xlabel("x / degrees")
ax2 = f.add_subplot(121)



im2 = ax2.imshow(s_af.T, extent = (120.0, 325.0, 0, 5.0), aspect= 'auto', cmap = cm, vmin = 0, vmax = 10)
ax2.set_ylabel("y / degrees")
ax2.set_xlabel("Frequency / MHz")

cb = f.colorbar(im2, ax=ax1)
cb.set_label("Flux / Jy")

f.savefig("foreground_ps.pdf")
f.clf()
