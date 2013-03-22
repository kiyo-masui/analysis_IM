#! /usr/bin/env python 

import os
import numpy as np
import data_paths
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid


#datadb = data_paths.DataPath()
#sim_t = datadb.fetch('sim_15hr_oldmap_str_temperature')
#sim_d = datadb.fetch('sim_15hr_oldmap_str_delta')
#
#for i in range(10):
#    maproot_t = sim_t[1]['%d'%i]
#    maproot_d = sim_d[1]['%d'%i]
#
#    map_t = np.load(maproot_t)
#    map_d = np.load(maproot_d)
#
#    tb =  map_t/map_d
#
#    print tb[50]

cmin = -9
cmax = -3


file_root = '/mnt/raid-project/gmrt/ycli/ps_result/reference_cros_15hr_ABCD_legendre_modes_0gwj_14conv_new_20/'

sim = np.load(file_root + "simmaps_p2_combined.npy")
sim = np.ma.array(sim)
sim[sim==0] = np.ma.masked
sim = np.ma.log10(sim)
fig = plt.figure(figsize=(25,15))
ax = ImageGrid(fig, 111,
               nrows_ncols = (1, 4),
               direction = "row",
               axes_pad = 0.3,
               add_all = True,
               label_mode = "L",
               share_all = True,
               cbar_location = "right",
               cbar_mode = "each",
               cbar_size = "5%",
               cbar_pad = 0.00,
               )

sim_deg = np.load(file_root + "simmaps_beam_p2_combined.npy")
sim_deg = np.ma.array(sim_deg)
sim_deg[sim_deg==0] = np.ma.masked
sim_deg = np.ma.log10(sim_deg)

im0 = ax[0].pcolormesh(sim)
im0.set_clim(cmin, cmax)
ax[0].cax.colorbar(im0)
ax[0].set_title('raw simulation maps \n $P(sim_{raw} x sim_{delta})$')

im1 = ax[1].pcolormesh(sim_deg)
im1.set_clim(cmin, cmax)
ax[1].cax.colorbar(im0)
ax[1].set_title('simulation maps\n1.4 common beam convolved\n$P(Conv(sim_{raw}) x sim_{delta})$')


file_root = '/mnt/raid-project/gmrt/ycli/ps_result/bias/cros_15hr_ABCD_legendre_modes_0gwj_14conv_new_40_subreal/'
sim_clean = np.load(file_root + "simmaps1_p2_combined.npy")
sim_clean = np.ma.array(sim_clean)
sim_clean[sim_clean==0] = np.ma.masked
sim_clean = np.ma.log10(sim_clean)

im2 = ax[2].pcolormesh(sim_clean)
im2.set_clim(cmin, cmax)
ax[2].cax.colorbar(im0)
ax[2].set_title('simulation maps\n1.4 common beam convolved\ncleaned subtract realmap\n$P((Cleaned(Conv(sim_{raw}))-Cleaned(real)) x sim_{delta})$')

bias = sim - sim_clean
im3 = ax[3].pcolormesh(bias)
#im3.set_clim(cmin, cmax)
ax[3].cax.colorbar(im3)
ax[3].set_title('simulation maps\n1.4 common beam convolved\ncleaned subtract realmap\n$\\frac{P(sim_{raw} x sim_{delta})}{P((Cleaned(Conv(sim_{raw}))-Cleaned(real)) x sim_{delta})}$')

#bias_old = np.load(file_root + "b2_bias.npy")
#bias_old = np.ma.array(bias_old)
#bias_old[bias_old==0] = np.ma.masked
#bias_old = np.ma.log10(bias_old)
#im4 = ax[4].pcolormesh(bias_old)
#ax[4].cax.colorbar(im4)

plt.savefig('./png/sim.png')
