#! /usr/bin/env python 
from matplotlib.transforms import Affine2D

from mpl_toolkits.axisartist.floating_axes import FloatingSubplot,\
     GridHelperCurveLinear

import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes

def setup_axes(fig):

    # rotate a bit for better orientation
    tr_rotate = Affine2D().translate(-17, 0)

    # scale degree to radians
    tr_scale = Affine2D().scale(np.pi/180., 1.)

    tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

    grid_locator1 = angle_helper.LocatorHMS(2)
    tick_formatter1 = angle_helper.FormatterHMS()

    #from mpl_toolkits.axes_grid.grid_finder import FixedLocator
    #grid_locator2 = FixedLocator([0., 5000, 10000, 15000])

    from mpl_toolkits.axisartist.grid_finder import MaxNLocator
    grid_locator2 = MaxNLocator(3)

    ra0, ra1 = 1.5*15, 2.5*15
    cz0, cz1 = 0.04, 0.1
    grid_helper = GridHelperCurveLinear(tr,
                                        extremes=(ra1, ra0, cz1, cz0),
                                        grid_locator1=grid_locator1,
                                        grid_locator2=grid_locator2,
                                        tick_formatter1=tick_formatter1,
                                        tick_formatter2=None,
                                        )

    ax1 = FloatingSubplot(fig, 111, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # adjust axis
    ax1.axis["left"].toggle(ticklabels=False)
    ax1.axis["right"].toggle(ticklabels=True)
    ax1.axis["right"].set_axis_direction("bottom")
    ax1.axis["right"].label.set_visible(True)
    #ax1.axis["right"].major_ticklabels.set_pad(5) #label.set_visible(True)

    ax1.axis["bottom"].major_ticklabels.set_axis_direction("top")
    ax1.axis["bottom"].label.set_axis_direction("top")

    ax1.axis["top"].set_visible(False)

    ax1.axis["right"].label.set_text(r"z ")
    ax1.axis["bottom"].label.set_text(r"$\alpha_{2000}$")

    #ax1.axis["right"].set_visible(False)
    #ax1.axis["bottom"].set_visible(False)
    #ax1.axis["left"].set_visible(False)

    # create a parasite axes whose transData in RA, cz
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                        # drawn twice, and possibly over some other
                        # artists. So, we decrease the zorder a bit to
                        # prevent this.

    return ax1, aux_ax

import matplotlib.pyplot as plt
from core import algebra
from utils import binning

root_map = '/Users/ycli/DATA/2df/map_2929.5/'
name_map = 'real_map_2df'
real_map = algebra.make_vect(algebra.load(root_map + name_map + '.npy'))
print real_map.shape
real_map_ra_z = np.ma.array(np.sum(real_map, axis=2))
real_map_ra_z[real_map_ra_z==0] = np.ma.masked
ra = real_map.get_axis('ra')
#dec = real_map.get_axis('dec')
freq = real_map.get_axis('freq')
z = 1420.e6/freq - 1.
ra_bin_edge = binning.find_edges(ra)
z_bin_edge = binning.find_edges(z)

real_map_coor = np.zeros(shape= (2,) + ra.shape + z.shape)
real_map_coor[0] = ra[:, None]
real_map_coor[1] =  z[None, :]
print real_map_coor.shape
real_map_coor = real_map_coor.reshape(2,-1)
real_map = np.sum(real_map, axis=-1)
print real_map.shape
real_map = real_map.T.flatten()
ras = []
decs = []
for i in range(real_map_coor.shape[1]):
    for j in range(int(real_map[i])):
        ras.append(real_map_coor[0][i])
        decs.append(real_map_coor[1][i])

fig = plt.figure(figsize=(10,5))
fig.clf()
ax, aux_ax = setup_axes(fig)
aux_ax.scatter(ras, decs, s=10, edgecolor='none', alpha=0.3)
#ax.pcolor(ra_bin_edge, z_bin_edge, real_map_ra_z, cmap='Greens')
#plt.xlim(xmin=ra_bin_edge.min(), xmax=ra_bin_edge.max())
#ax.set_clim(0, 2)
plt.savefig('./png/real_map.png')

root_cat = '/Users/ycli/DATA/2df/catalogue/'
name_cat = 'real_catalogue_2df'
real_cat = np.loadtxt(root_cat + name_cat + '.out')
fig = plt.figure(figsize=(10,5))
fig.clf()
ax, aux_ax = setup_axes(fig)
aux_ax.scatter(real_cat[:,0]*180./np.pi, real_cat[:,2], s=4, 
               edgecolor='none', alpha=0.5)
plt.savefig('./png/real_cat.png')
#plt.show()

#
##ax.plot(real_cat[:,0], real_cat[:,2], 'r.')
##im = ax.pcolor(ra_bin_edge, z_bin_edge,  real_map_ra_z, cmap='Greens')
#plt.scatter(real_cat[:,0], real_cat[:,2], s=10, c='0.5', marker='.', 
#           alpha=0.7,
#           edgecolor='none')
#plt.xlim(xmin=ra_bin_edge.min(), xmax=ra_bin_edge.max())
#
##ax.set_rlim(rmin=0., rmax=z_bin_edge[-1])
##ax.set_rlim(rmin=np.min(real_cat[:,2]), rmax=np.max(real_cat[:,2]))
#
##plt.colorbar(im)
#plt.show()
#
