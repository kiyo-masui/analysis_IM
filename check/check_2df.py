#! /usr/bin/env python 
from matplotlib.transforms import Affine2D

from mpl_toolkits.axisartist.floating_axes import FloatingSubplot,\
     GridHelperCurveLinear

import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from map import physical_gridding as gridding

def setup_axes(fig, ra0, ra1, cz0, cz1, label0=r'$\alpha[2000]$', label1='z'):

    rotate_angle = 90. - 0.5*(ra0 + ra1)
    # rotate a bit for better orientation
    tr_rotate = Affine2D().translate(rotate_angle, 0)

    # scale degree to radians
    tr_scale = Affine2D().scale(np.pi/180., 1.)

    tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

    #grid_locator1 = angle_helper.LocatorHMS(2)
    #tick_formatter1 = angle_helper.FormatterHMS()

    grid_locator1 = angle_helper.LocatorD(4)
    tick_formatter1 = None

    #from mpl_toolkits.axes_grid.grid_finder import FixedLocator
    #grid_locator2 = FixedLocator([0., 5000, 10000, 15000])

    from mpl_toolkits.axisartist.grid_finder import MaxNLocator
    grid_locator2 = MaxNLocator(3)

    #ra0, ra1 = 1.5*15, 2.5*15
    #cz0, cz1 = 0.04, 0.1
    ra0, ra1 = ra0, ra1
    cz0, cz1 = cz0, cz1 
    grid_helper = GridHelperCurveLinear(tr,
                                        extremes=(ra1, ra0, cz1, cz0),
                                        grid_locator1=grid_locator1,
                                        grid_locator2=grid_locator2,
                                        tick_formatter1=tick_formatter1,
                                        tick_formatter2=None,
                                        )

    ax1 = FloatingSubplot(fig, 131, grid_helper=grid_helper)
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

    ax1.axis["bottom"].label.set_text(label0)
    ax1.axis["right"].label.set_text(label1)

    #ax1.axis["right"].set_visible(False)
    #ax1.axis["bottom"].set_visible(False)
    #ax1.axis["left"].set_visible(False)

    # create a parasite axes whose transData in RA, cz
    aux_ax1 = ax1.get_aux_axes(tr)

    aux_ax1.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                        # drawn twice, and possibly over some other
                        # artists. So, we decrease the zorder a bit to
                        # prevent this.

    ax2 = FloatingSubplot(fig, 132, grid_helper=grid_helper)
    fig.add_subplot(ax2)

    # adjust axis
    ax2.axis["left"].toggle(ticklabels=False)
    ax2.axis["right"].toggle(ticklabels=True)
    ax2.axis["right"].set_axis_direction("bottom")
    ax2.axis["right"].label.set_visible(True)
    #ax2.axis["right"].major_ticklabels.set_pad(5) #label.set_visible(True)

    ax2.axis["bottom"].major_ticklabels.set_axis_direction("top")
    ax2.axis["bottom"].label.set_axis_direction("top")

    ax2.axis["top"].set_visible(False)

    #ax2.axis["bottom"].label.set_text(r"$\alpha_{2000}$")
    #ax2.axis["right"].label.set_text(r"z ")
    ax2.axis["bottom"].label.set_text(label0)
    ax2.axis["right"].label.set_text(label1)

    #ax2.axis["right"].set_visible(False)
    #ax2.axis["bottom"].set_visible(False)
    #ax2.axis["left"].set_visible(False)

    # create a parasite axes whose transData in RA, cz
    aux_ax2 = ax2.get_aux_axes(tr)

    aux_ax2.patch = ax2.patch # for aux_ax to have a clip path as in ax
    ax2.patch.zorder=0.9 # but this has a side effect that the patch is
                        # drawn twice, and possibly over some other
                        # artists. So, we decrease the zorder a bit to
                        # prevent this.

    ax3 = plt.subplot(133)

    return ax1, aux_ax1, ax2, aux_ax2, ax3

def average_map(map_root_list):

    map_list = []
    for root in map_root_list:
        map = algebra.make_vect(algebra.load(root))
        info = map.info
        map_list.append(map)
    map = np.array(map_list)
    map = np.mean(map, axis=0)
    map = algebra.make_vect(map, info['axes'])
    map.info = info

    return map

def plot_sky(cat_root, map_root_list, save_name='./png/real_real_cat.png'):

    real_cat = np.loadtxt(cat_root)
    
    if len(map_root_list) == 1:
        real_map = algebra.make_vect(algebra.load(map_root_list[0]))
    else:
        real_map = average_map(map_root_list)
    
    real_box, real_box_info = gridding.physical_grid(real_map, refinement=1)
    
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
    real_map = 10. * np.sum(real_map, axis=-1)
    print real_map.shape
    real_map = real_map.T.flatten()
    ras = []
    decs = []
    for i in range(real_map_coor.shape[1]):
        for j in range(int(real_map[i])):
            ras.append(real_map_coor[0][i])
            decs.append(real_map_coor[1][i])
    
    real_box = algebra.make_vect(real_box, real_box_info['axes'])
    real_box.info = real_box_info
    y = binning.find_edges(real_box.get_axis('freq'))
    x = binning.find_edges(real_box.get_axis('ra'))
    
    fig = plt.figure(figsize=(15,10))
    fig.clf()
    
    ax1, aux_ax1, ax2, aux_ax2, ax3 =\
        setup_axes(fig, 1.5*15, 2.5*15, 0.01, 0.3)
    
    aux_ax1.scatter(real_cat[:,0]*180./np.pi, real_cat[:,2], s=4, 
                   edgecolor='none', alpha=0.5)
    
    aux_ax2.scatter(ras, decs, s=10, edgecolor='none', alpha=0.03)
    
    real_box = np.ma.sum(real_box, axis=2)
    real_box = np.ma.array(real_box)
    real_box[real_box==0] = np.ma.masked
    im = ax3.pcolormesh(x, y, real_box, )
    plt.colorbar(im)
    
    plt.savefig(save_name)

if __name__=="__main__":    
    
    import matplotlib.pyplot as plt
    from core import algebra
    from utils import binning
    import os
    
    
    root_cat = '/mnt/scratch-gl/ycli/2df_catalog/catalog/'

    root_map = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5/'

    name_cat = 'real_catalogue_2df'
    name_map = 'real_map_2df_delta'

    plot_sky(root_cat + name_cat + '.out', 
             [root_map + name_map + '.npy', ], 
             save_name='./png/real_cat.png')


    name_cat = 'mock_catalogue_2df_090'
    name_map = 'mock_map_2df_delta_090'

    plot_sky(root_cat + name_cat + '.out', 
             [root_map + name_map + '.npy',], 
             save_name='./png/mock_cat_090.png')

    name_cat = 'mock_catalogue_2df_090'
    map = []
    for i in range(100):
        map.append(root_map + 'mock_map_2df_delta_%03d.npy'%i)

    plot_sky(root_cat + name_cat + '.out', 
             map, 
             save_name='./png/mock_cat_average.png')


    name_cat = 'mock_catalogue_2df_090'
    name_map = 'sele_map_2df'

    plot_sky(root_cat + name_cat + '.out', 
             [root_map + name_map + '.npy',], 
             save_name='./png/sele_cat.png')

    name_cat = 'mock_catalogue_2df_090'
    name_map = 'sele_map_2df_separable'

    plot_sky(root_cat + name_cat + '.out', 
             [root_map + name_map + '.npy',], 
             save_name='./png/sele_cat_separable.png')




