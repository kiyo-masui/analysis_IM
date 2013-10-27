#! /usr/bin/env python 
from matplotlib.transforms import Affine2D

from mpl_toolkits.axisartist.floating_axes import FloatingSubplot,\
     GridHelperCurveLinear

import numpy as np
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from map import physical_gridding as gridding

import matplotlib.pyplot as plt
from core import algebra
from utils import binning
import os

def setup_axes_ra_dec(fig, ra0, ra1, dec0, dec1, label0=r'$\alpha[2000]$', 
                      label1=r'$\delta$'):
    if dec0 < 0:
        rotate_angle = - 0.5*(ra0 + ra1) - 90.
    else:
        rotate_angle = 90. - 0.5*(ra0 + ra1)
    # rotate a bit for better orientation
    tr_rotate = Affine2D().translate(rotate_angle, 0)

    # scale degree to radians
    tr_scale = Affine2D().scale(np.pi/180., np.pi/180.)

    tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

    grid_locator1 = angle_helper.LocatorD(4)
    tick_formatter1 = None

    grid_locator2 = angle_helper.LocatorD(4)
    tick_formatter2 = None

    grid_helper = GridHelperCurveLinear(tr,
                                        extremes=(ra1, ra0, dec1, dec0),
                                        grid_locator1=grid_locator1,
                                        grid_locator2=grid_locator2,
                                        tick_formatter1=tick_formatter1,
                                        tick_formatter2=tick_formatter2,
                                        )

    ax1 = FloatingSubplot(fig, 211, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # adjust axis
    ax1.axis["left"].toggle(ticklabels=True)
    ax1.axis["right"].toggle(ticklabels=False)
    ax1.axis["top"].toggle(ticklabels=True)
    ax1.axis["bottom"].toggle(ticklabels=False)
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["right"].set_axis_direction("bottom")
    ax1.axis["left"].label.set_visible(True)
    ax1.axis["top"].label.set_visible(True)
    ax1.axis["right"].label.set_visible(False)
    #ax1.axis["right"].major_ticklabels.set_pad(5) #label.set_visible(True)

    ax1.axis["top"].major_ticklabels.set_axis_direction("bottom")
    ax1.axis["top"].label.set_axis_direction("bottom")
    #ax1.axis["bottom"].major_ticklabels.set_axis_direction("top")
    #ax1.axis["bottom"].major_ticklabels.set_axis_direction("top")
    #ax1.axis["bottom"].label.set_axis_direction("bottom")

    #ax1.axis["top"].set_visible(False)

    #ax1.axis["bottom"].label.set_text(label0)
    ax1.axis["top"].label.set_text(label0)
    #ax1.axis["right"].label.set_text(label1)
    ax1.axis["left"].label.set_text(label1)

    aux_ax1 = ax1.get_aux_axes(tr)

    aux_ax1.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                        # drawn twice, and possibly over some other
                        # artists. So, we decrease the zorder a bit to
                        # prevent this.

    ax2 = FloatingSubplot(fig, 212, grid_helper=grid_helper)
    fig.add_subplot(ax2)

    # adjust axis
    ax2.axis["left"].toggle(ticklabels=False)
    ax2.axis["right"].toggle(ticklabels=True)
    ax2.axis["right"].set_axis_direction("bottom")
    ax2.axis["right"].label.set_visible(True)
    #ax2.axis["right"].major_ticklabels.set_pad(5) #label.set_visible(True)

    ax2.axis["bottom"].major_ticklabels.set_axis_direction("top")
    ax2.axis["bottom"].label.set_axis_direction("top")

    #ax2.axis["top"].set_visible(False)

    ax2.axis["bottom"].label.set_text(label0)
    ax2.axis["right"].label.set_text(label1)

    aux_ax2 = ax2.get_aux_axes(tr)

    aux_ax2.patch = ax2.patch # for aux_ax to have a clip path as in ax
    ax2.patch.zorder=0.9 # but this has a side effect that the patch is
                        # drawn twice, and possibly over some other
                        # artists. So, we decrease the zorder a bit to
                        # prevent this.


    return ax1, aux_ax1 , ax2, aux_ax2

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

    #ax1.axis["top"].set_visible(False)

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

    #ax2.axis["top"].set_visible(False)

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

def project_to_2d(map_root_list, integrate_axis=2, get_box=False):

    if len(map_root_list) == 1:
        real_map = algebra.make_vect(algebra.load(map_root_list[0]))
    else:
        real_map = average_map(map_root_list)
    
    if get_box:
        real_box, real_box_info = gridding.physical_grid_largeangle(real_map, 
                                                                    refinement=1)
    
    real_map_2d = np.ma.array(np.sum(real_map, axis=integrate_axis))
    real_map_2d[real_map_2d==0] = np.ma.masked
    coor = np.zeros((2, real_map_2d.shape[0]+1, real_map_2d.shape[1]+1))
    coor_index = 0
    for axis_index in range(real_map.ndim):
        if axis_index == integrate_axis:
            continue
        axis_name = real_map.axes[axis_index]
        axis_cent = real_map.get_axis(axis_name)
        axis_delt = real_map.info['%s_delta'%axis_name]
        axis_edge = axis_cent - 0.5*axis_delt
        axis_edge = np.append(axis_edge, axis_edge[-1]+axis_delt)
        if axis_name == 'freq':
            axis_edge = 1.42e9/axis_edge - 1.
        if coor_index == 0:
            coor[coor_index,...] = axis_edge[:, None]
            coor_index += 1
        else:
            coor[coor_index,...] = axis_edge[None, :]

    coor_index_array = np.arange(np.product(coor.shape[1:])).reshape(coor.shape[1:])
    coor_tri_up = np.zeros(real_map_2d.shape + (3,))
    coor_tri_up[...,0] = coor_index_array[:-1, :-1]
    coor_tri_up[...,1] = coor_index_array[1: , :-1]
    coor_tri_up[...,2] = coor_index_array[:-1, 1: ]
    coor_tri_lo = np.zeros(real_map_2d.shape + (3,))
    coor_tri_lo[...,0] = coor_index_array[1: , :-1]
    coor_tri_lo[...,1] = coor_index_array[1: , 1: ]
    coor_tri_lo[...,2] = coor_index_array[:-1, 1: ]

    coor_tri_2d = np.append(coor_tri_up, coor_tri_lo, axis=0)
    real_map_2d = np.append(real_map_2d, real_map_2d, axis=0)

    if get_box:
        return real_map_2d, coor, coor_tri_2d, real_box
    else:
        return real_map_2d, coor, coor_tri_2d

def plot_sky_ra_dec(cat_root, map_root_list, save_name='./png/ra_dec.png'):

    real_cat = np.loadtxt(cat_root)

    real_map ,real_coor, triang= project_to_2d(map_root_list, integrate_axis=0)

    fig = plt.figure(figsize=(20,20))
    fig.clf()
    ax1, aux_ax1, ax2, aux_ax2 = setup_axes_ra_dec(fig, real_coor[0].min(),
                                                        real_coor[0].max(),
                                                        real_coor[1].min(),
                                                        real_coor[1].max())

    aux_ax1.scatter(real_cat[:,0]*180./np.pi, real_cat[:,1]*180./np.pi, s=4, 
                   edgecolor='none', alpha=0.5)

    real_coor = real_coor.reshape(2, -1)
    triang = triang.reshape(-1,3).astype('int')
    real_map = real_map.flatten()
    cmap = plt.get_cmap('Greens')
    #norm = plt.normalize(real_map.flatten().min(), real_map.flatten().max())
    im = aux_ax2.tripcolor(real_coor[0], real_coor[1], triang, 
                           facecolors=real_map, edgecolors='none',
                           mask= real_map==0,
                           #cmap=cmap,
                           )
    #plt.colorbar(im, orientation='horizontal')

    plt.savefig(save_name)

def plot_sky(cat_root, map_root_list, save_name='./png/real_real_cat.png'):

    real_cat = np.loadtxt(cat_root)
    
    real_map ,real_coor, triang, real_box = project_to_2d(map_root_list, 
                                                          integrate_axis=2,
                                                          get_box=True)
    
    #fig = plt.figure(figsize=(15,10))
    fig = plt.figure(figsize=(26,8))
    fig.clf()
    
    ax1, aux_ax1, ax2, aux_ax2, ax3 = setup_axes(fig, real_coor[1].min(),
                                                      real_coor[1].max(),
                                                      real_coor[0].min(),
                                                      real_coor[0].max())
        #setup_axes(fig, -3*15, 4*15, 0.01, 0.3)
    
    aux_ax1.scatter(real_cat[:,0]*180./np.pi, real_cat[:,2], s=4, 
                    edgecolor='none', alpha=0.5)
    
    #aux_ax2.scatter(ras, decs, s=4, edgecolor='none', alpha=0.0025)
    real_coor = real_coor.reshape(2, -1)
    triang = triang.reshape(-1,3).astype('int')
    real_map = real_map.flatten()
    cmap = plt.get_cmap('Greens')
    norm = plt.normalize(real_map.flatten().min(), real_map.flatten().max())
    im = aux_ax2.tripcolor(real_coor[1], real_coor[0], triang, 
                           facecolors=real_map, edgecolors='none',
                           #mask= real_map==0,
                           #cmap=cmap,
                           )
    
    real_box = np.ma.sum(real_box, axis=2)
    real_box = np.ma.array(real_box)
    real_box[real_box==0] = np.ma.masked
    im = ax3.pcolormesh(real_box, )
    #plt.colorbar(im)
    
    plt.savefig(save_name)

if __name__=="__main__":    
    
    
    
    root_cat = '/mnt/scratch-gl/ycli/2df_catalog/catalog/'

    root_map = '/mnt/raid-project/gmrt/anderson/first_parkes_pipe/maps/'

    #---------------------------------------------------------------------

    name_cat = 'real_catalogue_2df'
    name_map = 'test_allbeams_27n30_10by7_clean_map_I_1315'

    #plot_sky_ra_dec(root_cat + name_cat + '.out',
    #                [root_map + name_map + '.npy',], 
    #                save_name='./png/parkes_cat_ra_dec.png')

    #plot_sky(root_cat + name_cat + '.out', 
    #         [root_map + name_map + '.npy', ], 
    #         save_name='./png/parkes_cat_ra_z.png')

    #exit()
    #---------------------------------------------------------------------

    root_map = '/mnt/scratch-gl/ycli/2df_catalog/map/map_2929.5_full_selection_1000mock/'

    name_cat = 'real_catalogue_2df'
    name_map = 'real_map_2df'

    #plot_sky_ra_dec(root_cat + name_cat + '.out',
    #                [root_map + name_map + '.npy',], 
    #                save_name='./png/real_cat_ra_dec.png')

    plot_sky(root_cat + name_cat + '.out', 
             [root_map + name_map + '.npy', ], 
             save_name='./png/real_cat.png')

    exit()
    #---------------------------------------------------------------------


    #name_cat = 'mock_catalogue_2df_090'
    #name_map = 'mock_map_2df_delta_090'

    #plot_sky(root_cat + name_cat + '.out', 
    #         [root_map + name_map + '.npy',], 
    #         save_name='./png/mock_cat_090.png')

    #name_cat = 'mock_catalogue_2df_090'
    #map = []
    #for i in range(100):
    #    map.append(root_map + 'mock_map_2df_delta_%03d.npy'%i)

    #plot_sky(root_cat + name_cat + '.out', 
    #         map, 
    #         save_name='./png/mock_cat_average.png')


    name_cat = 'mock_catalogue_2df_0090'
    name_map = 'sele_map_2df'

    plot_sky_ra_dec(root_cat + name_cat + '.out',
                    [root_map + name_map + '.npy',], 
                    save_name='./png/sele_cat_ra_dec.png')

    #plot_sky(root_cat + name_cat + '.out', 
    #         [root_map + name_map + '.npy',], 
    #         save_name='./png/sele_cat.png')

    name_cat = 'mock_catalogue_2df_0090'
    name_map = 'sele_map_2df_separable'

    #plot_sky(root_cat + name_cat + '.out', 
    #         [root_map + name_map + '.npy',], 
    #         save_name='./png/sele_cat_separable.png')




