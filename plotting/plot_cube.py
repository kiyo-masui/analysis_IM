"""tools for plotting data cubes; this code is new and in developement"""
import subprocess
import os
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import algebra
import multiprocessing


def plot_single(plotitem):
    """plot a single map slice; helper function for process pool"""
    (index, cube_slice, xaxis, yaxis, vaxis, xylabels, aspect, \
            title, cbar_title, fileprefix) = plotitem
    print title, repr(cube_slice.shape)
    plt.figure(figsize=(7,3.3))
    cplot = plt.contourf(xaxis, yaxis, np.transpose(cube_slice), vaxis)
    plt.axis('scaled')
    plt.axes().set_aspect(aspect) # 'equal' also?
    plt.xlim((np.min(xaxis), np.max(xaxis)))
    plt.ylim((np.min(yaxis), np.max(yaxis)))
    plt.xlabel(xylabels[0])
    plt.ylabel(xylabels[1])
    plt.title(title, fontsize=9)
    cticks = np.linspace(min(vaxis), max(vaxis), 7, endpoint=True)
    cbar = plt.colorbar(cplot, ticks=cticks)
    cbar.ax.set_yticklabels([('%5.2g' % val) for val in cticks])
    #cbar = plt.colorbar(cplot)
    cbar.ax.set_ylabel(cbar_title)
    filename = fileprefix + str('%03d' % index) + '.png'
    plt.savefig(filename, dpi=200)
    plt.clf()


# TODO: delete frames once generated
def make_cube_movie(cubename, colorbar_title, cube_frame_dir,
                    outputdir="./", sigmarange=3., ignore=None, multiplier=1.,
                    transverse=False):
    """Make a stack of spatial slice maps and animate them
    transverse plots along RA and freq and image plane is in Dec
    """
    cube = algebra.make_vect(algebra.load(cubename)) * multiplier
    tag = ".".join(cubename.split(".")[:-1])  # extract root name
    tag = tag.split("/")[-1]
    fileprefix = cube_frame_dir + tag

    if transverse:
        orientation = "_freqRA"
    else:
        orientation = "_RADec"

    if ignore:
        cube[cube == ignore] = ma.masked
    cube_mean = ma.mean(cube)
    cube_std = ma.std(cube)
    try:
        len(sigmarange)
        color_axis = np.linspace(sigmarange[0], sigmarange[1],
                                500, endpoint=True)
    except TypeError:
        if (sigmarange > 0.):
            color_axis = np.linspace(cube_mean - sigmarange * cube_std,
                                    cube_mean + sigmarange * cube_std,
                                    100, endpoint=True)
        else:
            color_axis = np.linspace(np.min(cube),  np.max(cube),
                                    100, endpoint=True)

    freq_axis = cube.get_axis('freq')
    freq_axis /= 1.e6  # convert to MHz
    (ra_axis, dec_axis) = (cube.get_axis('ra'), cube.get_axis('dec'))

    runlist = []
    if transverse:
        for decind in range(cube.shape[2]):
            fulltitle = tag + " (dec = %3.1f)" % (dec_axis[decind])
            runlist.append((decind, cube[:, :, decind], freq_axis,
                            ra_axis, color_axis, ["Freq","Ra"], 20., fulltitle,
                            colorbar_title, fileprefix))
    else:
        for freqind in range(cube.shape[0]):
            fulltitle = tag + " (freq = %3.1f MHz)" % (freq_axis[freqind])
            runlist.append((freqind, cube[freqind, :, :], ra_axis,
                            dec_axis, color_axis, ["RA","Dec"], 1., fulltitle,
                            colorbar_title, fileprefix))

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(plot_single, runlist)

    subprocess.check_call(('ffmpeg', '-r', '10', '-y', '-i', cube_frame_dir + tag +\
               '%03d.png', outputdir + tag + orientation + '.mp4'))

    for fileindex in range(len(runlist)):
        os.remove(fileprefix + str('%03d' % fileindex) + '.png')
