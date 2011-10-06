"""tools for plotting data cubes; this code is new and in developement"""
import subprocess
import sys
import numpy as np
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import algebra
import multiprocessing


def plot_single(plotitem):
    """plot a single map slice; helper function for process pool"""
    (index, cube_slice, xaxis, yaxis, zaxis, \
            title, cbar_title, fileprefix) = plotitem
    print title
    plt.figure(figsize=(7,3.3))
    cplot = plt.contourf(xaxis, yaxis, np.transpose(cube_slice), zaxis)
    plt.axis('scaled')
    plt.xlim((np.min(xaxis), np.max(xaxis)))
    plt.ylim((np.min(yaxis), np.max(yaxis)))
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title(title, fontsize=9)
    cticks = np.linspace(min(zaxis), max(zaxis), 7, endpoint=True)
    cbar = plt.colorbar(cplot, ticks=cticks)
    cbar.ax.set_yticklabels([('%3.1f' % val) for val in cticks])
    #cbar = plt.colorbar(cplot)
    cbar.ax.set_ylabel(cbar_title)
    filename = fileprefix + str('%03d' % index) + '.png'
    plt.savefig(filename, dpi=200)
    plt.clf()


# TODO: delete frames once generated
def make_cube_movie(cubename, colorbar_title, cube_frame_dir,
                    outputdir="", sigmarange=3., ignore=None, multiplier=1.):
    """Make a stack of spatial slice maps and animate them"""
    cube = algebra.make_vect(algebra.load(cubename)) * multiplier
    tag = ".".join(cubename.split(".")[:-1])  # extract root name
    fileprefix = cube_frame_dir + tag

    if ignore:
        cube[cube == ignore] = ma.masked
    cube_mean = ma.mean(cube)
    cube_std = ma.std(cube)
    try:
        len(sigmarange)
        coloraxis = np.linspace(sigmarange[0], sigmarange[1],
                                500, endpoint=True)
    except TypeError:
        if (sigmarange > 0.):
            coloraxis = np.linspace(cube_mean - sigmarange * cube_std,
                                    cube_mean + sigmarange * cube_std,
                                    100, endpoint=True)
        else:
            coloraxis = np.linspace(np.min(cube),  np.max(cube),
                                    100, endpoint=True)

    freq_axis = cube.get_axis('freq')
    space_axis = (cube.get_axis('ra'), cube.get_axis('dec'))

    runlist = []
    for freqind in range(cube.shape[0]):
        fulltitle = tag + " (freq = %3.1f MHz)" % (freq_axis[freqind] / 1.e6)
        runlist.append((freqind, cube[freqind, :, :], space_axis[0],
                        space_axis[1], coloraxis, fulltitle,
                        colorbar_title, fileprefix))

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(plot_single, runlist)

    subprocess.check_call(('ffmpeg', '-r', '10', '-y', '-i', cube_frame_dir + tag +\
               '%03d.png', output_dir + tag + '.mp4'))


