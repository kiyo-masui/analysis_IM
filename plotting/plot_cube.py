"""tools for plotting data cubes; this code is new and in developement"""
import subprocess
import os
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import numpy.ma as ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import algebra
import multiprocessing
from utils import data_paths
import tempfile
import subprocess
import sys
from utils.cosmology import Cosmology
from utils import units
from core import constants as cc

#cube_frame_dir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/cube_frames/"
cube_frame_dir = "/scratch/eswitzer/cube_frames/"

def gnuplot_single(plotitem):
    """plot a single map slice; helper function for process pool"""
    (index, cube_slice, xaxis, yaxis, vaxis, xylabels, aspect, \
            title, cbar_title, fileprefix, freq) = plotitem
    print title, repr(cube_slice.shape)

    # set up the size calc
    cosmology = Cosmology()
    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    beam_FWHM = interp1d(freq_data, beam_data)
    FWHM = beam_FWHM(800.)
    z_here = cc.freq_21cm_MHz/freq
    angular_scale = 20. / units.degree / cosmology.proper_distance(z_here)
    FWHM = beam_FWHM(freq)


    input_data_file = tempfile.NamedTemporaryFile()
    #input_data_file = open("/tmp/plotdata.out", "w")
    gplfile = tempfile.NamedTemporaryFile(suffix=".gpl")
    outplot_file = tempfile.NamedTemporaryFile(suffix=".eps")

    xminmax = (np.min(xaxis), np.max(xaxis))
    yminmax = (np.min(yaxis), np.max(yaxis))
    cminmax = (np.min(vaxis), np.max(vaxis))
    aspect = (yminmax[1] - yminmax[0]) / (xminmax[1] - xminmax[0])
    print FWHM, angular_scale

    # TODO: optimize this
    for xind in range(len(xaxis)):
        for yind in range(len(yaxis)):
            outstring = "%g %g %g\n" % \
                        (xaxis[xind], yaxis[yind], cube_slice[xind, yind])
            input_data_file.write(outstring)

    input_data_file.flush()

    # TODO: set uniform colobar range
    gplfile.write("set view map\n")
    gplfile.write("set xrange [%g:%g]\n" % xminmax)
    gplfile.write("set yrange [%g:%g]\n" % yminmax)
    gplfile.write("set cbrange [%g:%g]\n" % cminmax)
    #gplfile.write("set size %g, %g\n" % (0.9, 0.9*aspect))
    gplfile.write("set size ratio -1\n")
    #gplfile.write("set size square\n")
    gplfile.write("set nokey\n")
    gplfile.write("set tics out\n")
    gplfile.write('set xlabel "%s"\n' % xylabels[0])
    gplfile.write('set ylabel "%s"\n' % xylabels[1])
    gplfile.write('set title "%s"\n' % title)
    gplfile.write('set cblabel "%s"\n' % cbar_title)
    gplfile.write('set rmargin 10\n')
    #gplfile.write("set terminal postscript eps enhanced color\n")
    gplfile.write("set terminal postscript eps color size %g, %g\n" % (5, 5.*aspect*1.1))
    gplfile.write('set output "%s"\n' % outplot_file.name)

    gplfile.write('set obj 10 circle at graph 0.9,.1 size %g front\n' % \
                    (FWHM/2.))
    gplfile.write('set obj 10 lw 5 fill empty border rgb "purple"\n')
    gplfile.write('set obj 11 rect at graph 0.9,.1 size %g, %g front\n' % \
                   (angular_scale, angular_scale))
    gplfile.write('set obj 11 lw 3 fill empty border rgb "black"\n')
    #gplfile.write('set obj 10 fc rgb "blue"\n')

    #gplfile.write("set palette positive nops_allcF maxcolors 0 gamma 1.5 gray\n")
    #gplfile.write("set palette color\n")
    gplfile.write("set palette rgbformulae 22, 13, -31\n")

    #gplfile.write("set pm3d map\n")
    #gplfile.write("set pm3d interpolate 5,5\n")
    #gplfile.write('splot "%s" matrix using 1:2:3 with pm3d\n' % input_data_file.name)
    gplfile.write('plot "%s" using 1:2:3 with image\n' % input_data_file.name)
    gplfile.flush()

    subprocess.check_call(('/cita/h/home-2/eswitzer/local/bin/gnuplot', gplfile.name))
    gplfile.close()
    input_data_file.close()

    #filename = fileprefix + str('.%03d' % index) + '.png'
    filename = fileprefix + str('.%03d' % index) + '.jpeg'
    # consider adding -flatten to avoid transparency
    subprocess.check_call(('convert', '-density', '300', '-trim', '+repage',
                           '-border', '40x40', '-bordercolor', 'white',
                            outplot_file.name, filename))
    outplot_file.close()

def plot_single(plotitem):
    """plot a single map slice; helper function for process pool"""
    (index, cube_slice, xaxis, yaxis, vaxis, xylabels, aspect, \
            title, cbar_title, fileprefix, freq) = plotitem
    print title, repr(cube_slice.shape)
    plt.figure(figsize=(7, 3.3))
    cplot = plt.contourf(xaxis, yaxis, np.transpose(cube_slice), vaxis)
    plt.axis('scaled')
    plt.axes().set_aspect(aspect)  # 'equal' also?
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


def make_cube_movie(cubename, colorbar_title, frame_dir,
                    filetag_suffix="",
                    outputdir="./", sigmarange=6., ignore=None, multiplier=1.,
                    transverse=False, title=None, sigmacut=None,
                    logscale=False):
    """Make a stack of spatial slice maps and animate them
    transverse plots along RA and freq and image plane is in Dec
    First mask any points that exceed `sigmacut`, and then report the extent of
    `sigmarange` away from the mean
    """
    # set up the labels:
    tag = ".".join(cubename.split(".")[:-1])  # extract root name
    tag = tag.split("/")[-1]
    fileprefix = frame_dir + tag
    nlevels = 500

    if transverse:
        orientation = "_freqRA"
    else:
        orientation = "_RADec"

    if not title:
        title = tag

    # prepare the data
    cube = algebra.make_vect(algebra.load(cubename)) * multiplier
    if logscale:
        cube = np.log10(cube)

    isnan = np.isnan(cube)
    isinf = np.isinf(cube)
    maskarray = ma.mask_or(isnan, isinf)

    if ignore is not None:
        maskarray = ma.mask_or(maskarray, (cube == ignore))

    if sigmacut:
        #np.set_printoptions(threshold=np.nan, precision=4)
        deviation = np.abs((cube - np.mean(cube)) / np.std(cube))
        extend_maskarray = (cube > (sigmacut * deviation))
        maskarray = ma.mask_or(extend_maskarray, maskarray)

    mcube = ma.masked_array(cube, mask=maskarray)

    try:
        whmaskarray = np.where(maskarray)[0]
        mask_fraction = float(len(whmaskarray)) / float(cube.size)
    except:
        mask_fraction = 0.

    print "fraction of map clipped: %f" % mask_fraction
    (cube_mean, cube_std) = (mcube.mean(), mcube.std())
    print "cube mean=%g std=%g" % (cube_mean, cube_std)

    try:
        len(sigmarange)
        color_axis = np.linspace(sigmarange[0], sigmarange[1],
                                nlevels, endpoint=True)
    except TypeError:
        if (sigmarange > 0.):
            color_axis = np.linspace(cube_mean - sigmarange * cube_std,
                                    cube_mean + sigmarange * cube_std,
                                    nlevels, endpoint=True)
        else:
            color_axis = np.linspace(ma.min(mcube),  ma.max(mcube),
                                    nlevels, endpoint=True)

    print "using range: [%g, %g]" % (np.min(color_axis), np.max(color_axis))

    freq_axis = cube.get_axis('freq')
    freq_axis /= 1.e6  # convert to MHz
    (ra_axis, dec_axis) = (cube.get_axis('ra'), cube.get_axis('dec'))

    runlist = []
    if transverse:
        for decind in range(cube.shape[2]):
            fulltitle = title + " (dec = %3.1f)" % (dec_axis[decind])
            runlist.append((decind, cube[:, :, decind], freq_axis,
                            ra_axis, color_axis, ["Freq", "Ra"], 20.,
                            fulltitle, colorbar_title, fileprefix, 800.))
    else:
        for freqind in range(cube.shape[0]):
            fulltitle = title + " (freq = %3.1f MHz)" % (freq_axis[freqind])
            runlist.append((freqind, cube[freqind, :, :], ra_axis,
                            dec_axis, color_axis, ["RA", "Dec"], 1.,
                            fulltitle, colorbar_title, fileprefix,
                            freq_axis[freqind]))

    pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count() - 4))
    #pool.map(plot_single, runlist)
    pool.map(gnuplot_single, runlist)
    #gnuplot_single(runlist[0])

    #argument = frame_dir + tag + '.%03d.png'
    argument = frame_dir + tag + '.%03d.jpeg'
    outfile = outputdir + tag + filetag_suffix + orientation + '.mp4'
    subprocess.check_call(('ffmpeg', '-r', '10', '-y', '-i',
                            argument, outfile))

    for fileindex in range(len(runlist)):
        #os.remove(fileprefix + str('.%03d' % fileindex) + '.png')
        os.remove(fileprefix + str('.%03d' % fileindex) + '.jpeg')


def plot_gbt_maps(keyname, transverse=False, skip_noise=False, skip_map=False,
                  outputdir="./", sigmarange=[0., 0.001]):
    r"""plot the 15hr, 22hr and 1hr real maps"""
    datapath_db = data_paths.DataPath()

    section_list = ['A', 'B', 'C', 'D']
    for section in section_list:
        if not skip_map:
            filename = datapath_db.fetch(keyname, intend_read=True,
                                         pick=(section + ';clean_map'))
            title = "Sec. %s, %s" % (section, keyname)
            make_cube_movie(filename,
                               "Temperature (mK)", cube_frame_dir,
                               sigmarange=6.,
                               outputdir=outputdir, multiplier=1000.,
                               transverse=transverse,
                               title=title)

            filename = datapath_db.fetch(keyname, intend_read=True,
                                         pick=(section + ';dirty_map'))
            title = "Sec. %s, %s (dirty)" % (section, keyname)
            make_cube_movie(filename,
                               "Temperature (mK)", cube_frame_dir,
                               sigmarange=6.,
                               outputdir=outputdir, multiplier=1000.,
                               transverse=transverse,
                               title=title)

        if not skip_noise:
            filename = datapath_db.fetch(keyname, intend_read=True,
                                         pick=(section + ';noise_diag'))
            title = "Sec. %s, %s (noise)" % (section, keyname)
            # sigmacut=0.008
            make_cube_movie(filename,
                               "Covariance", cube_frame_dir,
                               sigmarange=sigmarange,
                               outputdir=outputdir, multiplier=1.,
                               logscale=False,
                               transverse=transverse,
                               title=title)


def plot_simulations(keyname, transverse=False, outputdir="./"):
    """make movies of the 15hr simulations
    permutations: with or without streaming, including beam, adding real data
    """
    datapath_db = data_paths.DataPath()

    filename = datapath_db.fetch(keyname, pick='1')
    make_cube_movie(filename, "Temperature (mK)", cube_frame_dir,
                    sigmarange=5., outputdir=outputdir, multiplier=1000.,
                    transverse=transverse)
    #make_cube_movie(root_directory + "sim_beam_" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=5.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #make_cube_movie(root_directory + "sim_beam_plus_data" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=10.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #make_cube_movie(root_directory + "simvel_" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=5.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #make_cube_movie(root_directory + "simvel_beam_" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=5.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #make_cube_movie(root_directory + "simvel_beam_plus_data" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=10.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)


def plot_difference(filename1, filename2, title, sigmarange=6., sigmacut=None,
                    transverse=False, outputdir="./", multiplier=1000.,
                    logscale=False, fractional=False,
                    ignore=None, diff_filename="./difference.npy"):
    """make movies of the difference of two maps (assuming same dimensions)"""
    map1 = algebra.make_vect(algebra.load(filename1))
    map2 = algebra.make_vect(algebra.load(filename2))

    if fractional:
        difftitle = "fractional diff."
        dmap = (map1 - map2) / map1 * 100.
    else:
        difftitle = "difference"
        dmap = map1 - map2

    algebra.save(diff_filename, dmap)

    make_cube_movie(diff_filename,
                       difftitle, cube_frame_dir, sigmarange=6.,
                       sigmacut=sigmacut, outputdir=outputdir, ignore=ignore,
                       multiplier=multiplier, transverse=transverse,
                       logscale=False)

    make_cube_movie(filename1,
                       title, cube_frame_dir, sigmarange=sigmarange,
                       sigmacut=sigmacut, outputdir=outputdir, ignore=ignore,
                       multiplier=multiplier, transverse=transverse,
                       logscale=logscale, filetag_suffix="_1")

    make_cube_movie(filename2,
                       title, cube_frame_dir, sigmarange=sigmarange,
                       sigmacut=sigmacut, outputdir=outputdir, ignore=ignore,
                       multiplier=multiplier, transverse=transverse,
                       logscale=logscale, filetag_suffix="_2")
