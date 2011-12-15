"""tools for plotting data cubes; this code is new and in developement"""
import subprocess
import os
import numpy as np
import numpy.ma as ma
from core import algebra
import multiprocessing
from utils import data_paths
from plotting import plot_slice

cube_frame_dir = "/scratch/eswitzer/cube_frames/"


def gnuplot_radec_single(plotitem):
    """plot a single map slice; helper function for process pool"""
    (index, cube_slice, xaxis, yaxis, vaxis, xylabels, aspect, \
            title, cbar_title, fileprefix, freq, physical) = plotitem

    outfilename = fileprefix + str('.%03d' % index) + '.jpeg'

    plot_slice.gnuplot_radec_slice(outfilename, cube_slice, xaxis, yaxis,
                            vaxis, xylabels, aspect, title, cbar_title,
                            freq, index, physical=physical)


def make_cube_movie(cubename, colorbar_title, frame_dir,
                    filetag_suffix="",
                    outputdir="./", sigmarange=6., ignore=None, multiplier=1.,
                    transverse=False, title=None, sigmacut=None,
                    logscale=False, physical=False):
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
    (ra_axis, dec_axis) = (cube.get_axis('ra'), cube.get_axis('dec'))

    runlist = []
    # TODO: make transverse work with gnuplot
    if transverse:
        for decind in range(cube.shape[2]):
            fulltitle = title + " (dec = %3.1f)" % (dec_axis[decind])
            runlist.append((decind, cube[:, :, decind], freq_axis,
                            ra_axis, color_axis, ["Freq", "Ra"], 20.,
                            fulltitle, colorbar_title, fileprefix, 800.))
    else:
        for freqind in range(cube.shape[0]):
            if physical:
                runlist.append((freqind, cube[freqind, :, :], ra_axis,
                                dec_axis, color_axis,
                                ["x (RA, cMpc/h)", "y (Dec, cMpc/h)"], 1.,
                                title, colorbar_title, fileprefix,
                                freq_axis[freqind], physical))
            else:
                # convert to MHz
                runlist.append((freqind, cube[freqind, :, :], ra_axis,
                                dec_axis, color_axis, ["RA", "Dec"], 1.,
                                title, colorbar_title, fileprefix,
                                freq_axis[freqind] / 1.e6, physical))

    pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count() - 4))
    #pool.map(plot_single, runlist)
    pool.map(gnuplot_radec_single, runlist)
    #gnuplot_radec_single(runlist[0])  # for troubleshooting

    #argument = frame_dir + tag + '.%03d.png'
    argument = frame_dir + tag + '.%03d.jpeg'
    outfile = outputdir + tag + filetag_suffix + orientation + '.mp4'
    subprocess.check_call(('ffmpeg', '-vb', '2500000', '-r', '10',
                           '-y', '-i',
                            argument, outfile))

    for fileindex in range(len(runlist)):
        #os.remove(fileprefix + str('.%03d' % fileindex) + '.png')
        os.remove(fileprefix + str('.%03d' % fileindex) + '.jpeg')


def make_cube_slice(cubename, outfilename, slice_index, colorbar_title,
                    sigmarange=6., ignore=None, multiplier=1.,
                    title=None, sigmacut=None,
                    logscale=False):
    """Make a stack of spatial slice maps and animate them
    transverse plots along RA and freq and image plane is in Dec
    First mask any points that exceed `sigmacut`, and then report the extent of
    `sigmarange` away from the mean
    """
    # set up the labels:
    tag = ".".join(cubename.split(".")[:-1])  # extract root name
    tag = tag.split("/")[-1]
    nlevels = 500

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

    plot_slice.gnuplot_radec_slice("placeholder.jpeg", cube[slice_index, :, :],
                            ra_axis, dec_axis, color_axis,
                            ["RA", "Dec"], 1., title, colorbar_title,
                            freq_axis[slice_index], slice_index,
                            eps_outfile=outfilename)


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
