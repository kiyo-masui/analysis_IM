"""tools for plotting data cubes; this code is new and in developement"""
import subprocess
import os
import numpy as np
import numpy.ma as ma
from core import algebra
import multiprocessing
from utils import data_paths
from plotting import plot_slice
import map.beam as beam

cube_frame_dir = "/scratch/eswitzer/cube_frames/"


def gnuplot_radec_single(plotitem):
    """plot a single map slice; helper function for process pool"""
    (index, cube_slice, xaxis, yaxis, vaxis, xylabels, aspect, \
            title, cbar_title, fileprefix, freq, physical) = plotitem

    outfilename = fileprefix + str('.%03d' % index) + '.jpeg'

    plot_slice.gnuplot_radec_slice(outfilename, cube_slice, xaxis, yaxis,
                            vaxis, xylabels, aspect, title, cbar_title,
                            freq, index, physical=physical)


def make_cube_movie(source_key, colorbar_title, frame_dir,
                    filetag_suffix="",
                    outputdir="./", sigmarange=6., ignore=None, multiplier=1.,
                    transverse=False, title=None, sigmacut=None,
                    logscale=False, physical=False, convolve=False, tag=None):
    """Make a stack of spatial slice maps and animate them
    transverse plots along RA and freq and image plane is in Dec
    First mask any points that exceed `sigmacut`, and then report the extent of
    `sigmarange` away from the mean
    """
    datapath_db = data_paths.DataPath()

    if tag is None:
        tag = '_'.join(source_key.split(";"))
        tag = '-'.join(tag.split(":"))

        # for a given path
        #tag = ".".join(source_key.split(".")[:-1])  # extract root name
        #tag = tag.split("/")[-1]

    print tag
    fileprefix = frame_dir + tag
    nlevels = 500

    if transverse:
        orientation = "_freqRA"
    else:
        orientation = "_RADec"

    if not title:
        title = tag

    # prepare the data
    #cube = algebra.make_vect(algebra.load(source_key)) * multiplier
    cube =  datapath_db.fetch_multi(source_key) * multiplier
    if logscale:
        cube = np.log10(cube)

    isnan = np.isnan(cube)
    isinf = np.isinf(cube)
    maskarray = ma.mask_or(isnan, isinf)

    if ignore is not None:
        maskarray = ma.mask_or(maskarray, (cube == ignore))

    convolved = ""
    if convolve:
        convolved = "_convolved"
        beam_data = np.array([0.316148488246, 0.306805630985, 0.293729620792,
                     0.281176247549, 0.270856788455, 0.26745856078,
                     0.258910010848, 0.249188429031])
        freq_data = np.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
        freq_data *= 1.0e6

        beamobj = beam.GaussianBeam(beam_data, freq_data)
        cube = beamobj.apply(cube)

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
    outfile = outputdir + tag + filetag_suffix + orientation + convolved + '.mp4'
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


def plot_gbt_maps(keyname, transverse=False,
                  make_noise_diag=False, make_map=True,
                  make_dirty_map=False, make_noise_inv=False,
                  outputdir="./", sigmarange=[0., 0.001]):
    r"""plot the 15hr, 22hr and 1hr real maps"""
    #datapath_db = data_paths.DataPath()

    section_list = ['A', 'B', 'C', 'D']
    for section in section_list:
        if make_map:
            title = "Sec. %s, %s" % (section, keyname)
            make_cube_movie("db:%s:%s;clean_map" % (keyname, section),
                               "Temperature (mK)", cube_frame_dir,
                               sigmarange=3.,
                               outputdir=outputdir, multiplier=1000.,
                               transverse=transverse,
                               title=title)

        if make_dirty_map:
            title = "Sec. %s, %s (dirty)" % (section, keyname)
            make_cube_movie("db:%s:%s;dirty_map" % (keyname, section),
                               "Temperature (mK)", cube_frame_dir,
                               sigmarange=3.,
                               outputdir=outputdir, multiplier=1000.,
                               transverse=transverse,
                               title=title)

        if make_noise_diag:
            title = "Sec. %s, %s (noise diag)" % (section, keyname)
            make_cube_movie("db:%s:%s;noise_diag" % (keyname, section),
                               "Covariance", cube_frame_dir,
                               sigmarange=sigmarange,
                               outputdir=outputdir, multiplier=1.,
                               logscale=False,
                               transverse=transverse,
                               title=title)

        if make_noise_inv:
            title = "Sec. %s, %s (noise inv)" % (section, keyname)
            make_cube_movie("db:%s:%s;noise_inv" % (keyname, section),
                               "Covariance inverse", cube_frame_dir,
                               sigmarange=-1,
                               outputdir=outputdir, multiplier=1.,
                               logscale=False,
                               transverse=transverse,
                               title=title)


def plot_cleaned_maps(source_key, alt_weight=None,
                 signal='map', weight='noise_inv', divider_token=";",
                 outputdir="/cita/d/www/home/eswitzer/movies/",
                 convolve=False,
                 transverse=False):
    r"""
    `source_key` is the file db key for the maps to combine
    `signal` is the tag in the file db entry for the signal maps
    `weight` is the tag in the file db entry for the N^-1 weights
    `divider_token` is the token that divides the map section name
            from the data type e.g. "A_with_B;noise_inv"
    """
    datapath_db = data_paths.DataPath()
    source_fdb = datapath_db.fetch(source_key, intend_read=True,
                                   silent=True)
    source_fdict = source_fdb[1]

    # accumulate all the files to combine
    weightkeys = {}
    signalkeys = {}
    for filekey in source_fdb[0]:
        if divider_token in filekey:
            data_type = filekey.split(divider_token)[1]
            map_section = filekey.split(divider_token)[0]

            if data_type == signal:
                signalkeys[map_section] = source_fdict[filekey]

            if data_type == weight:
                weightkeys[map_section] = source_fdict[filekey]

    for mapkey in signalkeys:
        signalfile = signalkeys[mapkey]
        weightfile = weightkeys[mapkey]
        print "loading pair: %s %s" % (signalfile, weightfile)

        make_cube_movie(signalfile, "Temperature (mK)", cube_frame_dir,
                        sigmarange=2.5, outputdir=outputdir, multiplier=1000.,
                        transverse=transverse, convolve=convolve)

        make_cube_movie(weightfile, "Inverse variance weight", cube_frame_dir,
                        sigmarange=2.5, outputdir=outputdir, multiplier=1.,
                        transverse=transverse)

        #algebra.save(signal_out, newmap)
        #make_cube_movie(filename, "Cleaned map times weights", cube_frame_dir,
        #                sigmarange=-1, outputdir=outputdir, multiplier=1000.,
        #                transverse=transverse, convolve=convolve)


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
