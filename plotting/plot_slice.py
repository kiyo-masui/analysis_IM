r"""functions to plot slices of the cubes"""
import scipy as sp
import numpy as np
import tempfile
import subprocess
from scipy.interpolate import interp1d
from utils.cosmology import Cosmology
from utils import units
from core import constants as cc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def gnuplot_radec_slice(outfilename, cube_slice, xaxis, yaxis, vaxis, xylabels,
                        aspect, title, cbar_title, slice_label, index,
                        eps_outfile=None, physical=False):
    r"""plot a single map slice
    Here, slice_label is the value of the value over which the slice is fixed
        for physical=False, this is the frequency of the slice in MHz
        for physical=True, this is the comoving distance to the slice
    """

    if not physical:
        z_here = cc.freq_21cm_MHz / slice_label - 1.
        cosmology = Cosmology()
        littleh = (cosmology.H0 / 100.0)
        comoving_distance = cosmology.comoving_distance(z_here) / littleh
        proper_distance = cosmology.proper_distance(z_here) / littleh
        angular_scale = 20. / units.degree / proper_distance

        fulltitle = "%s (i = %d, freq = %3.1f MHz, z = %3.3f, Dc=%3.0f cMpc)" % \
                (title, index, slice_label, z_here, comoving_distance)

        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                     0.281176247549, 0.270856788455, 0.26745856078,
                     0.258910010848, 0.249188429031])
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                                 dtype=float)

        beam_fwhm = interp1d(freq_data, beam_data)
        if (slice_label <= freq_data.min()) or \
           (slice_label >= freq_data.max()):
            fwhm = 0.
        else:
            fwhm = beam_fwhm(slice_label)

    else:
        fulltitle = "%s (i = %d, Dc = %10.3f cMpc)" % \
                (title, index, slice_label)

    print fulltitle, repr(cube_slice.shape)

    input_data_file = tempfile.NamedTemporaryFile()
    gplfile = tempfile.NamedTemporaryFile(suffix=".gpl")
    outplot_file = tempfile.NamedTemporaryFile(suffix=".eps")

    deltax = abs(xaxis[1] - xaxis[0])
    deltay = abs(yaxis[1] - yaxis[0])
    xminmax = (np.min(xaxis) - 0.5 * deltax, np.max(xaxis) + 0.5 * deltax)
    yminmax = (np.min(yaxis) - 0.5 * deltay, np.max(yaxis) + 0.5 * deltay)
    cminmax = (np.min(vaxis), np.max(vaxis))
    aspect = (yminmax[1] - yminmax[0]) / (xminmax[1] - xminmax[0])

    # TODO: optimize this
    for xind in range(len(xaxis)):
        for yind in range(len(yaxis)):
            outstring = "%g %g %g\n" % \
                        (xaxis[xind], yaxis[yind], cube_slice[xind, yind])
            input_data_file.write(outstring)

    input_data_file.flush()

    gplfile.write("set view map\n")
    gplfile.write("set xrange [%g:%g]\n" % xminmax)
    gplfile.write("set yrange [%g:%g]\n" % yminmax)
    gplfile.write("set cbrange [%g:%g]\n" % cminmax)
    gplfile.write("set size ratio -1\n")
    gplfile.write("set nokey\n")
    gplfile.write("set tics out\n")
    gplfile.write('set xlabel "%s"\n' % xylabels[0])
    gplfile.write('set ylabel "%s"\n' % xylabels[1])
    gplfile.write('set title "%s"\n' % fulltitle)
    gplfile.write('set cblabel "%s"\n' % cbar_title)
    gplfile.write('set rmargin 10\n')
    gplfile.write("set terminal postscript eps color size %g, %g\n" % \
                  (5, 5. * aspect * 1.1))

    if eps_outfile is None:
        gplfile.write('set output "%s"\n' % outplot_file.name)
    else:
        gplfile.write('set output "%s"\n' % eps_outfile)

    if not physical:
        gplfile.write('set obj 10 circle at graph 0.9,.15 size %g front\n' % \
                       (fwhm / 2.))
        gplfile.write('set obj 10 lw 5 fill empty border rgb "purple"\n')
        gplfile.write('set obj 11 rect at graph 0.9,.15 size %g, %g front\n' % \
                      (angular_scale, angular_scale))
        gplfile.write('set obj 11 lw 3 fill empty border rgb "black"\n')

    gplfile.write("set palette rgbformulae 22, 13, -31\n")

    gplfile.write('plot "%s" using 1:2:3 with image\n' % input_data_file.name)
    gplfile.flush()

    gnuplot = '/cita/h/home-2/eswitzer/local/bin/gnuplot'
    subprocess.check_call((gnuplot, gplfile.name))
    gplfile.close()
    input_data_file.close()

    if eps_outfile is None:
        subprocess.check_call(('convert', '-density', '300', '-trim', '+repage',
                               '-border', '40x40', '-bordercolor', 'white',
                                outplot_file.name, outfilename))

    outplot_file.close()


def plot_single(plotitem):
    """plot a single map slice; helper function for process pool"""
    (index, cube_slice, xaxis, yaxis, vaxis, xylabels, aspect, \
            title, cbar_title, fileprefix, slice_label) = plotitem
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