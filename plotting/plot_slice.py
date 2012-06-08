r"""functions to plot slices of the cubes"""
import numpy as np
import tempfile
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def gnuplot_2D(outfilename, region, xaxis, yaxis, vaxis, xylabels,
               aspect, fulltitle, cbar_title, eps_outfile=None,
               draw_objects=[], density='300'):
    r"""plot a 2D matrix
    """

    print fulltitle, repr(region.shape)

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
                        (xaxis[xind], yaxis[yind], region[xind, yind])
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

    objnum = 10
    for drawing in draw_objects:
        drawing["objnum"] = objnum

        objcmd = "set obj %(objnum)d %(primitive)s " % drawing
        objcmd += "at graph %(center_x)g, %(center_y)g " % drawing

        if drawing["primitive"] == "circle":
            objcmd += "size %(radius)g front\n" % drawing

        if drawing["primitive"] == "rect":
            objcmd += "size %(size_x)g, %(size_y)g front\n" % drawing

        gplfile.write(objcmd)

        objcmd = "set obj %(objnum)d lw %(width)d " % drawing
        objcmd += "fill empty border rgb \"%(color)s\"\n" % drawing
        gplfile.write(objcmd)

        objnum += 1

    gplfile.write("set palette rgbformulae 22, 13, -31\n")

    gplfile.write('plot "%s" using 1:2:3 with image\n' % input_data_file.name)
    gplfile.flush()

    gnuplot = '/cita/h/home-2/eswitzer/local/bin/gnuplot'
    subprocess.check_call((gnuplot, gplfile.name))
    gplfile.close()
    input_data_file.close()

    if eps_outfile is None:
        subprocess.check_call(('convert', '-density', density, '-trim', '+repage',
                               '-border', '40x40', '-bordercolor', 'white',
                                outplot_file.name, outfilename))

    outplot_file.close()


def simpleplot_2D(outfilename, region, xaxis, yaxis, xylabels,
               aspect, fulltitle, cbar_title, eps_outfile=None,
               draw_objects=[], density='300', size=5., logscale=False,
               debug_datafile=None):
    r"""similar to above
    TODO: find options that make this call identical
    """

    print fulltitle, repr(region.shape)

    if logscale:
        # TODO: this is a kludge
        region += np.max(region)/1.e6
        try:
            region = np.log10(region)
        except FloatingPointError:
            print "plotting error (all zeros, log requested?)"
            region = np.zeros_like(region)

    if debug_datafile:
        input_data_file = open(debug_datafile, "w")
    else:
        input_data_file = tempfile.NamedTemporaryFile()

    gplfile = tempfile.NamedTemporaryFile(suffix=".gpl")
    outplot_file = tempfile.NamedTemporaryFile(suffix=".eps")

    deltax = abs(xaxis[1] - xaxis[0])
    deltay = abs(yaxis[1] - yaxis[0])
    xminmax = (np.min(xaxis) - 0.5 * deltax, np.max(xaxis) + 0.5 * deltax)
    yminmax = (np.min(yaxis) - 0.5 * deltay, np.max(yaxis) + 0.5 * deltay)
    aspect = (yminmax[1] - yminmax[0]) / (xminmax[1] - xminmax[0])

    # TODO: optimize this
    for xind in range(len(xaxis)):
        for yind in range(len(yaxis)):
            outstring = "%g %g %g\n" % \
                        (xaxis[xind], yaxis[yind], region[xind, yind])
            input_data_file.write(outstring)

    input_data_file.flush()

    gplfile.write("set view map\n")
    gplfile.write("set xrange [%g:%g]\n" % xminmax)
    gplfile.write("set yrange [%g:%g]\n" % yminmax)
    #gplfile.write("set cbrange [%g:%g]\n" % cminmax)
    gplfile.write("set size ratio -1\n")
    gplfile.write("set nokey\n")
    gplfile.write("set tics out\n")
    gplfile.write('set xlabel "%s"\n' % xylabels[0])
    gplfile.write('set ylabel "%s"\n' % xylabels[1])
    gplfile.write('set title "%s"\n' % fulltitle)
    gplfile.write('set cblabel "%s"\n' % cbar_title)
    gplfile.write('set rmargin 10\n')
    gplfile.write("set terminal postscript eps color size %g, %g\n" % \
                  (size, size * aspect * 1.1))

    if eps_outfile is None:
        gplfile.write('set output "%s"\n' % outplot_file.name)
    else:
        gplfile.write('set output "%s"\n' % eps_outfile)

    objnum = 10
    for drawing in draw_objects:
        drawing["objnum"] = objnum

        objcmd = "set obj %(objnum)d %(primitive)s " % drawing
        objcmd += "at graph %(center_x)g, %(center_y)g " % drawing

        if drawing["primitive"] == "circle":
            objcmd += "size %(radius)g front\n" % drawing

        if drawing["primitive"] == "rect":
            objcmd += "size %(size_x)g, %(size_y)g front\n" % drawing

        gplfile.write(objcmd)

        objcmd = "set obj %(objnum)d lw %(width)d " % drawing
        objcmd += "fill empty border rgb \"%(color)s\"\n" % drawing
        gplfile.write(objcmd)

        objnum += 1

    gplfile.write("set palette rgbformulae 22, 13, -31\n")

    if debug_datafile:
        gplfile.write('plot "%s" using 1:2:3 with image\n' % debug_datafile)
    else:
        gplfile.write('plot "%s" using 1:2:3 with image\n' % input_data_file.name)

    gplfile.flush()

    gnuplot = '/cita/h/home-2/eswitzer/local/bin/gnuplot'
    subprocess.check_call((gnuplot, gplfile.name))
    gplfile.close()
    input_data_file.close()

    if eps_outfile is None:
        subprocess.check_call(('convert', '-density', density, '-trim', '+repage',
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
