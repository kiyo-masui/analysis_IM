"""tools for plotting data cubes"""
import subprocess
import os
import sys
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from core import algebra
import multiprocessing
matplotlib.use('Agg')

print 'Executing on', os.uname()
print 'Python version', sys.version
print 'matplotlib version', matplotlib.__version__

not_found_msg = """
The mencoder command was not found;
mencoder is used by this script to make an avi file from a set of pngs.
It is typically not installed by default on linux distros because of
legal restrictions, but it is widely available.
"""

try:
    subprocess.check_call(['mencoder'])
except subprocess.CalledProcessError:
    print "mencoder command was found"
    pass # mencoder is found, but returns non-zero exit as expected
    # This is a quick and dirty check; it leaves some spurious output
    # for the user to puzzle over.
except OSError:
    print not_found_msg
    sys.exit("quitting\n")

def make_cube_movie(filename, title, colorbar_title):
    cube = algebra.make_vect(algebra.load(filename))
    cube_min = np.min(cube)
    cube_max = np.max(cube)
    coloraxis = np.linspace(cube_min, cube_max, 100, endpoint=True)

    nfreq = cube.shape[0]
    print nfreq, cube_min, cube_max
    freq_axis = cube.get_axis('freq')
    ra_axis = cube.get_axis('ra')
    dec_axis = cube.get_axis('dec')
    #for freqind in range(nfreq):
    for freqind in range(5):
        cube_slice = cube[freqind, :, :]
        f = plt.contourf(ra_axis, dec_axis, np.transpose(cube_slice), coloraxis)
        freq = freq_axis[freqind]
        plt.axis('scaled')
        plt.xlim((np.min(ra_axis), np.max(ra_axis)))
        plt.ylim((np.min(dec_axis), np.max(dec_axis)))
        plt.xlabel("RA")
        plt.ylabel("Dec")
        fulltitle = title + " (freq = "+repr(freq/1e6)+" MHz)"
        print fulltitle
        plt.title(fulltitle)
        ##c = plt.colorbar(f, ticks=coloraxis)
        c = plt.colorbar(f)

        c.ax.set_ylabel(colorbar_title)
        filename = str('%03d' % freqind) + '.png'
        plt.savefig(filename, dpi=100)
        plt.clf()

make_cube_movie("data/reg15selection.npy", "Selection function", 
                "# galaxies/pixel")

#command = ('mencoder mf://*.png -mf type=png:w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -of mpeg1video -o output.mpg')

#os.spawnvp(os.P_WAIT, 'mencoder', command)
#subprocess.check_call(command)


