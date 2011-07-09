"""tools for plotting data cubes; this code is new and in developement"""
import subprocess
import os
import sys
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import algebra
import multiprocessing

print 'Executing on', os.uname()
print 'Python version', sys.version
print 'matplotlib version', matplotlib.__version__

not_found_msg = """
The mencoder command was not found;
mencoder is used by this script to make an avi file from a set of pngs.
It is typically not installed by default on linux distros because of
legal restrictions, but it is widely available.
"""

#try:
#    subprocess.check_call(['mencoder'])
#except subprocess.CalledProcessError:
#    print "mencoder command was found"
#    pass # mencoder is found, but returns non-zero exit as expected
#except OSError:
#    print not_found_msg
#    sys.exit("quitting\n")

def plot_single(plotitem):
    (index, cube_slice, xaxis, yaxis, zaxis, title, cbar_title, fileprefix) = plotitem
    print title
    f = plt.contourf(xaxis, yaxis, np.transpose(cube_slice), zaxis)
    plt.axis('scaled')
    plt.xlim((np.min(xaxis), np.max(xaxis)))
    plt.ylim((np.min(yaxis), np.max(yaxis)))
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title(title, fontsize=10)
    ##c = plt.colorbar(f, ticks=zaxis)
    c = plt.colorbar(f)
    c.ax.set_ylabel(cbar_title)
    filename = fileprefix + str('%03d' % index) + '.png'
    plt.savefig(filename, dpi=100)
    plt.clf()

def make_cube_movie(filename, tagname, colorbar_title, fileprefix):
    cube = algebra.make_vect(algebra.load(filename))
    cube_min = np.min(cube)
    cube_max = np.max(cube)
    coloraxis = np.linspace(cube_min, cube_max, 100, endpoint=True)

    nfreq = cube.shape[0]
    print nfreq, cube_min, cube_max
    freq_axis = cube.get_axis('freq')
    ra_axis = cube.get_axis('ra')
    dec_axis = cube.get_axis('dec')

    runlist = []
    for freqind in range(nfreq):
        cube_slice = cube[freqind, :, :]
        freq = freq_axis[freqind]
        fulltitle = tagname + " (freq = "+repr(freq/1e6)+" MHz)"
        runlist.append((freqind, cube_slice, ra_axis, dec_axis, coloraxis, fulltitle,
                        colorbar_title, fileprefix))

    count = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=count)
    pool.map(plot_single, runlist)

    command = ('ffmpeg','-y','-i','movies/'+tagname+'%03d.png',tagname+'.mp4') 
    subprocess.check_call(command)

root_selection = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
make_cube_movie(root_selection+"reg15separable.npy", 
                "15hr_Wigglez_separable_selection", 
                "# galaxies/pixel", "movies/15hr_Wigglez_separable_selection")
root_selection = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
make_cube_movie(root_selection+"reg15selection.npy", 
                "15hr_Wigglez_selection", 
                "# galaxies/pixel", "movies/15hr_Wigglez_selection")

root_selection = "/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/test/"
maplist = ["sec_A_15hr_41-73_cleaned_clean_map_I_with_B",
        "sec_A_15hr_41-73_cleaned_clean_map_I_with_C",
        "sec_A_15hr_41-73_cleaned_clean_map_I_with_D",
        "sec_B_15hr_41-73_cleaned_clean_map_I_with_A",
        "sec_B_15hr_41-73_cleaned_clean_map_I_with_C",
        "sec_B_15hr_41-73_cleaned_clean_map_I_with_D",
        "sec_C_15hr_41-73_cleaned_clean_map_I_with_A",
        "sec_C_15hr_41-73_cleaned_clean_map_I_with_B",
        "sec_C_15hr_41-73_cleaned_clean_map_I_with_D",
        "sec_D_15hr_41-73_cleaned_clean_map_I_with_A",
        "sec_D_15hr_41-73_cleaned_clean_map_I_with_B",
        "sec_D_15hr_41-73_cleaned_clean_map_I_with_C",]
covlist = ["sec_A_15hr_41-73_cleaned_noise_inv_I_with_B",
        "sec_A_15hr_41-73_cleaned_noise_inv_I_with_C",
        "sec_A_15hr_41-73_cleaned_noise_inv_I_with_D",
        "sec_B_15hr_41-73_cleaned_noise_inv_I_with_A",
        "sec_B_15hr_41-73_cleaned_noise_inv_I_with_C",
        "sec_B_15hr_41-73_cleaned_noise_inv_I_with_D",
        "sec_C_15hr_41-73_cleaned_noise_inv_I_with_A",
        "sec_C_15hr_41-73_cleaned_noise_inv_I_with_B",
        "sec_C_15hr_41-73_cleaned_noise_inv_I_with_D",
        "sec_D_15hr_41-73_cleaned_noise_inv_I_with_A",
        "sec_D_15hr_41-73_cleaned_noise_inv_I_with_B",
        "sec_D_15hr_41-73_cleaned_noise_inv_I_with_C"]

for tagname in maplist:
    make_cube_movie(root_selection+tagname+".npy",
                tagname, "Temperature", "movies/"+tagname)
for tagname in covlist:
    make_cube_movie(root_selection+tagname+".npy",
                tagname, "Covariance", "movies/"+tagname)

#command = ('mencoder mf://*.png -mf type=png:w=800:h=600:fps=25 -ovc lavc -lavcopts vcodec=mpeg4 -oac copy -of mpeg1video -o output.mpg')
#os.spawnvp(os.P_WAIT, 'mencoder', command)


