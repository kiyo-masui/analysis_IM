"""tools for plotting data cubes; this code is new and in developement"""
import subprocess
import sys
import numpy as np
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
    cplot = plt.contourf(xaxis, yaxis, np.transpose(cube_slice), zaxis)
    plt.axis('scaled')
    plt.xlim((np.min(xaxis), np.max(xaxis)))
    plt.ylim((np.min(yaxis), np.max(yaxis)))
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title(title, fontsize=10)
    ##cbar = plt.colorbar(cplot, ticks=zaxis)
    cbar = plt.colorbar(cplot)
    cbar.ax.set_ylabel(cbar_title)
    filename = fileprefix + str('%03d' % index) + '.png'
    plt.savefig(filename, dpi=100)
    plt.clf()


def make_cube_movie(filename, tag, colorbar_title, fileprefix,
                    sigmarange=3.):
    """Make a stack of spatial slice maps and animate them"""
    cube = algebra.make_vect(algebra.load(filename))
    cube_mean = np.mean(cube)
    cube_std = np.std(cube)
    try:
        len(sigmarange)
        coloraxis = np.linspace(sigmarange[0], sigmarange[1],
                                500, endpoint=True)
    except TypeError:
        if (sigmarange > 0.):
            coloraxis = np.linspace(cube_mean - sigmarange * cube_std,
                                    cube_mean + sigmarange * cube_std,
                                    500, endpoint=True)
        else:
            coloraxis = np.linspace(np.min(cube),  np.max(cube),
                                    500, endpoint=True)

    freq_axis = cube.get_axis('freq')
    space_axis = (cube.get_axis('ra'), cube.get_axis('dec'))

    runlist = []
    for freqind in range(cube.shape[0]):
        fulltitle = tag + " (freq = " + \
                    repr(freq_axis[freqind] / 1.e6) + " MHz)"
        runlist.append((freqind, cube[freqind, :, :], space_axis[0],
                        space_axis[1], coloraxis, fulltitle,
                        colorbar_title, fileprefix))

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(plot_single, runlist)

    subprocess.check_call(('ffmpeg', '-y', '-i', 'movies/' + tag +\
               '%03d.png', tag + '.mp4'))


root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
make_cube_movie(root_directory + "combined_41-73_cleaned_clean_test.npy",
                "combined_41-73_cleaned_clean_test",
                "Temperature", "movies/combined_41-73_cleaned_clean_test",
                sigmarange=4.)
make_cube_movie(root_directory + "combined_41-73_cleaned_noise_inv_test.npy",
                "combined_41-73_cleaned_noise_inv_test",
                "Covariance", "movies/combined_41-73_cleaned_noise_inv_test",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/"
make_cube_movie(root_directory + \
                "sec_B_15hr_41-69_cleaned_clean_map_I_with_A.npy",
                "sec_B_15hr_41-69_cleaned_clean_map_I_with_A",
                "Temperature",
                "movies/sec_B_15hr_41-69_cleaned_clean_map_I_with_A",
                sigmarange=[-0.001, 0.001])
make_cube_movie(root_directory + \
                "sec_A_15hr_41-69_cleaned_clean_map_I_with_B.npy",
                "sec_A_15hr_41-69_cleaned_clean_map_I_with_B",
                "Temperature",
                "movies/sec_A_15hr_41-69_cleaned_clean_map_I_with_B",
                sigmarange=[-0.001, 0.001])
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
make_cube_movie(root_directory + "combined_41-69_cleaned_product.npy",
                "combined_41-69_cleaned_product",
                "Temperature", "movies/combined_41-69_cleaned_product",
                sigmarange=20.)
make_cube_movie(root_directory + "combined_41-69_cleaned_clean.npy",
                "combined_41-69_cleaned_clean",
                "Temperature", "movies/combined_41-69_cleaned_clean",
                sigmarange=[-0.01, 0.01])
make_cube_movie(root_directory + "combined_41-69_cleaned_noise_inv.npy",
                "combined_41-69_cleaned_noise_inv",
                "Covariance", "movies/combined_41-69_cleaned_noise_inv",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/"
make_cube_movie(root_directory + "sec_A_15hr_41-69_cleaned_clean_map_I.npy",
                "sec_A_15hr_41-69_cleaned_clean_map_I",
                "Temperature", "movies/sec_A_15hr_41-69_cleaned_clean_map_I",
                sigmarange=[-0.01, 0.01])
make_cube_movie(root_directory + "sec_A_15hr_41-69_cleaned_noise_inv_I.npy",
                "sec_A_15hr_41-69_cleaned_noise_inv_I",
                "Covariance", "movies/sec_A_15hr_41-69_cleaned_noise_inv_I",
                sigmarange=-1)
make_cube_movie(root_directory + "sec_B_15hr_41-69_cleaned_clean_map_I.npy",
                "sec_B_15hr_41-69_cleaned_clean_map_I",
                "Temperature", "movies/sec_B_15hr_41-69_cleaned_clean_map_I",
                sigmarange=[-0.01, 0.01])
make_cube_movie(root_directory + "sec_B_15hr_41-69_cleaned_noise_inv_I.npy",
                "sec_B_15hr_41-69_cleaned_noise_inv_I",
                "Covariance", "movies/sec_B_15hr_41-69_cleaned_noise_inv_I",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
make_cube_movie(root_directory + "reg15separable.npy",
                "15hr_Wigglez_separable_selection",
                "# galaxies/pixel", "movies/15hr_Wigglez_separable_selection")
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
make_cube_movie(root_directory + "reg15selection.npy",
                "15hr_Wigglez_selection",
                "# galaxies/pixel", "movies/15hr_Wigglez_selection")

root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/test/"
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
        "sec_D_15hr_41-73_cleaned_clean_map_I_with_C"]
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
    make_cube_movie(root_directory + tagname + ".npy",
                tagname, "Temperature", "movies/" + tagname)
for tagname in covlist:
    make_cube_movie(root_directory + tagname + ".npy",
                tagname, "Covariance", "movies/" + tagname)
