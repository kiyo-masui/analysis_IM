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

cube_root = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/cube_frames/"

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


def make_cube_movie(filename, tag, colorbar_title, fileprefix,
                    sigmarange=3., ignore=None, multiplier=1.):
    """Make a stack of spatial slice maps and animate them"""
    cube = algebra.make_vect(algebra.load(filename)) * multiplier
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

    subprocess.check_call(('ffmpeg', '-r', '10', '-y', '-i', cube_root + tag +\
               '%03d.png', tag + '.mp4'))

# make plots of the 22hr selection
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
make_cube_movie(root_directory + "reg22selection.npy",
                "reg22selection",
                "Selection", cube_root + "reg22selection.npy",
                sigmarange=None)
make_cube_movie(root_directory + "reg22separable.npy",
                "reg22separable",
                "Separable Selection", cube_root + "reg22separable.npy",
                sigmarange=None)
sys.exit()

# make plots of the 15hr field
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
make_cube_movie(root_directory + "combined_41-73_cleaned_clean_test.npy",
                "combined_41-73_cleaned_clean_test",
                "Temperature (mK)", cube_root + "combined_41-73_cleaned_clean_test",
                sigmarange=6., multiplier = 1000.)
make_cube_movie(root_directory + "combined_41-73_cleaned_noise_inv_test.npy",
                "combined_41-73_cleaned_noise_inv_test",
                "Covariance", cube_root + "combined_41-73_cleaned_noise_inv_test",
                sigmarange=-1)

# make plots of the 22hr field
#root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/corr/84_ABCD_all_15_modes/"
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
make_cube_movie(root_directory + "combined_22hr_41-84_cleaned_clean.npy",
                "combined_22hr_41-84_cleaned_clean.npy",
                "Temperature (mK)", cube_root + "combined_22hr_41-84_cleaned_clean.npy",
                sigmarange=6, multiplier = 1000.)

make_cube_movie(root_directory + "combined_22hr_41-84_cleaned_noise_inv.npy",
                "combined_22hr_41-84_cleaned_noise_inv.npy",
                "Covariance", cube_root + "combined_22hr_41-84_cleaned_noise_inv.npy",
                sigmarange=-1)

sys.exit()
root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/test100/"
make_cube_movie(root_directory + "simulated_signal_map_1.npy",
                "simulated_signal_map_1",
                "Temperature", cube_root + "simulated_signal_map_1",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/calinliv/wiggleZ/simulations/test100/"
make_cube_movie(root_directory + "simulated_signal_map_1_with_beam.npy",
                "simulated_signal_map_1_with_beam",
                "Temperature", cube_root + "simulated_signal_map_1_with_beam",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/modetest/73_ABCD_all_0_modes_sim_maponly_NOCONV/"
make_cube_movie(root_directory + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy",
                "sec_A_15hr_41-73_cleaned_clean_map_I_with_B_noconv",
                "Temperature", cube_root + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B_noconv",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/modetest/73_ABCD_all_0_modes_sim_maponly/"
make_cube_movie(root_directory + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B.npy",
                "sec_A_15hr_41-73_cleaned_clean_map_I_with_B",
                "Temperature", cube_root + "sec_A_15hr_41-73_cleaned_clean_map_I_with_B",
                sigmarange=[-1.,1.])
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/modetest_combined_maps_0_50/"
make_cube_movie(root_directory + "combined_sim_41-73_cleaned_clean_0.npy",
                "combined_sim_41-73_cleaned_clean_0",
                "Temperature", cube_root + "combined_sim_41-73_cleaned_clean_0",
                sigmarange=[-1.,1.])
sys.exit()

make_cube_movie("delta_selection.npy",
                "delta_selection",
                "Difference in selection", cube_root + "delta_selection")
sys.exit()
make_cube_movie("reg15selection_est.npy",
                "reg15selection_est",
                "Full selection", cube_root + "reg15selection_est",
                sigmarange=None)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/noise_model/"
make_cube_movie(root_directory + "completeness_model_41-73_15modes.npy",
                "completeness_model_41-73_15modes",
                "Completeness", cube_root + "completeness_model_41-73_15modes",
                sigmarange=5., ignore=-1.)
sys.exit()
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/noise_model/"
make_cube_movie(root_directory + "noise_model_41-73_15modes.npy",
                "noise_model_41-73_15modes",
                "Temperature", cube_root + "noise_model_41-73_15modes",
                sigmarange=5., ignore=-1.)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
make_cube_movie(root_directory + "combined_41-73_clean_test.npy",
                "combined_41-73_clean_test",
                "Temperature", cube_root + "combined_41-73_clean_test",
                sigmarange=8.)
make_cube_movie(root_directory + "combined_41-73_noise_inv_test.npy",
                "combined_41-73_noise_inv_test",
                "Covariance", cube_root + "combined_41-73_noise_inv_test",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
make_cube_movie(root_directory + "combined_41-73_cleaned_clean_test.npy",
                "combined_41-73_cleaned_clean_test",
                "Temperature", cube_root + "combined_41-73_cleaned_clean_test",
                sigmarange=4.)
make_cube_movie(root_directory + "combined_41-73_cleaned_noise_inv_test.npy",
                "combined_41-73_cleaned_noise_inv_test",
                "Covariance", cube_root + "combined_41-73_cleaned_noise_inv_test",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/"
make_cube_movie(root_directory + \
                "sec_B_15hr_41-69_cleaned_clean_map_I_with_A.npy",
                "sec_B_15hr_41-69_cleaned_clean_map_I_with_A",
                "Temperature",
                cube_root + "sec_B_15hr_41-69_cleaned_clean_map_I_with_A",
                sigmarange=[-0.001, 0.001])
make_cube_movie(root_directory + \
                "sec_A_15hr_41-69_cleaned_clean_map_I_with_B.npy",
                "sec_A_15hr_41-69_cleaned_clean_map_I_with_B",
                "Temperature",
                cube_root + "sec_A_15hr_41-69_cleaned_clean_map_I_with_B",
                sigmarange=[-0.001, 0.001])
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/combined_maps/"
make_cube_movie(root_directory + "combined_41-69_cleaned_product.npy",
                "combined_41-69_cleaned_product",
                "Temperature", cube_root + "combined_41-69_cleaned_product",
                sigmarange=20.)
make_cube_movie(root_directory + "combined_41-69_cleaned_clean.npy",
                "combined_41-69_cleaned_clean",
                "Temperature", cube_root + "combined_41-69_cleaned_clean",
                sigmarange=[-0.01, 0.01])
make_cube_movie(root_directory + "combined_41-69_cleaned_noise_inv.npy",
                "combined_41-69_cleaned_noise_inv",
                "Covariance", cube_root + "combined_41-69_cleaned_noise_inv",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/kiyo/wiggleZ/corr/"
make_cube_movie(root_directory + "sec_A_15hr_41-69_cleaned_clean_map_I.npy",
                "sec_A_15hr_41-69_cleaned_clean_map_I",
                "Temperature", cube_root + "sec_A_15hr_41-69_cleaned_clean_map_I",
                sigmarange=[-0.01, 0.01])
make_cube_movie(root_directory + "sec_A_15hr_41-69_cleaned_noise_inv_I.npy",
                "sec_A_15hr_41-69_cleaned_noise_inv_I",
                "Covariance", cube_root + "sec_A_15hr_41-69_cleaned_noise_inv_I",
                sigmarange=-1)
make_cube_movie(root_directory + "sec_B_15hr_41-69_cleaned_clean_map_I.npy",
                "sec_B_15hr_41-69_cleaned_clean_map_I",
                "Temperature", cube_root + "sec_B_15hr_41-69_cleaned_clean_map_I",
                sigmarange=[-0.01, 0.01])
make_cube_movie(root_directory + "sec_B_15hr_41-69_cleaned_noise_inv_I.npy",
                "sec_B_15hr_41-69_cleaned_noise_inv_I",
                "Covariance", cube_root + "sec_B_15hr_41-69_cleaned_noise_inv_I",
                sigmarange=-1)
sys.exit()

root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
make_cube_movie(root_directory + "reg15separable.npy",
                "15hr_Wigglez_separable_selection",
                "# galaxies/pixel", cube_root + "15hr_Wigglez_separable_selection")
root_directory = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
make_cube_movie(root_directory + "reg15selection.npy",
                "15hr_Wigglez_selection",
                "# galaxies/pixel", cube_root + "15hr_Wigglez_selection")

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
                tagname, "Temperature", cube_root + tagname)
for tagname in covlist:
    make_cube_movie(root_directory + tagname + ".npy",
                tagname, "Covariance", cube_root + tagname)
