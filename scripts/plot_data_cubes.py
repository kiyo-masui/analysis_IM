r"""Several scripts to call plot_cube_movie and make movies of GBT data,
simulations, etc.
"""
import numpy as np
from utils import data_paths
from plotting import plot_cube as pc
from core import algebra
import sys

cube_frame_dir = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/cube_frames/"
output_dir = ""

def plot_GBT_maps(keyname, transverse=False, skip_noise=False):
    datapath_db = data_paths.DataPath()

    section_list = ['A', 'B', 'C', 'D']
    for section in section_list:
        filename = datapath_db.fetch(keyname,
                                     pick=(section + '_clean_map'))
        title = "Sec. %s, %s" % (section, keyname)
        pc.make_cube_movie(filename,
                           "Temperature (mK)", cube_frame_dir, sigmarange=5.,
                           outputdir="./", multiplier=1000.,
                           transverse=transverse,
                           title=title)

        filename = datapath_db.fetch(keyname,
                                     pick=(section + '_dirty_map'))
        title = "Sec. %s, %s (dirty)" % (section, keyname)
        pc.make_cube_movie(filename,
                           "Temperature (mK)", cube_frame_dir, sigmarange=5.,
                           outputdir="./", multiplier=1000.,
                           transverse=transverse,
                           title=title)

        if not skip_noise:
            filename = datapath_db.fetch(keyname,
                                         pick=(section + '_noise_diag'))
            title = "Sec. %s, %s (noise)" % (section, keyname)
            pc.make_cube_movie(filename,
                               "Covariance", cube_frame_dir, sigmarange=-1,
                               outputdir="./", multiplier=1000.,
                               transverse=transverse,
                               title=title)


def plot_GBT_mapset():
    plot_GBT_maps('GBT_15hr_map_proposal', transverse=True, skip_noise=True)
    plot_GBT_maps('GBT_15hr_map_proposal', transverse=False, skip_noise=True)
    plot_GBT_maps('GBT_15hr_map', transverse=True)
    plot_GBT_maps('GBT_15hr_map', transverse=False)
    plot_GBT_maps('GBT_22hr_map', transverse=True)
    plot_GBT_maps('GBT_22hr_map', transverse=False)
    plot_GBT_maps('GBT_1hr_map', transverse=True)
    plot_GBT_maps('GBT_1hr_map', transverse=False)


def plot_simulations(keyname, transverse=False):
    """make movies of the 15hr simulations
    permutations: with or without streaming, including beam, adding real data
    """
    datapath_db = data_paths.DataPath()

    filename = datapath_db.fetch(keyname, pick='15')
    pc.make_cube_movie(filename,
                       "Temperature (mK)", cube_frame_dir, sigmarange=5.,
                       outputdir="./", multiplier=1000., transverse=transverse)

    #pc.make_cube_movie(root_directory + "sim_beam_" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=5.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #pc.make_cube_movie(root_directory + "sim_beam_plus_data" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=10.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #pc.make_cube_movie(root_directory + "simvel_" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=5.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #pc.make_cube_movie(root_directory + "simvel_beam_" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=5.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)
    #pc.make_cube_movie(root_directory + "simvel_beam_plus_data" + suffix,
    #                   "Temperature (mK)", cube_frame_dir, sigmarange=10.,
    #                   outputdir="./", multiplier=1000., transverse=transverse)


def plot_GBT_simset():
    plot_simulations('sim_15hr', transverse=True)
    plot_simulations('sim_15hr', transverse=False)


def plot_difference(filename1, filename2, title, sigmarange=6., sigmacut=None,
                    transverse=False, outputdir="./", multiplier=1000.,
                    logscale=False, fractional=False, filetag_suffix="",
                    ignore=None, diff_filename="./difference.npy"):
    """make movies of the difference of two maps (assuming same dimensions)"""
    map1 = algebra.make_vect(algebra.load(filename1))
    map2 = algebra.make_vect(algebra.load(filename2))

    if fractional:
        difftitle = "fractional diff."
        dmap = (map1-map2)/map1*100.
    else:
        difftitle = "difference"
        dmap = map1-map2

    algebra.save(diff_filename, dmap)

    pc.make_cube_movie(diff_filename,
                       difftitle, cube_frame_dir, sigmarange=6,
                       sigmacut=sigmacut, outputdir=outputdir, ignore=ignore,
                       multiplier=multiplier, transverse=transverse,
                       logscale=False)

    pc.make_cube_movie(filename1,
                       title, cube_frame_dir, sigmarange=sigmarange,
                       sigmacut=sigmacut, outputdir=outputdir, ignore=ignore,
                       multiplier=multiplier, transverse=transverse,
                       logscale=logscale, filetag_suffix="_1")

    pc.make_cube_movie(filename2,
                       title, cube_frame_dir, sigmarange=sigmarange,
                       sigmacut=sigmacut, outputdir=outputdir, ignore=ignore,
                       multiplier=multiplier, transverse=transverse,
                       logscale=logscale, filetag_suffix="_2")


def plot_GBT_diff_tests():
    tcv_15root = "/mnt/raid-project/gmrt/tcv/"
    tcv_15root += "modetest/73_ABCD_all_15_modes_real_maponly/"
    tcv_15map = tcv_15root + "sec_A_15hr_41-90_cleaned_clean_map_I_with_B.npy"
    tcv_15noise = tcv_15root + "sec_A_15hr_41-90_cleaned_noise_inv_I_with_B.npy"
    ers_15root = "/mnt/raid-project/gmrt/eswitzer/GBT/"
    ers_15root += "cleaned_maps/freq_slices_refactor_tests_15modes/"
    ers_15map = ers_15root + "sec_A_15hr_41-90_cleaned_clean_map_I_with_B.npy"
    ers_15noise = ers_15root + "sec_A_15hr_41-90_cleaned_noise_inv_I_with_B.npy"

    plot_difference(tcv_15map, ers_15map, "Temperature (mK)", sigmarange=6.,
                    fractional=False, diff_filename="./map_difference.npy",
                    transverse=False)

    plot_difference(tcv_15map, ers_15map, "Temperature (mK)", sigmarange=6.,
                    fractional=False, diff_filename="./map_difference.npy",
                    transverse=True)

    plot_difference(tcv_15noise, ers_15noise, "log inv. covariance", sigmarange=-1.,
                    multiplier=1., logscale=True, fractional=True,
                    diff_filename="./noise_inv_fractional_difference.npy",
                    transverse=False)

    plot_difference(tcv_15noise, ers_15noise, "log inv. covariance", sigmarange=-1.,
                    multiplier=1., logscale=True, fractional=True,
                    diff_filename="./noise_inv_fractional_difference.npy",
                    transverse=True)


if __name__ == "__main__":
    #plot_GBT_mapset():
    #plot_GBT_simset():
    plot_GBT_diff_tests()
