import scipy as sp
import numpy as np
import core.algebra as algebra
import map.beam as beam
from core import constants as cc
import multiprocessing
from utils import data_paths
import sys
from numpy import random
import struct
from kiyopy import parse_ini
import kiyopy.utils
from utils import units
from plotting import plot_cube as pc


#filelist = ['15hr_41-90_flux_clean_map_V.npy',
#            '15hr_41-90_flux_clean_map_U.npy',
#            '15hr_41-90_flux_clean_map_Q.npy',
#            '15hr_41-90_flux_clean_map_I.npy',
#            '15hr_41-90_tc_clean_map_V.npy',
#            '15hr_41-90_tc_clean_map_U.npy',
#            '15hr_41-90_tc_clean_map_Q.npy',
#            '15hr_41-90_tc_clean_map_I.npy',
#            '15hr_41-90_fdg_clean_map_V.npy',
#            '15hr_41-90_fdg_clean_map_U.npy',
#            '15hr_41-90_fdg_clean_map_Q.npy',
#            '15hr_41-90_fdg_clean_map_I.npy']
#filelist = ['15hr_41-90_fdg_clean_map_V.npy',
#            '15hr_41-90_fdg_clean_map_U.npy',
#            '15hr_41-90_fdg_clean_map_Q.npy',
#            '15hr_41-90_fdg_clean_map_I.npy']
#filelist = ['15hr_41-90_RM_clean_map_I.npy',
#            '15hr_41-90_RM_clean_map_Q.npy',
#            '15hr_41-90_RM_clean_map_U.npy',
#            '15hr_41-90_RM_clean_map_V.npy']
filelist = ['secA_15hr_41-90_clean_map_I_762.npy',
            'secB_15hr_41-90_clean_map_I_762.npy',
            'secC_15hr_41-90_clean_map_I_762.npy',
            'secD_15hr_41-90_clean_map_I_762.npy']
#dirname = '/mnt/raid-project/gmrt/tcv/maps/'
dirname = "/mnt/raid-project/gmrt/kiyo/gbt_out_new/maps/may04.2012/"

for filenam in filelist:
    outputdir = "/cita/d/www/home/eswitzer/movies/"
    tag = ".".join(filenam.split(".")[:-1])  # extract root name
    tag = tag.split("/")[-1]
    tag += "_v2"
    pc.make_cube_movie(dirname + filenam, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="", tag=tag)

