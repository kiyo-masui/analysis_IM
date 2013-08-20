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

live_sessions = range(41, 81+1)
live_sessions.remove(54)
live_sessions.remove(58)
live_sessions.remove(60)
live_sessions.remove(63)
live_sessions.remove(67)
live_sessions.remove(70)
live_sessions.remove(72)
live_sessions.remove(78)
live_sessions.remove(79)

maplist = ["sess_%d_calib_15hr_41-73_clean_map_I.npy" % item for item in live_sessions] 
input_root = '/mnt/raid-project/gmrt/kiyo/gbt_out/maps/apr16.2012/'
fullpath_maplist = [ input_root + item for item in maplist] 

weightlist = ["sess_%d_calib_15hr_41-73_weight_I.npy" % item for item in live_sessions]
output_root = "/mnt/raid-project/gmrt/eswitzer/GBT/maps/self_calibrated_2/"
weightlist_fullpath = [ output_root + item for item in weightlist]

reference_clean = [input_root + "sess_most_calib_15hr_41-73_clean_map_I.npy"]
reference_weight = [output_root + "sess_most_calib_15hr_41-73_weight_I.npy"]

#for filename in weightlist_fullpath:
for filename in reference_weight:
    outputdir = "/cita/d/www/home/eswitzer/movies/"
    tag = ".".join(filename.split(".")[:-1])  # extract root name
    tag = tag.split("/")[-1]
    tag += "_v2"
    pc.make_cube_movie(filename, "Weight", pc.cube_frame_dir,
                        sigmarange=-1, outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="", tag=tag)

#for filename in fullpath_maplist:
for filename in reference_clean:
    outputdir = "/cita/d/www/home/eswitzer/movies/"
    tag = ".".join(filename.split(".")[:-1])  # extract root name
    tag = tag.split("/")[-1]
    tag += "_v2"
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="", tag=tag)

