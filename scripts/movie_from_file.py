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
#filelist = ['secA_15hr_41-90_clean_map_I_762.npy',
#            'secB_15hr_41-90_clean_map_I_762.npy',
#            'secC_15hr_41-90_clean_map_I_762.npy',
#            'secD_15hr_41-90_clean_map_I_762.npy']
#dirname = '/mnt/raid-project/gmrt/tcv/maps/'
#dirname = "/mnt/raid-project/gmrt/kiyo/gbt_out_new/maps/may04.2012/"
#dirname = "/mnt/raid-project/gmrt/kiyo/gbt_out_new/maps/june25/"
#filelist = ["secA_15hr_41-90_clean_map_I_762_1024_225.npy",
#            "secA_15hr_41-90_clean_map_I_762_1024_225.npy",
#            "secA_15hr_41-90_clean_map_I_762_1024_225.npy",
#            "secA_15hr_41-90_clean_map_I_762_1024_225.npy"]
#filelist = ["15hr_41-90_fdgp_RM_clean_map_V.npy",
#            "15hr_41-90_fdgp_RM_clean_map_U.npy",
#            "15hr_41-90_fdgp_RM_clean_map_Q.npy",
#            "15hr_41-90_fdgp_RM_clean_map_I.npy",
#            "15hr_41-90_avg_fdg_clean_map_V.npy",
#            "15hr_41-90_avg_fdg_clean_map_U.npy",
#            "15hr_41-90_avg_fdg_clean_map_Q.npy",
#            "15hr_41-90_avg_fdg_clean_map_I.npy",
#            "15hr_41-90_fdgp_clean_map_V.npy",
#            "15hr_41-90_fdgp_clean_map_U.npy",
#            "15hr_41-90_fdgp_clean_map_Q.npy",
#            "15hr_41-90_fdgp_clean_map_I.npy"]

#dirname = "/mnt/raid-project/gmrt/kiyo/gbt_out_fg/maps/july9/"
#filelist = ["secA_15hr_41-90_clean_map_I_762.npy",
#            "secB_15hr_41-90_clean_map_I_762.npy",
#            "secC_15hr_41-90_clean_map_I_762.npy",
#            "secD_15hr_41-90_clean_map_I_762.npy"]
#full_list = [ dirname + filename for filename in filelist ]
#cbtitle = "Temperature (mK)"
#multiplier = 1000.
#sigmarange = 3.

full_list = ['/mnt/raid-project/gmrt/tcv/maps/15hr_41-80_avg_fdgp_new/secA_15hr_41-80_avg_fdgp_new_clean_map_I_800.npy',
             '/mnt/raid-project/gmrt/tcv/maps/15hr_41-80_avg_fdgp_new/secA_15hr_41-80_avg_fdgp_new_clean_map_Q_800.npy',
             '/mnt/raid-project/gmrt/tcv/maps/15hr_41-80_avg_fdgp_new/secA_15hr_41-80_avg_fdgp_new_clean_map_U_800.npy',
             '/mnt/raid-project/gmrt/tcv/maps/15hr_41-80_avg_fdgp_new/secA_15hr_41-80_avg_fdgp_new_clean_map_V_800.npy']

#full_list = ['/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/3903866168_simulation/gaussian_signal_simulation.npy',
#             '/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/3903866168_simulation/sec_A_15hr_41-90_clean_map_I.npy',
#             '/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_oldcal/3903866168_simulation/sec_A_15hr_41-90_clean_map_I_thermal.npy']

#full_list = ["/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_delta/15hr_v2/reg15data.npy",
#             "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned/15hr_v2/reg15separable.npy",
#             "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned/15hr_v2/reg15selection.npy"]
#full_list = ["/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_delta/1hr_v2/reg01data.npy",
#             "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned/1hr_v2/reg01separable.npy",
#             "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned/1hr_v2/reg01selection.npy"]
#full_list = ["/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secA_1hr_41-90_clean_map_I_800.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal/secA_1hr_41-18_clean_map_I_800.npy"]
#full_list = ["/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secA_1hr_41-90_clean_map_I_800.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secA_1hr_41-90_noise_weight_I_800.npy"]
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secA_1hr_41-90_noise_weight_I_800.npy"]
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secA_1hr_41-90_noise_weight_I_800.npy"]
#full_list = ["/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secB_1hr_41-90_noise_weight_I_800.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secC_1hr_41-90_noise_weight_I_800.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/1hr_oldcal_july16/secD_1hr_41-90_noise_weight_I_800.npy"]
#full_list = ["/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secA_15hr_41-90_clean_map_I_all.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secB_15hr_41-90_clean_map_I_all.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secC_15hr_41-90_clean_map_I_all.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secD_15hr_41-90_clean_map_I_all.npy",
#            "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secA_15hr_41-90_noise_inv_diag_I_all.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secB_15hr_41-90_noise_inv_diag_I_all.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secC_15hr_41-90_noise_inv_diag_I_all.npy",
#             "/mnt/raid-project/gmrt/eswitzer/GBT/maps/15hr_optimal_july11/secD_15hr_41-90_noise_inv_diag_I_all.npy"]

cbtitle = ""
multiplier = 1.
#sigmarange = -1.
sigmarange = 4.

for filename in full_list:
    outputdir = "/cita/d/www/home/eswitzer/movies/"
    tag = ".".join(filename.split(".")[:-1])  # extract root name
    tag = tag.split("/")[-1]
    pc.make_cube_movie(filename, cbtitle, pc.cube_frame_dir,
                        sigmarange=sigmarange, outputdir=outputdir,
                        multiplier=multiplier, logscale=False,
                        transverse=False, filetag_suffix="", tag=tag)

