#! /usr/bin/env python 

import os
from plotting import plot_cube


#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_4goodmode/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_25modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_4gwj'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_5goodmode_weightedjump/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_25modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gwj'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_20goodmode_weightedjump/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_15modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_20gwj'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_5goodmode_jump/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_25modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gj'

#maproot = "/mnt/raid-project/gmrt/ycli/map_result/maps/may10.2012/"
#mapname = "secD_1hr_41-90_clean_map_I_762.npy"
#tag = mapname.split('.')[0]

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_5gwj/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_90modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_norlegendre_5gwj'

#maproot = "/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_optimalmap_selfcal_762/"
#mapname = "sec_A_cleaned_clean_map_I_with_B_45modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_optimalmap_selfcal_762'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_15hr_map_oldcal_legendre_modes_5gwj/simmapmode_simmap/"
#mapname = "combined_clean_map_15modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gwj0_simmapmode_simmap'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_15hr_map_oldcal_legendre_modes_5gwj/mapmode_map/"
#mapname = "combined_clean_map_15modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gwj0_mapmode_map'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_15hr_map_oldcal_legendre_modes_5gwj/simmapmode_sim/"
#mapname = "combined_clean_map_15modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gwj0_simmapmode_sim'

#maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/GBT_1hr_map_oldcal_legendre_modes_0gwj/mapmode_map/"
#mapname = "combined_clean_map_50modes.npy"
#tag = mapname.split('.')[0]
#tag = tag + '_legendre_5gwj0_simmapmode_simmap'

maproot = "/mnt/data-pen3/ycli/map_result/maps/parkes/"
#mapname = "fir_dirty_map_I_1315"
mapname = "fir_parkes_2008_09_12_west_clean_map_I_1315.npy"
tag = "P080912W_clean_map_I"

#movieroot = "/mnt/raid-project/gmrt/ycli/movies/"
movieroot = "/cita/d/www/home/ycli/movies/"
cube_frame_dir = "/mnt/raid-project/gmrt/ycli/cube_frames/"
plot_cube.make_cube_movie(maproot + mapname, "Temperature (mK)", cube_frame_dir,
                          sigmarange=3., outputdir=movieroot, multiplier=1.,
                          transverse=False, tag=tag, convolve=True,
                          physical=False, beam_data=[0.233, 0.233],
                          freq_data=[1280,1350])
