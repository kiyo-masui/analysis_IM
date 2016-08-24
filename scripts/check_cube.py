#! /usr/bin/env python 

import os
from plotting import plot_cube


maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_4goodmode/"
mapname = "sec_A_cleaned_clean_map_I_with_B_25modes.npy"

maproot = "/mnt/raid-project/gmrt/ycli/foreground_cleand/mode_changed_5goodmode_weightedjump/"
mapname = "sec_A_cleaned_clean_map_I_with_B_10modes.npy"

movieroot = "/mnt/raid-project/gmrt/ycli/movies/"
tag = mapname.split('.')[0]
tag = tag + '_legendre_5gwj'

cube_frame_dir = "/mnt/raid-project/gmrt/ycli/cube_frames/"
plot_cube.make_cube_movie(maproot + mapname, "Temperature (mK)", cube_frame_dir,
                          sigmarange=3., outputdir=movieroot, multiplier=1000.,
                          transverse=False, tag=tag,
                          physical=False)
