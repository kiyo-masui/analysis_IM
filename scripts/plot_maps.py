import numpy as np
from core import algebra
from utils import data_paths
from plotting import plot_cube

if __name__ == "__main__":
    datapath_db = data_paths.DataPath()
    frame_dir = "/scratch/eswitzer/cube_frames/"

    plot_cube.make_cube_movie("db:GBT_15hr_map_oldcal:A;clean_map",
                    "Temperature (mK)", frame_dir,
                    saveslice=127, saveslice_file="GBT_15hr_map_secA.eps",
                    sigmarange=3., multiplier=1000.,
                    title="GBT 15hr field (%(freq)3.1f MHz, z = %(redshift)3.3f)")

    plot_cube.make_cube_movie("db:GBT_15hr_map_oldcal_cleaned_noconv_combined:map;25modes",
                    "Temperature (mK)", frame_dir,
                    saveslice=127, saveslice_file="cleaned_15hr_map_25modes.eps",
                    sigmarange=0.4, multiplier=1000.,
                    title="GBT 15hr field, cleaned (%(freq)3.1f MHz, z = %(redshift)3.3f)")

    plot_cube.make_cube_movie("db:GBT_15hr_map_oldcal_cleaned_noconv_combined:map;25modes",
                    "Temperature (mK)", frame_dir,
                    saveslice=127, saveslice_file="cleaned_15hr_map_25modes_convolved.eps",
                    convolve=True,
                    sigmarange=0.4, multiplier=1000.,
                    title="GBT 15hr field, cleaned, beam convolved (%(freq)3.1f MHz, z = %(redshift)3.3f)")

