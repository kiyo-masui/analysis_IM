import numpy as np
from core import algebra
from utils import data_paths
from plotting import plot_cube

def plot_first_few_amps(dbitem, modelist):
    datapath_db = data_paths.DataPath()
    frame_dir = "/scratch/eswitzer/cube_frames/"

    dbname = "db:%s:A_with_B;modes;100modes" % dbitem

    for modeindex in modelist:
        filename = "%s_%d.eps" % (dbitem, modeindex)
        plot_cube.make_cube_movie(dbname,
                        "Amplitude", frame_dir,
                        saveslice=modeindex, saveslice_file=filename,
                        sigmarange=-1, multiplier=1.,
                        title="%s" % filename)


if __name__ == "__main__":
    plot_first_few_amps("GBT_15hr_map_oldcal_cleaned", range(5))
    #plot_first_few_amps("GBT_15hr_map_mapcal2_cleaned", range(5))
