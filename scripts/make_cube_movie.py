from optparse import OptionParser
import os
import numpy as np
from core import algebra
from plotting import plot_cube
from utils import data_paths
import string
#for saving slice
#saveslice=127,
#saveslice_file="/cita/d/www/home/mufma/Instrumental.eps"

#add this to title: + ",(%(freq)3.1f MHz,z=%(redshift)3.3f)"

# i also changed plot_cube.py RA went to X Dec went to Y (when plotting slices)


def plot_3D_map(filename, file_directory="/tmp/mufma/data/",
             title_name="Video", scale_name="Temperature(mK)",
             frame_directory="/tmp/mufma/data/", data_type="npy",
             output_directory="/cita/d/www/home/mufma/movies/",
             multiplier_num=1000.):
    r"""
    (data_set) -> jpeg-video

    Returns a 3D_map to output_directory using filename

    >>> plot_3D_map(data.npy)
    """

    # TODO: determine how to specify title
    datapath_db = data_paths.DataPath()
    saveslice_file = "/cita/d/www/home/mufma/Halos_40_bins.eps"
    saveslice = None
    frame_dir = frame_directory
    given_tag = string.rstrip(filename, data_type)
    plot_cube.make_cube_movie(file_directory+filename, scale_name,
                              frame_dir,
                              saveslice=saveslice,
                              saveslice_file=saveslice_file,
                              sigmarange=[-0.05, 0.05],
                              multiplier=multiplier_num,
                              title=title_name,
                              outputdir=output_directory,
                              tag=given_tag)


def main():
    r"""can call the plotting function directly from the shell"""
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("--od", "--output_directory",
                      action="store",
                      type="string",
                      default="/cita/d/www/home/mufma/movies/",
                      dest="output_directory",
                      help="give an output_directory")
    parser.add_option("-f", "--frame_dir",
                      action="store",
                      type="string",
                      default="/tmp/mufma/",
                      dest="frame_directory",
                      help="give a frame_directory")
    parser.add_option("-s", "--scale",
                      action="store",
                      type="string",
                      default="Temperature(mK)",
                      dest="scale_name",
                      help="give a scale(Temperature or Volume or etc)")
    parser.add_option("-t", "--title",
                      action="store",
                      type="string",
                      default="Movie_slices_3Dmap",
                      dest="title_name",
                      help="give a title")
    parser.add_option("-d", "--destination",
                      action="store",
                      type="string",
                      default="/tmp/mufma/data/",
                      dest="file_directory",
                      help="give a file_directory")
    parser.add_option("--dt", "--data_type",
                      action="store",
                      type="string",
                      default=".npy",
                      dest="data_type",
                      help="give a data_type")
    parser.add_option("--multiplier",
                      action="store",
                      type="float",
                      default=1.,
                      dest="multiplier_num",
                      help="multiply the units of scale")

    (options, args) = parser.parse_args()
    optdict = vars(options)
    print args
    print optdict
    if len(args) != 1:
        parser.error("wrong number of arguments")
    plot_3D_map(*args, **optdict)


if __name__ == '__main__':
    main()
