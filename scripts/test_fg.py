"""Make simulated signal realizations in the GBT survey volume"""
import scipy as sp
import numpy as np
from simulations import corr21cm
from simulations import foregroundsck
from simulations import pointsource
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

def save_and_plot(array, template, filename):
    array = algebra.make_vect(array, axis_names=('freq', 'ra', 'dec'))
    array.copy_axis_info(template)

    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
    freq_data *= 1.0e6

    beamobj = beam.GaussianBeam(beam_data, freq_data)
    array_beam = beamobj.apply(array)

    algebra.save(filename, array_beam)
    outputdir = "/cita/d/www/home/eswitzer/movies/"
    pc.make_cube_movie(filename, "Temperature (mK)", pc.cube_frame_dir,
                        sigmarange=3., outputdir=outputdir, multiplier=1000.,
                        transverse=False, filetag_suffix="_trial")


if __name__ == '__main__':
    template_file = '/mnt/raid-project/gmrt/eswitzer/GBT/simulations/15hr_oldmap_ideal/sim_000.npy'
    template_map = algebra.make_vect(algebra.load(template_file))

    #simobj = corr21cm.Corr21cm.like_kiyo_map(template_map)
    #(gbtsim, gbtphys, physdim) = simobj.get_kiyo_field_physical(refinement=2)
    #save_and_plot(gbtsim, template_map, "skysim.npy")

    syncobj = foregroundsck.Synchrotron.like_kiyo_map(template_map)
    sync_field = syncobj.getfield() * 0.001
    save_and_plot(sync_field, template_map, "synsim.npy")

    ptsrcobj = pointsource.DiMatteo.like_kiyo_map(template_map)
    ptsrc_field = ptsrcobj.getfield()
    save_and_plot(ptsrc_field, template_map, "pssim.npy")
