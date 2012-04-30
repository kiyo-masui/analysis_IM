"""Make simulated signal realizations in the GBT survey volume"""
import scipy as sp
from simulations import corr21cm
#from simulations import foregroundsck
#from simulations import pointsource
import core.algebra as algebra
import map.beam as beam
from core import constants as cc
import multiprocessing
from utils import data_paths
import sys
from numpy import random
import struct


def test_scheme(template_file, sim_filename1, sim_filename2):
    r"""look at some differences between maps"""
    template_map = algebra.make_vect(algebra.load(template_file))
    gbtsim1 = realize_simulation(template_map, scenario='streaming',
                                seed=5489, refinement=1.)
    gbtsim2 = realize_simulation(template_map,
                                    seed=5489, refinement=1.)

    sim_map1 = algebra.make_vect(gbtsim1, axis_names=('freq', 'ra', 'dec'))
    sim_map2 = algebra.make_vect(gbtsim2, axis_names=('freq', 'ra', 'dec'))
    sim_map1.copy_axis_info(template_map)
    sim_map2.copy_axis_info(template_map)
    algebra.save(sim_filename1, sim_map1)
    algebra.save(sim_filename2, sim_map2)


def run_scheme_test():
    template_file = "/mnt/raid-project/gmrt/tcv/maps/sec_A_15hr_41-90_clean_map_I.npy"
    sim_filename1 = "sim_streaming1.npy"
    sim_filename2 = "sim_streaming2.npy"
    test_scheme(template_file, sim_filename1, sim_filename2)


if __name__ == '__main__':
    run_scheme_test()
