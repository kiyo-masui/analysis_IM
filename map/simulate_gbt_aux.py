import scipy as sp
import numpy as np
from simulations import corr21cm
import core.algebra as algebra
import map.beam as beam
from core import constants as cc
from utils import units
import multiprocessing
from utils import data_paths
import sys
from numpy import random
import struct


def generate_delta_sim(input_file, output_file):
    r"""make the map with the temperature divided out (delta)"""
    print "reading %s -> %s (dividing by T_b(z))" % (input_file, output_file)

    simmap = algebra.make_vect(algebra.load(input_file))
    freq_axis = simmap.get_axis('freq') / 1.e6
    z_axis = units.nu21 / freq_axis - 1.0

    simobj = corr21cm.Corr21cm()
    T_b = simobj.T_b(z_axis)*1e-3

    simmap /= T_b[:, np.newaxis, np.newaxis]

    print "saving to" + output_file
    algebra.save(output_file, simmap)


def generate_proc_sim(input_file, weightfile, output_file,
                      meansub=False, degrade=False):
    r"""make the maps with various combinations of beam conv/meansub"""
    print "%s -> %s (beam, etc.)" % (input_file, output_file)
    simmap = algebra.make_vect(algebra.load(input_file))

    if meansub:
        print "performing mean subtraction"
        noise_inv = algebra.make_vect(algebra.load(weightfile))
        means = sp.sum(sp.sum(noise_inv * simmap, -1), -1)
        means /= sp.sum(sp.sum(noise_inv, -1), -1)
        means.shape += (1, 1)
        simmap -= means

    if degrade:
        print "performing common resolution convolution"
        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
        freq_data *= 1.0e6
        beam_diff = sp.sqrt(max(1.1 * beam_data) ** 2 - (beam_data) ** 2)
        common_resolution = beam.GaussianBeam(beam_diff, freq_data)
        # Convolve to a common resolution.
        simmap = common_resolution.apply(simmap)

    print "saving to" + output_file
    algebra.save(output_file, simmap)


def generate_aux_simset(fieldname, silent=False):

    datapath_db = data_paths.DataPath()
    weightfile = datapath_db.fetch("GBT_15hr_map_combined_cleaned_0mode_weight",
                                   intend_read=True, silent=silent)

    input_rawsimset = datapath_db.fetch("%s" % fieldname, intend_read=True,
                                        silent=silent)

    output_deltasimset = datapath_db.fetch("%s_delta" % fieldname,
                                          intend_write=True, silent=silent)

    input_beamsimset = datapath_db.fetch("%s_beam" % fieldname,
                                         intend_read=True, silent=silent)

    output_meansubsimset = datapath_db.fetch("%s_beam_meansub" % fieldname,
                                         intend_write=True, silent=silent)

    output_convsimset = datapath_db.fetch("%s_beam_conv" % fieldname,
                                         intend_write=True, silent=silent)

    output_meansubconvsimset = datapath_db.fetch("%s_beam_meansubconv" % fieldname,
                                         intend_write=True, silent=silent)

    for index in input_rawsimset[0]:
        generate_delta_sim(input_rawsimset[1][index],
                           output_deltasimset[1][index])

        generate_proc_sim(input_beamsimset[1][index], weightfile,
                          output_meansubsimset[1][index],
                          meansub=True, degrade=False)

        generate_proc_sim(input_beamsimset[1][index], weightfile,
                          output_convsimset[1][index],
                          meansub=False, degrade=True)

        generate_proc_sim(input_beamsimset[1][index], weightfile,
                          output_meansubconvsimset[1][index],
                          meansub=True, degrade=True)


def extend_full_simset(fieldlist):
    for fieldname in fieldlist:
        generate_aux_simset("sim_%s" % fieldname)
        generate_aux_simset("simideal_%s" % fieldname)
        generate_aux_simset("simvel_%s" % fieldname)
        #generate_simset(fieldname, scenario="ideal")
        #generate_simset(fieldname)
        #generate_simset(fieldname, scenario="streaming")


if __name__ == '__main__':
    #extend_full_simset(['15hr', '22hr', '1hr'])
    extend_full_simset(['15hr'])

