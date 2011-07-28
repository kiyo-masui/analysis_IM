import numpy as np
import sys
from core import algebra
from correlate import wigglez_xcorr as wxc
import copy
import shelve
# ugly odds and ends that are needed to read the pkl file
import cPickle
from correlate.freq_slices import *
import multiprocessing

def process_optical_to_delta(optical_file, optical_selection_file, outfile):
    print "-"*80
    print "in: " + optical_file
    print "nbar: " + optical_selection_file
    print "out: " + outfile
    map_opt = algebra.make_vect(algebra.load(optical_file))
    map_nbar = algebra.make_vect(algebra.load(optical_selection_file))

    # convert to delta-overdensity
    map_opt = map_opt / map_nbar - 1.
    #algebra.compressed_array_summary(map_opt, "opt after conversion to delta")

    # set the NaNs and infs to zero in data and weights
    nan_array = np.isnan(map_opt)
    map_opt[nan_array] = 0.
    map_nbar[nan_array] = 0.
    inf_array = np.isinf(map_opt)
    map_opt[inf_array] = 0.
    map_nbar[inf_array] = 0.

    algebra.save(outfile, map_opt)

optical_root = "/mnt/raid-project/gmrt/eswitzer/wiggleZ/binned_wiggleZ/"
optical_signal = "reg15data.npy"
selection_function = "reg15separable.npy"

process_optical_to_delta(optical_root + optical_signal,
                         optical_root + selection_function,
                         optical_root + "delta15data.npy")
random_list = [optical_root + "reg15rand%03d.npy" % index for index in range(0,100)]
output_list = [optical_root + "delta15rand%03d.npy" % index for index in range(0,100)]

[process_optical_to_delta(infile, optical_root + selection_function, outfile) for (infile, outfile) in zip(random_list, output_list)]

