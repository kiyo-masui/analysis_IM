#Takes raw fits files.  Divides by scaled bandpass file, saved
# in .npy format as a (b, f) array, where b is the number of beams
# and f is the number of frequency channels.  
# The fits files should each have 13 beams ('BEAM') and 2 polarizations ('CRVAL4').

import core.algebra as al
import numpy as np
import h5py
from itertools import combinations
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

def rescale_fits_list(orig_dir, new_dir, x_bps, y_bps, x_bp_beams, y_bp_beams):
    

def rescale_data(data, x_bp_dat, y_bp_dat, x_beams, y_beams):
    x_beams = np.array(x_beams) - 1
    y_beams = np.array(y_beams) - 1
