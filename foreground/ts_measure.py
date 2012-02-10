"""Module to measure the foreground spectral components in the time stream."""

import shelve
import multiprocessing as mp
import time as time_module

import numpy as np
import scipy.linalg as linalg
import numpy.ma as ma

from kiyopy import parse_ini, utils
import kiyopy.pickle_method
import kiyopy.utils
import kiyopy.custom_exceptions as ce
import core.fitsGBT

# XXX
import matplotlib.pyplot as plt


params_init = {
               # IO.
               "input_root" : "./testdata/",
               "file_middles" : ("testfile_guppi_combined",),
               "input_end" : ".fits",
               "output_root" : "./",
               "output_filename" : "foreground_modes.shelve",
               "scans" : (),
               "IFs" : ()
               }

prefix = 'tf_'

class Measure(object) :
    """Measures the noise of data files.
    """

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                      prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1) :
                
        params = self.params
        kiyopy.utils.mkparents(params['output_root'] + 
                               params['output_filename'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        file_middles = params['file_middles']
        covar_dict = {}
        counts_dict = {}

        for middle in file_middles:
            key, covar, counts = self.process_file(middle)
            if covar_dict.haskey(key):
                if covar_dict[key].shape != covar.shape:
                    msg = ("All data needs to have the same band and"
                           "polarization structure.")
                    raise ce.DataError(msg)
                covar_dict[key] += covar
                counts_dict[key] += counts
            else:
                covar_dict[key] = covar
                counts_dict[key] = counts

    def self.process_file(self, middle):
        params = self.params
        file_name = (params['input_root'] + file_middle
                     + params['input_end'])
        band_inds = params["IFs"]
        Reader = core.fitsGBT.Reader(file_name, feedback=self.feedback)
        n_bands = len(Reader.IF_set)
        if not band_inds:
            band_inds = range(n_bands)
        # Number of bands we acctually process.
        n_bands_proc = len(band_inds)
        if not band_inds:
            band_inds = range(n_bands)
        # Number of bands we acctually process.
        n_bands_proc = len(band_inds)
        # Read one block to figure out how many polarizations and channels
        # there are.
        Data = Reader.read(0,0)
        n_pol = Data.dims[1]
        n_cal = Data.dims[2]
        n_chan = Data.dims[3]
        # For now just use the session number as the key.
        separate = file_middle.split('/')[-1]
        sess_num = separate.split('_')[0]
        key = int(sess_num)
        # Allowcate memory for the outputs.
        covar = np.zeros((n_bands_proc, n_pol, n_cal, n_chan, n_chan),
                         dtype=float)
        counts = np.zeros(covar.shape, dtype=int)
        for ii in range(n_bands_proc):
            Blocks = Reader.read((), ii)
            for Data in Blocks:
                this_covar, this_counts = get_covar(Data)
                covar[ii,...] += this_covar
                this_counts[ii,...] += this_counts
        return key, covar, counts

def get_covar(Data):
    return 1.







