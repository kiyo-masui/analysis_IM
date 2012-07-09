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
#import matplotlib.pyplot as plt


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
        # Store the covariance and counts session by session in a dictionary.
        covar_dict = {}
        counts_dict = {}
        # Also store the frequency axis.
        freq_dict = {}
        # Loop though all the files and accumulate the covariance and the
        # counts.  Files in the same session are summed together.
        for middle in file_middles:
            key, covar, counts, freq = self.process_file(middle)
            if covar_dict.has_key(key):
                if covar_dict[key].shape != covar.shape:
                    msg = ("All data needs to have the same band and"
                           "polarization structure.")
                    raise ce.DataError(msg)
                covar_dict[key] += covar
                counts_dict[key] += counts
                if not np.allclose(freq_dict[key], freq):
                    raise ce.DataError("Frequency structure not consistant.")
            else:
                covar_dict[key] = covar
                counts_dict[key] = counts
                freq_dict[key] = freq
        # Now that we have the covariance, factor it into eigen-vectors and
        # store it in a data base.
        output_fname = params['output_root'] + params["output_filename"]        
        out_db = shelve.open(output_fname)
        # Loop through all the data divisions and processes them one at a
        # time.
        for key in covar_dict.iterkeys():
            covar = covar_dict[key]
            counts = counts_dict[key]
            # Normalize.
            counts[counts==0] = 1
            covar /= counts
            # Loop to each matrix, decompose it and save it.
            eigen_vects = np.empty_like(covar)
            for band_ii in range(covar.shape[0]):
                for pol_jj in range(covar.shape[1]):
                    for cal_kk in range(covar.shape[2]):
                        # Factor
                        h, v = linalg.eigh(covar[band_ii,pol_jj,cal_kk])
                        eigen_vects[band_ii,pol_jj,cal_kk] = v
                        #plt.semilogy(h, '.')
                        #plt.figure()
                        #for ii in range(1,5):
                        #    plt.plot(v[:,-ii])
                        #plt.show()
            out_db[key + '.vects'] = eigen_vects
            out_db[key + '.freq'] = freq_dict[key]
        out_db.close()


    def process_file(self, middle):
        params = self.params
        file_name = (params['input_root'] + middle
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
        # Get the key that will group this file with other files.
        key = get_key(middle)
        # Read one block to figure out how many polarizations and channels
        # there are.
        Data = Reader.read(0,0)
        n_pol = Data.dims[1]
        n_cal = Data.dims[2]
        n_chan = Data.dims[3]
        # Allowcate memory for the outputs.
        covar = np.zeros((n_bands_proc, n_pol, n_cal, n_chan, n_chan),
                         dtype=float)
        counts = np.zeros(covar.shape, dtype=int)
        freq = np.empty((n_bands_proc, n_chan))
        for ii in range(n_bands_proc):
            Blocks = Reader.read((), ii)
            Blocks[0].calc_freq()
            freq[ii,:] = Blocks[0].freq
            for Data in Blocks:
                this_covar, this_counts = get_covar(Data)
                covar[ii,...] += this_covar
                counts[ii,...] += this_counts
        return key, covar, counts, freq


def get_key(middle):
    # For now just use the session number as the key.
    separate = middle.split('/')[-1]
    sess_num = separate.split('_')[0]
    key = sess_num
    return key

def get_covar(Data):
    
    data = Data.data.filled(0)
    mask = ma.getmaskarray(Data.data)
    weights = np.logical_not(mask)
    # Works much better if I take out the first moment.  Just remember to take
    # it out data when I take out the modes as well.
    data -= ma.mean(data, 0)

    #for ii in range(weights.shape[0]):
    #    if np.any(weights[ii,0,0,:] == 0):
    #        weights[ii,0,0,:] = 0
    
    data *= weights

    covar = np.sum(data[...,None,:] * data[...,:,None], 0)
    counts = np.sum(weights[...,None,:] * weights[...,:,None], 0)

    return covar, counts








