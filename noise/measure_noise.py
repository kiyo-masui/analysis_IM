"""Measure the noise parameters of the data."""

import shelve

import scipy as sp
import numpy.ma as ma

import noise_power as np
from scipy import optimize
from kiyopy import parse_ini, utils
import kiyopy.utils
import core.fitsGBT

# XXX for testing.
import matplotlib.pyplot as plt


params_init = {
               # IO.
               "input_root" : "./testdata/",
               "file_middles" : ("testfile_GBTfits",),
               "input_end" : ".fits",
               "output_root" : "./",
               "output_filename" : "noise_parameters.shelf",
               # Algorithm.
               "model" : "variance"
               }

prefix = 'mn_'

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
        model = params["model"]
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        output_fname = params["output_root"] + params["output_filename"]
        out_db = shelve.open(output_fname)
        # Loop over files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            Reader = core.fitsGBT.Reader(input_fname, feedback=self.feedback)
            Blocks = Reader.read()
            if model == "variance" :
                parameters = get_var(Blocks)
            else :
                raise ValueError("Invalid noise model: " + model)
            out_db["file_middle"] = parameters

        out_db.close()
        if self.feedback > 1 :
            print ("Wrote noise parameters to file: " 
                   + utils.abbreviate_file_path(output_fname))

def get_var(Blocks) :
    """Measures the variance of a set of scans."""
    
    # These all become arrays on the first iteration.
    var = 0.0
    mean = 0.0
    counts = 0
    for Data in Blocks:
        var += ma.sum(Data.data**2, 0).filled(0)
        mean += ma.sum(Data.data, 0).filled(0)
        counts += ma.count(Data.data, 0)
    # If we didn't get at least 5 good hits, throw aways the
    # scan.
    counts[counts < 5] = -1
    var = var/counts - (mean/counts)**2
    var[counts < 5] = 1.0e10
    
    return var

def get_correlated_overf(Blocks, f_0=1.0) :
    """Measures noise parameters of a set of scans assuming correlated 1/f.
    """
    
    full_data, mask, dt = np.make_masked_time_stream(Blocks)
    n_time = full_data.shape[0]
    # Broadcast to shape such that all pairs of channels are calculated.
    full_data1 = full_data[:,0,0,:,None]
    full_data2 = full_data[:,0,0,None,:]
    mask1 = mask[:,0,0,:,None]
    mask2 = mask[:,0,0,None,:]
    # Frequencies for the power spectrum (Hz)
    frequency = np.ps_freq_axis(dt, n_time)
    n_chan = full_data.shape[-1]
    # Loop to do only a bit of the calculation at a time.  Reduces
    # memory use by a factor of a few.
    power_mat = sp.empty((n_time//2, n_chan, n_chan),
                              dtype=float)
    for ii in xrange(n_chan):
        # `out` argument must be a view, not a copy.
        np.windowed_power(full_data1[:,[ii],:], mask1[:,[ii],:], full_data2,
                       mask2, axis=0, out=(power_mat[:,ii,:])[:,None,:])
    # Discard the DC 0 frequency.
    frequency = frequency[1:]
    power_mat = power_mat[1:,...]

    # Now that we have the power spectrum, fit the noise model to it.
    # Make the residual function.
    n_f = len(frequency)
    
    def correlated_overf_spec(params) :
        amp = params[-2]
        index = params[-1]
        chan_thermal = params[:-2] # len = n_chan.
        power_theory = sp.zeros_like(power_mat)
        # Add the thermal compontent to the diagonal only (channel1 =
        # channel2).
        chan_diag = power_theory.reshape((n_f, n_chan*n_chan))[:,::n_chan+1]
        chan_diag += chan_thermal
        # Add the 1/f component to all channel combinations.
        over_f_spec = amp*(frequency/f_0)**-index
        power_theory += over_f_spec[:,None,None]
        return power_theory
    
    def get_residuals(params) :
        # Return residuals.
        power_theory = correlated_overf_spec(params)
        # The windowed power spectrum really has a full covariance.  Could use
        # that as weights...
        residuals = (power_mat - power_theory)/weights
        return residuals.flat
    
    # On first fit, use uniform errors. On each following iteration, get the
    # errors from the best fit solution.
    weights = 1.0
    x = sp.ones(n_chan + 2)
    for ii in range(3) :
        x, ier = sp.optimize.leastsq(get_residuals, x)
        weights = 2*correlated_overf_spec(x)

    if False :
        #p = sp.empty(n_chan+2)
        #p[:n_chan] = 1.0 + 1.0/n_chan*sp.arange(n_chan)
        #p[-2] = 2.0
        #p[-1] = 1.5
        p = x
        plt.loglog(frequency, power_mat[:,5,5])
        plt.loglog(frequency, correlated_overf_spec(p)[:,5,5])
        plt.figure()
        plt.loglog(frequency, power_mat[:,5,4])
        plt.loglog(frequency, correlated_overf_spec(p)[:,5,4])
        plt.show()
        
        
    # Unpack answers.
    over_f = (x[-2], x[-1])
    thermal = x[:-2]
    
    return thermal, over_f


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Measure(str(sys.argv[1])).execute()
               
