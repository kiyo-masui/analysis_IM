"""Measure the noise parameters of the data."""

import shelve
import multiprocessing as mp

import scipy as sp
import numpy.ma as ma

import noise_power as npow
from scipy import optimize
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
        kiyopy.utils.mkparents(params['output_root'] + 
                               params['output_filename'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        output_fname = params['output_root'] + params["output_filename"]
        out_db = shelve.open(output_fname)
        file_middles = params['file_middles']
        n_files = len(file_middles)
        
        n_new = nprocesses-1  # How many new processes to spawn at once.
        if n_new > 0:
            # Loop over files and spawn processes to deal with them, but make 
            # sure that only n_new processes are going at once.
            process_list = range(n_new)
            pipe_list = range(n_new)
            for ii in xrange(n_files + n_new) :
                if ii >= n_new :
                    out_db[file_middles[ii-n_new]] = pipe_list[ii%n_new].recv()
                    process_list[ii%n_new].join()
                    if process_list[ii%n_new].exitcode != 0 : 
                        raise RuntimeError("A thread failed with exit code: "
                                        + str(process_list[ii%n_new].exitcode))
                if ii < n_files :
                    input_fname = (params['input_root'] + file_middles[ii] +
                               params['input_end'])
                    Here, Far = mp.Pipe()
                    pipe_list[ii%n_new] = Here
                    process_list[ii%n_new] = mp.Process(
                        target=self.process_file, args=(input_fname, Far))
                    process_list[ii%n_new].start()
        else :
            for middle in file_middles:
                input_fname = (params['input_root'] + middle +
                               params['input_end'])
                out_db[middle] = self.process_file(input_fname)
        out_db.close()
        if self.feedback > 1 :
            print ("Wrote noise parameters to file: " 
                   + utils.abbreviate_file_path(output_fname))

    def process_file(self, file_name, Pipe=None) :
        
        try :
            model = self.params["model"]
            Reader = core.fitsGBT.Reader(file_name, feedback=self.feedback)
            Blocks = Reader.read()
            if model == "variance" :
                parameters = get_var(Blocks)
            elif model == "correlated_overf":
                parameters = get_correlated_overf(Blocks, 1.0, window='hanning')
            else :
                raise ValueError("Invalid noise model: " + model)
            if Pipe:
                Pipe.send(parameters)
            else:
                return parameters
        except :
            if Pipe:
                Pipe.send(-1)
            raise
        

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

def get_correlated_overf(Blocks, f_0=1.0, window='hanning') :
    """Measures noise parameters of a set of scans assuming correlated 1/f.

    Fits model to full f,f' cross power spectrum matrix.  Ignores the fact that
    the entries of this matrix are horribly correlated.
    
    Parameters
    ----------
    """
    
    # ---- First measre the noise power spectrum from the data. ----
    power_mat, window_function, dt, channel_means = \
            npow.full_power_mat(Blocks, window=window, deconvolve=False)
    # Check that we got finite answers.  If we didn't, need to update the code.
    if not sp.alltrue(sp.isfinite(power_mat)) :
        msg = ("Non finite power spectrum calculated.  Offending data in "
               "file starting with scan %d." % (Blocks[0].field['SCAN']))
        raise ce.DataError(msg)
    n_time = power_mat.shape[0]
    frequency = npow.ps_freq_axis(dt, n_time)
    power_mat = npow.prune_power(power_mat, 0)
    power_mat = npow.make_power_physical_units(power_mat, dt)
    # Discard the DC 0 frequency.
    frequency = frequency[1:]
    power_mat = power_mat[1:,...]
    n_chan = power_mat.shape[-1]
    n_f = len(frequency)
    
    # ---- Now fit to the measured spectrum. ----
    # First fit to the thermal free correlated part (cross terms).
    # Take mean excluding auto correlation terms.  Also take mean of the window
    # function so we can convolve the theory (this is valid since the
    # convolution is linear in the window function).
    correlated_part = sp.sum(sp.sum(power_mat, -1), -1)
    mean_window = sp.sum(sp.sum(window_function, -1), -1)
    for ii in xrange(n_chan) :
        correlated_part -= power_mat[:, ii, ii]
        mean_window -= window_function[:, ii, ii]
    correlated_part /= (n_chan - 1)*n_chan
    mean_window /= (n_chan - 1)*n_chan
    # Fit power law to this.
    def over_f_spec(params) :
        spec = npow.overf_power_spectrum(params[0], params[1],
                                                f_0, dt, n_time)
        spec = npow.convolve_power(spec, mean_window, 0)
        spec = npow.prune_power(spec)
        spec = spec[1:].real
        return spec

    def correlated_part_residuals(params):
        return (correlated_part - over_f_spec(params))/weights

    # Initial parameter guess.
    over_f_params = sp.zeros(2)
    over_f_params[0] = sp.mean(correlated_part * f_0 / frequency)
    over_f_params[1] = -1.0
    # Initial weights.
    weights = abs(correlated_part)
    # Perform fit iteratively, updating the weights.
    for ii in range(6) :
        # Memory to eliminate oscillations.  Seems to be nessisary.
        weights = 2*sp.sqrt(abs(over_f_spec(over_f_params)*weights/2))    
        over_f_params, ier = sp.optimize.leastsq(correlated_part_residuals, 
                                                 over_f_params)
        # XXX
        #plt.figure()
        #plt.loglog(frequency, correlated_part, 'b.')
        #plt.loglog(frequency, -correlated_part, 'r.')
        #plt.loglog(frequency, over_f_spec(over_f_params), 'g')
    #print over_f_params
    #plt.show() 
    # XXX
    
    # Now fit to the thermal parts one channel at a time.
    # Subtract the correlated part out of the channel autocorrelation data.
    thermal_part = power_mat
    thermal_part.shape = (n_f, n_chan*n_chan)
    thermal_part = thermal_part[:,::n_chan+1]
    thermal_part -= over_f_spec(over_f_params)[:,None]
    # Do the fit. This is just a linear fit to the diagonal parts of the power
    # matrix.  There is no need to window the theory since for a normalized
    # window, there is not affect on a constant (would of course affect the
    # weight matrix, but that's too complicated).    
    thermal_params = sp.mean(thermal_part, 0)
    for ii in xrange(3):
        weights = 2*(over_f_spec(over_f_params)[:,None] + thermal_params)
        thermal_params = (sp.sum(thermal_part/weights**2, 0)
                          / sp.sum(1.0/weights**2, 0))
    # Unpack and repack answers.
    over_f = (over_f_params[0], over_f_params[1], f_0)

    return thermal_params, over_f


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Measure(str(sys.argv[1])).execute()
               
