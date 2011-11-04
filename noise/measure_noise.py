"""Measure the noise parameters of the data."""

import shelve
import multiprocessing as mp
import time as time_module

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
               "file_middles" : ("testfile_guppi_combined",),
               "input_end" : ".fits",
               "output_root" : "./",
               "output_filename" : "noise_parameters.shelve",
               "scans" : (),
               "bands" : (),
               # What parameters to measure.
               "parameters" : ["channel_var", "mean_over_f"]
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
            band_inds = self.params["bands"]
            parameter_names = self.params["parameters"]
            Reader = core.fitsGBT.Reader(file_name, feedback=self.feedback)
            n_bands = len(Reader.IF_set)
            if not band_inds:
                band_inds = range(n_bands)
            measured_parameters = {}
            for ii in range(n_bands):
                if ii in band_inds:
                    Blocks = Reader.read(self.params["scans"], ii)
                    Blocks[0].calc_freq()
                    n_chan = Blocks[0].dims[-1]
                    band = (int(round(Blocks[0].freq[n_chan//2]/1e6)))
                    measured_parameters[band] = measure_noise_parameters(
                            Blocks, parameter_names)
            if Pipe:
                Pipe.send(measured_parameters)
            else:
                return measured_parameters
        except :
            if Pipe:
                Pipe.send(-1)
            raise

def measure_noise_parameters(Blocks, parameters):
    """Given a set of data blocks, measure noise parameters.

    Measurement done for all polarizations but only the first cal state.
    """
    
    # Initialize the output.
    out_parameters = {}
    # Calculate the full correlated power spectrum.
    power_mat, window_function, dt, channel_means = npow.full_power_mat(
            Blocks, window="hanning", deconvolve=False, n_time=-1.05)
    # This shouldn't be nessisary, since I've tried to keep thing finite in the
    # above function.  However, leave it in for now just in case.
    if not sp.alltrue(sp.isfinite(power_mat)) :
        msg = ("Non finite power spectrum calculated.  Offending data in "
               "file starting with scan %d." % (Blocks[0].field['SCAN']))
        raise ce.DataError(msg)
    # Get frequency axis and do unit conversions.
    n_time = power_mat.shape[0]
    n_chan = power_mat.shape[-1]
    frequency = npow.ps_freq_axis(dt, n_time)
    power_mat = npow.prune_power(power_mat, 0)
    power_mat = npow.make_power_physical_units(power_mat, dt)
    # Discard the mean mode.
    frequency = frequency[1:]
    power_mat = power_mat[1:,...]
    n_f = len(frequency)
    # Loop over polarizations.
    cal_ind = 0
    n_pols = power_mat.shape[1]
    for ii in range(n_pols):
        this_pol_power = power_mat[:,ii,cal_ind,:]
        this_pol_window = window_function[:,ii,cal_ind,:]
        this_pol = Blocks[0].field['CRVAL4'][ii]
        this_pol_parameters = {}
        # Now figure out what we want to measure and measure it.
        if "channel_var" in parameters:
            power_diag = this_pol_power.view()
            power_diag.shape = (n_f, n_chan**2)
            power_diag = power_diag[:,::n_chan + 1].real
            this_pol_parameters["channel_var"] = sp.mean(power_diag, 0)/(2*dt)
        if "mean_over_f" in parameters:
            this_pol_parameters["mean_over_f"] = get_mean_over_f(
                    this_pol_power, this_pol_window, frequency)
        out_parameters[this_pol] = this_pol_parameters
    return out_parameters

def get_mean_over_f(power_mat, window_function, frequency) :
    """Measures noise parameters of a set of scans assuming correlated 1/f.

    Fits model to full f,f' cross power spectrum matrix.  Ignores the fact that
    the entries of this matrix are horribly correlated.
    
    Parameters
    ----------
    power_mat : array
        Pruned of egitive frequencies and 0 frequency.
    """
    
    n_f = len(frequency)
    dt = 1.0 / (frequency[-1] * 2)
    n_chan = power_mat.shape[-1]
    f_0 = 1.0
    n_time = window_function.shape[0]
    # First find the auto correlation part.
    auto_corr = power_mat.view()
    auto_corr.shape = (n_f, n_chan*n_chan)
    auto_corr = auto_corr[:,::n_chan + 1].real    
    # First fit to the thermal free correlated part (cross terms).
    # Take mean excluding auto correlation terms.
    correlated_part = sp.sum(sp.sum(power_mat, -1), -1)
    mean_window = sp.sum(sp.sum(window_function, -1), -1)
    for ii in xrange(n_chan) :
        correlated_part -= power_mat[:, ii, ii]
        mean_window -= window_function[:, ii, ii]
    correlated_part /= (n_chan - 1)*n_chan
    correlated_part = correlated_part.real
    mean_window /= (n_chan - 1)*n_chan
    # Fit power law to this.
    def over_f_spec(params, window) :
        spec = npow.overf_power_spectrum(params[0], params[1],
                                                f_0, dt, n_time)
        spec = npow.convolve_power(spec, window, 0)
        spec = npow.prune_power(spec)
        spec = spec[1:].real
        return spec
    def correlated_part_residuals(params):
        return (correlated_part - over_f_spec(params, mean_window))/weights
    # Initial parameter guesses.
    over_f_params = sp.zeros(2)
    over_f_params[0] = sp.mean(correlated_part * f_0 / frequency)
    over_f_params[1] = -1.0
    thermal_params = sp.mean(auto_corr, 0)
    # Initial weights.
    old_weights = abs(correlated_part)
    old_weights[old_weights < 1e-16] = 1
    old_thermal_weights = abs(auto_corr)
    old_thermal_weights[old_thermal_weights < 1e-16] = 1
    # Perform fit iteratively, updating the weights.  We do a two step fit,
    # doing the thermal part after the correlated part (much faster).
    auto_over_f_specs = sp.empty((n_f, n_chan))
    for ii in range(6) :
        # Update the weights for the correlated part.
        # Memory to eliminate oscillations.  Seems to be nessisary.
        new_weights = (abs(over_f_spec(over_f_params, mean_window)) +
                       abs(sp.mean(thermal_params)))
        new_weights[new_weights < 1e-16] = 1
        weights = 2*sp.sqrt(old_weights * new_weights)
        old_weights = new_weights
        # Fit for the correlated part.
        over_f_params, ier = sp.optimize.leastsq(correlated_part_residuals, 
                                                 over_f_params)
        # Update the weights for the thermal part.
        for ii in xrange(n_chan):
            auto_over_f_specs[:,ii] = over_f_spec(over_f_params,
                                                  window_function[:,ii,ii])
        new_thermal_weights = auto_over_f_specs + thermal_params
        new_thermal_weights[new_thermal_weights < 1e-16] = 1
        thermal_weights = 2*sp.sqrt(old_thermal_weights * new_thermal_weights)
        old_thermal_weights = new_thermal_weights
        # Fit for the thermal part.
        thermal_part = auto_corr - auto_over_f_specs
        thermal_params = (sp.sum(thermal_part/thermal_weights**2, 0)
                          / sp.sum(1.0/thermal_weights**2, 0))
        # XXX
        #plt.figure()
        #plt.loglog(frequency, correlated_part, 'b.')
        #plt.loglog(frequency, -correlated_part, 'r.')
        #plt.loglog(frequency, over_f_spec(over_f_params), 'g')
        #plt.loglog(frequency, sp.mean(thermal_params)*sp.ones(n_f), 'y')
    #print over_f_params, thermal_params
    #plt.show() 
    # XXX
    # Unpack and repack answers.
    over_f = (over_f_params[0], over_f_params[1], f_0)
    return thermal_params, over_f



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Measure(str(sys.argv[1])).execute()
               
