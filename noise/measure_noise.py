"""Measure the noise parameters of the data."""

import shelve
import multiprocessing as mp
import time as time_module

import scipy as sp
import scipy.linalg as linalg
import numpy.ma as ma

import noise_power as npow
from scipy import optimize
from kiyopy import parse_ini
import kiyopy.pickle_method
import kiyopy.utils
import kiyopy.custom_exceptions as ce
import core.fitsGBT
import utils.misc
from map.constants import T_infinity, T_small

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
               "IFs" : (),
               "save_spectra_plots" : False,
               # What parameters to measure.
               "parameters" : ["channel_var"],
               "time_block" : "scan"
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
                    Here, Far = mp.Pipe()
                    pipe_list[ii%n_new] = Here
                    process_list[ii%n_new] = mp.Process(
                        target=self.process_file, args=(file_middles[ii], Far))
                    process_list[ii%n_new].start()
        else :
            for middle in file_middles:
                out_db[middle] = self.process_file(middle)
        out_db.close()
        if self.feedback > 1 :
            print ("Wrote noise parameters to file: " 
                   + kiyopy.utils.abbreviate_file_path(output_fname))

    def process_file(self, file_middle, Pipe=None) :
        
        try :
            # Open a data file.
            params = self.params
            file_name = (params['input_root'] + file_middle
                         + params['input_end'])
            band_inds = params["IFs"]
            parameter_names = params["parameters"]
            if params["time_block"] == "scan":
                split_scans = True
            elif params["time_block"] == "file":
                split_scans = False
            else:
                raise ValueError("time_block must be 'scan' or 'file'.")
            Reader = core.fitsGBT.Reader(file_name, feedback=self.feedback)
            n_bands = len(Reader.IF_set)
            if not band_inds:
                band_inds = range(n_bands)
            # Number of bands we acctually process.
            n_bands_proc = len(band_inds)
            # Read one block to figure out how many polarizations there are.
            Data = Reader.read(0,0)
            pols = Data.field["CRVAL4"]
            n_pols = len(pols)
            # Initialize a figure that for saving the power spectra.
            if params["save_spectra_plots"]:
                h = plt.figure(figsize=(5.*n_pols, 3.5*n_bands_proc))
                h.add_subplot(n_bands_proc, n_pols, 1)
                # We will store the current subplot as an atribute of the
                # figure.
                h.current_subplot = (n_bands_proc, n_pols, 0)
            measured_parameters = {}
            band_centres = []
            for ii in range(n_bands):
                if ii in band_inds:
                    Blocks = Reader.read(params["scans"], ii)
                    Blocks[0].calc_freq()
                    n_chan = Blocks[0].dims[-1]
                    band = (int(round(Blocks[0].freq[n_chan//2]/1e6)))
                    band_centres.append(band)
                    measured_parameters[band] = measure_noise_parameters(
                            Blocks, parameter_names, split_scans=split_scans, 
                            plots=params["save_spectra_plots"])
            if params["save_spectra_plots"]:
                # Set all the plot lables.
                for ii in range(n_bands_proc):
                    for jj in range(n_pols):
                        a = h.add_subplot(n_bands_proc, n_pols, 
                                          ii * n_pols + jj + 1)
                        a.set_xlabel("frequency (Hz)")
                        a.set_ylabel("power (K^2)")
                        pol_str = utils.misc.polint2str(pols[jj])
                        band = band_centres[ii]
                        a.set_title("Band %dMHz, polarization " % band
                                    + pol_str)
                        a.autoscale_view(tight=True)
                h.subplots_adjust(hspace=0.4, left=0.2)
                fig_f_name = params['output_root'] + file_middle + ".pdf"
                kiyopy.utils.mkparents(fig_f_name)
                h.savefig(fig_f_name)
                plt.close(h.number)
            if Pipe:
                Pipe.send(measured_parameters)
            else:
                return measured_parameters
        except :
            if Pipe:
                Pipe.send(-1)
            raise

def measure_noise_parameters(Blocks, parameters, split_scans=False, 
                             plots=False):
    """Given a set of data blocks, measure noise parameters.

    Measurement done for all polarizations but only the first cal state.
    """
    
    # Initialize the output.
    out_parameters = {}
    # Calculate the full correlated power spectrum.
    power_mat, window_function, dt, channel_means = npow.full_power_mat(
            Blocks, window="hanning", deconvolve=False, n_time=-1.05,
            normalize=False, split_scans=split_scans, subtract_slope=False)
    # This shouldn't be nessisary, since I've tried to keep things finite in
    # the above function.  However, leave it in for now just in case.
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
        this_pol_power = power_mat[:,ii,cal_ind,:,:]
        this_pol_window = window_function[:,ii,cal_ind,:,:]
        this_pol = Blocks[0].field['CRVAL4'][ii]
        this_pol_parameters = {}
        # If we are plotting, activate the subplot for this polarization and
        # this band.
        if plots:
            h = plt.gcf()
            # The last subplot should be hidden in the figure object.
            current_subplot = h.current_subplot
            current_subplot = current_subplot[:2] + (current_subplot[2] + 1,)
            h.current_subplot = current_subplot
        # Now figure out what we want to measure and measure it.
        if "channel_var" in parameters:
            power_diag = this_pol_power.view()
            power_diag.shape = (n_f, n_chan**2)
            power_diag = power_diag[:,::n_chan + 1].real
            window_function_diag = this_pol_window.view()
            window_function_diag.shape = (n_time, n_chan**2)
            window_function_diag = window_function_diag[:,::n_chan + 1]
            # Integral of the power spectrum from -BW to BW.
            channel_var = sp.mean(power_diag, 0) / dt
            # Normalize for the window.
            norms = sp.mean(window_function_diag, 0).real
            bad_inds = norms < 10./n_time
            norms[bad_inds] = 1
            channel_var /= norms
            # If a channel is completly masked Deweight it by giving a high
            # variance
            channel_var[bad_inds] = T_infinity**2
            this_pol_parameters["channel_var"] = channel_var
        for noise_model in parameters:
            if noise_model[:18] == "freq_modes_over_f_":
                n_modes = int(noise_model[18:])
                this_pol_parameters[noise_model] = \
                        get_freq_modes_over_f(this_pol_power, this_pol_window,
                                              frequency, n_modes, plots=plots)
        out_parameters[this_pol] = this_pol_parameters
    return out_parameters

def get_freq_modes_over_f(power_mat, window_function, frequency, n_modes,
                          plots=False):
    """Fines the most correlated frequency modes and fits thier noise."""
    
    n_f = len(frequency)
    d_f = sp.mean(sp.diff(frequency))
    dt = 1. / 2. / frequency[-1]
    n_chan = power_mat.shape[-1]
    n_time = window_function.shape[0]
    # The threshold for assuming there isn't enough data to measure anything.
    no_data_thres = 10./n_time
    # Initialize the dictionary that will hold all the parameters.
    output_params = {}
    # First take the low frequency part of the spetrum matrix and average over
    # enough bins to get a well conditioned matrix.
    low_f_mat = sp.mean(power_mat[:4*n_chan,:,:].real, 0)
    # Factor the matrix to get the most correlated modes.
    e, v = linalg.eigh(low_f_mat)
    # Make sure they are sorted.
    if not sp.alltrue(sp.diff(e) >= 0):
        raise RuntimeError("Eigenvalues not sorted")
    # Power matrix striped of the biggest modes.
    reduced_power = sp.copy(power_mat)
    mode_list = []
    # Solve for the spectra of these modes.
    for ii in range(n_modes):
        this_mode_params = {}
        # Get power spectrum and window function for this mode.
        mode = v[:,-1 - ii]
        mode_power = sp.sum(mode * power_mat.real, -1)
        mode_power = sp.sum(mode * mode_power, -1)
        mode_window = sp.sum(mode[:,None]**2 * window_function, 1)
        mode_window = sp.sum(mode_window * mode[None,:]**2, 1)
        # Protect against no data.
        if sp.mean(mode_window).real < no_data_thres:
            this_mode_params['amplitude'] = 0.
            this_mode_params['index'] = 0.
            this_mode_params['f_0'] = 1.
            this_mode_params['thermal'] = T_infinity**2
        else:
            # Fit the spectrum.
            p = fit_overf_const(mode_power, mode_window, frequency)
            # Put all the parameters we measured into the output.
            this_mode_params['amplitude'] = p[0]
            this_mode_params['index'] = p[1]
            this_mode_params['f_0'] = p[2]
            this_mode_params['thermal'] = p[3]
        this_mode_params['mode'] = mode
        output_params['over_f_mode_' + str(ii)] = this_mode_params
        # Remove the mode from the power matrix.
        tmp_amp = sp.sum(reduced_power * mode, -1)
        reduced_power -= tmp_amp[:,:,None] * mode
        tmp_amp = sp.sum(reduced_power * mode[:,None], -2)
        reduced_power -= tmp_amp[:,None,:] * mode[:,None]
        tmp_amp = sp.sum(tmp_amp * mode, -1)
        reduced_power += tmp_amp[:,None,None] * mode[:,None] * mode
        mode_list.append(mode)
    # Initialize the compensation matrix, that will be used to restore thermal
    # noise that gets subtracted out.  See Jan 29, Feb 17th, 2012 of Kiyo's
    # notes.
    compensation = sp.eye(n_chan, dtype=float)
    for mode1 in mode_list:
        compensation.flat[::n_chan + 1] -= 2 * mode1**2
        for mode2 in mode_list:
            mode_prod = mode1 * mode2
            compensation += mode_prod[:,None] * mode_prod[None,:]
    # Now that we've striped the noisiest modes, measure the auto power
    # spectrum, averaged over channels.
    auto_spec_mean = reduced_power.view()
    auto_spec_mean.shape = (n_f, n_chan**2)
    auto_spec_mean = auto_spec_mean[:,::n_chan + 1].real
    auto_spec_mean = sp.mean(auto_spec_mean, -1)
    diag_window = window_function.view()
    diag_window.shape = (n_time, n_chan**2)
    diag_window = diag_window[:,::n_chan + 1]
    auto_spec_window = sp.mean(diag_window, -1)
    if sp.mean(auto_spec_window).real < no_data_thres:
        auto_cross_over = 0.
        auto_index = 0.
        auto_thermal = 0
    else:
        auto_spec_params = fit_overf_const(auto_spec_mean, auto_spec_window,
                                           frequency)
        auto_thermal = auto_spec_params[3]
        if (auto_spec_params[0] <= 0 or auto_spec_params[3] <= 0 or
            auto_spec_params[1] > -0.599):
            auto_cross_over = 0.
            auto_index = 0.
        else:
            auto_index = auto_spec_params[1]
            auto_cross_over = auto_spec_params[2] * (auto_spec_params[0]
                                     / auto_spec_params[3])**(-1./auto_index)
            #if auto_cross_over < d_f:
            #    auto_index = 0.
            #    auto_cross_over = 0.
    # Plot the mean auto spectrum if desired.
    if plots:
        h = plt.gcf()
        a = h.add_subplot(*h.current_subplot)
        norm = sp.mean(auto_spec_window).real
        auto_plot = auto_spec_mean / norm
        plotable = auto_plot > 0
        lines = a.loglog(frequency[plotable], auto_plot[plotable])
        c = lines[-1].get_color()
        # And plot the fit in a light color.
        if auto_cross_over > d_f / 4.:
            spec = npow.overf_power_spectrum(auto_thermal, auto_index, 
                                             auto_cross_over, dt, n_time)
        else:
            spec = sp.zeros(n_time, dtype=float)
        spec += auto_thermal
        spec[0] = 0
        spec = npow.convolve_power(spec, auto_spec_window)
        spec = npow.prune_power(spec)
        spec = spec[1:].real
        if norm > no_data_thres:
            spec /= norm
        plotable = spec > 0
        a.loglog(frequency[plotable], spec[plotable], c=c, alpha=0.4,
                 linestyle=':')
    output_params['all_channel_index'] = auto_index
    output_params['all_channel_corner_f'] = auto_cross_over
    # Finally measure the thermal part of the noise in each channel.
    cross_over_ind = sp.digitize([auto_cross_over * 4], frequency)[0]
    cross_over_ind = max(cross_over_ind, n_f // 2)
    cross_over_ind = min(cross_over_ind, int(9. * n_f / 10.))
    thermal = reduced_power[cross_over_ind:,:,:].real
    n_high_f = thermal.shape[0]
    thermal.shape = (n_high_f, n_chan**2)
    thermal = sp.mean(thermal[:,::n_chan + 1], 0)
    thermal_norms = sp.mean(diag_window, 0).real
    bad_inds = thermal_norms < no_data_thres
    thermal_norms[bad_inds] = 1.
    # Compensate for power lost in mode subtraction.
    compensation[:,bad_inds] = 0
    compensation[bad_inds,:] = 0
    for ii in xrange(n_chan):
        if bad_inds[ii]:
            compensation[ii,ii] = 1.
    thermal = linalg.solve(compensation, thermal)
    # Normalize
    thermal /= thermal_norms
    thermal[bad_inds] = T_infinity**2
    # Occationally the compensation fails horribly on a few channels.
    # When this happens, zero out the offending indices.
    thermal[thermal<0] = 0
    output_params['thermal'] = thermal
    # Now that we know what thermal is, we can subtract it out of the modes we
    # already measured.
    for ii in range(n_modes):
        mode_params = output_params['over_f_mode_' + str(ii)]
        thermal_contribution = sp.sum(mode_params['mode']**2 * thermal)
        # Subtract a maximum of 90% of the white noise to keep things positive
        # definate.
        new_white = max(mode_params['thermal'] - thermal_contribution, 
                        0.1 * mode_params['thermal'] )
        if mode_params['thermal'] < 0.5 * T_infinity**2:
            mode_params['thermal'] = new_white
    return output_params


def fit_overf_const(power, window, freq):
    """Fit $A*f**alpha + C$ to a 1D power spectrum.
    
    Power spectrum should be real and contain only positive fequencies.  Window
    function should be complex and have both positive and negitive frequencies.
    """
    
    n_f = len(freq)
    f_0 = 1.0
    dt = 1./2./freq[-1]
    n_time = len(window)
    # Minimum spectral index.
    min_index = 0.6
    if n_f != n_time // 2 - 1:
        raise RuntimeError("Power spectrum and window sizes incompatible.")
    # Make sure the window function is well behaved.
    if ((not window[0].real > 0)
        or (not abs(window[0].imag) < 1.e-6 * window[0].real)):
        raise ValueError("Badly behaved window function.")
    # Spectrum function we will fit too.
    # Instead of fitting for the index directly, we sue another parameter that
    # does not admit index >= -0.2.  This avoids a degeneracy with thermal when
    # index = 0.
    # Also force the amplitudes to be positive by squaring and give them a
    # minimum value.
    def model(params):
        a = params[0]**2 + T_small**2
        i = -(params[1]**2 + min_index)
        t = params[2]**2 + T_small**2
        spec = npow.overf_power_spectrum(a, i, f_0, dt, n_time)
        spec += t
        spec = npow.convolve_power(spec, window)
        spec = npow.prune_power(spec)
        spec = spec[1:].real
        return spec
    # Residuals function.
    def residuals(params):
        return (power - model(params))/weights
    # A function for the Jacobian matrix.
    def jacobian(params):
        a = params[0]
        i = params[1]
        t = params[2]
        # Get the frequencies, including the negitive ones.
        df = 1.0/dt/n_time
        f = sp.arange(n_time, dtype=float)
        f[n_time//2+1:] -= f[-1] + 1
        f = abs(f)*df
        # 0th mode is meaningless.  Set to unity to avoid errors.
        f[0] = 1
        # The power law part.
        p_law = (f/f_0)**(-(params[1]**2 + min_index))
        # Memory for the three derivative functions.
        spec = sp.empty((3, n_time), dtype=float)
        # The formula for the following derivatives derived in Kiyo's notes,
        # Nov. 28, 2011.
        # Derivative with respect to the amplitude parameter.
        spec[0,:] = -2.0 * a * p_law
        # Derivative with respect to the index parameter.
        spec[1,:] = (a**2 + T_small**2) * p_law * sp.log(f/f_0) * 2.0 * i
        # Derivative with respect to the thermal parameter.
        spec[2,:] = -2.0 * t
        # Get rid of the mean mode.
        spec[:,0] = 0
        # Convolve with the window function.
        spec = npow.convolve_power(spec, window[None,:], -1)
        # Prune to the same size as the data.
        spec = npow.prune_power(spec, -1)
        spec = spec[:,1:].real
        return spec/weights
    # Get good initial guesses for the parameters.  It is extremely important
    # to do well here.
    params = sp.empty(3, dtype=float)
    # First guess the thermal level by taking the mean at high frequencies.
    norm = n_time / sp.sum(window.real)
    #print norm
    params[2] = sp.sqrt(abs(sp.mean(power[-n_f//10:]) * norm))
    params[1] = sp.sqrt(1 - min_index) # Corresponds to index = -1
    params[0] = sp.sqrt(abs(sp.mean((power * freq)[:n_f//10]) * norm))
    old_weights = abs(power)
    # Iteratively fit then update the weights.
    for ii in range(4):
        new_weights = abs(model(params))
        weights = old_weights + new_weights
        old_weights = new_weights
        # XXX
        #plt.figure()
        #plt.loglog(freq, power)
        #plt.loglog(freq, model(params))
        #plt.figure()
        #plt.semilogx(freq, residuals(params))
        #print params
        ###
        params, cov_x, info, mesg, ier = sp.optimize.leastsq(residuals, 
                    params, Dfun=jacobian, col_deriv=True, xtol=0.001,
                    full_output=True)
    #plt.figure()
    #plt.loglog(freq, power)
    #plt.loglog(freq, model(params))
    #print params
    #plt.show()
    # Check that a solution was found.
    #if ier not in (1, 2, 3, 4):
    #    raise RuntimeError("Could not find a solution. " + repr(params))
    # Unpack results and return.
    amp = params[0]**2 + T_small**2
    index = -(params[1]**2 + min_index)
    thermal = params[2]**2 + T_small**2
    return amp, index, f_0, thermal



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Measure(str(sys.argv[1])).execute()
               
