"""Dirty map making module.

Module converts data in the time domain into noise weighted data in the map
domain, i.e. it creats the dirty map.  Module also contains many utilities like
the pointing operator (`Pointing`) and the time domain noise operator
(`Noise`).
"""

import math
import threading
from Queue import Queue
import shelve
import sys
import time as time_mod

import scipy as sp
import numpy.ma as ma
import scipy.fftpack as fft
from scipy import linalg
from scipy import interpolate
import numpy as np
#import matplotlib.pyplot as plt

import core.algebra as al
from core import fitsGBT
import utils.misc as utils
import tools
from noise import noise_power
from foreground import ts_measure
import kiyopy.custom_exceptions as ce
import kiyopy.utils
from kiyopy import parse_ini
import _mapmaker as _mapmaker_c
from constants import T_infinity, T_huge, T_large, T_medium, T_small, T_sys
from constants import f_medium, f_large

prefix ='dm_'
params_init = {# IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./testoutput_",
               # Shelve file that holds measurements of the noise.
               # What data to include from each file.
               'scans' : (),
               'IFs' : (0,),
               'polarizations' : ('I',),
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               # In pixels.
               'map_shape' : (5, 5),
               'pixel_spacing' : 0.5, # degrees
               # Interpolation between pixel points.  Options are 'nearest',
               # 'linear' and 'cubic'.
               'interpolation' : 'linear',
               # How to treat the data.
               # How much data to include at a time (scan by scan or file by
               # file)
               'time_block' : 'file',
               # Number of files to process at a time.  Increasing this number
               # reduces IO but takes more memory.
               'n_files_group' : 0,
               # What kind of frequency correlations to use.  Options are
               # 'None', 'mean' and 'measured'.
               'frequency_correlations' : 'None',
               # Ignored unless 'frequency_correlations' is 'measured'.
               'number_frequency_modes' : 1,
               # Of the above modes, how many to completely discard instead of
               # noise weighting.
               'number_frequency_modes_discard' : 0,
               # Where to get noise parameters from.
               'noise_parameter_file' : '',
               'deweight_time_mean' : True,
               'deweight_time_slope' : False,
               # If there where any foregrounds subtracted in the time stream,
               # let the noise know about it.
               'n_ts_foreground_modes' : 0,
               'ts_foreground_mode_file' : ''
               }

class DirtyMapMaker(object):
    """Dirty map maker.
    """

    def __init__(self, parameter_file_or_dict=None, feedback=2):
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback
    
    def iterate_data(self, file_middles):
        """An iterator over the input data.
        
        This can either iterate over the file names, processing the whole file
        at once, or it can iterate over both the file and the scans.
        """
    
        params = self.params
        for middle in file_middles:
            self.file_middle = middle
            fname = params["input_root"] + middle + params["input_end"]
            Reader = fitsGBT.Reader(fname, feedback=self.feedback)
            Blocks = Reader.read(self.params['scans'], self.band_ind,
                                 force_tuple=True)
            if params['time_block'] == 'scan':
                for Data in Blocks:
                    yield self.preprocess_data((Data,))
            elif params['time_block'] == 'file':
                yield self.preprocess_data(Blocks)
            else:
                msg = "time_block parameter must be 'scan' or 'file'."
                raise ValueError(msg)

    def preprocess_data(self, Blocks):
        """The method converts data to a mapmaker friendly format."""

        # On the first loop through the data, figure out the total number 
        # of time elements and check that the
        # data has expected properties (right number of channels, right band,
        # right polarization axis, etc.).
        # Also calculate the mean in each channel and subtract it off.
        nt = 0
        channel_sums = 0
        channel_sq_sums = 0
        channel_counts = 0
        for Data in Blocks:
            nt += Data.dims[0]
            if Data.dims[3] != self.n_chan:
                msg = "Data doesn't all have the same number of channels."
                raise ce.DataError(msg)
            Data.calc_freq()
            if (round(Data.freq[self.n_chan//2])
                != round(self.band_centres[self.band_ind])):
                msg = "Data doesn't all have the same band structure."
                raise ce.DataError(msg)
            if (self.pols[self.pol_ind] != Data.field['CRVAL4'][self.pol_ind]):
                msg = "Data doesn't all have the same polarizations."
                raise ce.DataError(msg)
            filled_data = Data.data[:,self.pol_ind,0,:].filled(0)
            channel_sums += sp.sum(filled_data, 0)
            channel_sq_sums += sp.sum(filled_data**2, 0)
            channel_counts += ma.count(Data.data[:,self.pol_ind,0,:], 0)
        channel_counts[channel_counts == 0] = 1
        channel_means = channel_sums / channel_counts
        channel_vars = channel_sq_sums / channel_counts
        channel_vars = channel_vars - channel_means**2
        # Get rid of any wierd channels.
        bad_channels =  (channel_vars < 1.e-4 
                        * sp.median(channel_vars[channel_counts > 5]))
        bad_channels = sp.logical_or(channel_counts < 5, bad_channels)
        channel_vars[bad_channels] = T_infinity**2
        # Store the variances in case they are needed as noise weights.
        self.channel_vars = channel_vars
        # Allocate memory for the outputs.
        time_stream = sp.empty((self.n_chan, nt), dtype=float)
        time_stream = al.make_vect(time_stream, axis_names=('freq', 'time'))
        ra = sp.empty(nt, dtype=float)
        dec = sp.empty(nt, dtype=float)
        az = sp.empty(nt, dtype=float)
        el = sp.empty(nt, dtype=float)
        time = sp.empty(nt, dtype=float)
        mask_inds_time = sp.zeros(0, dtype=int)
        mask_inds_chan = sp.zeros(0, dtype=int)
        # Now loop through and copy the data.
        tmp_time_ind = 0
        tmp_dt = []
        for Data in Blocks:
            # First deal with the data.
            this_nt = Data.dims[0]
            this_data = Data.data[:,self.pol_ind,0,:] - channel_means
            time_stream[:,tmp_time_ind:tmp_time_ind + this_nt] = \
                    this_data.filled(0).transpose()
            # Now figure out if any of the data is masked.
            mask = ma.getmaskarray(Data.data[:,self.pol_ind,0,:])
            masked = sp.where(mask)
            mask_inds_time = sp.concatenate((mask_inds_time, masked[0]), 0)
            mask_inds_chan = sp.concatenate((mask_inds_chan, masked[1]), 0)
            # Finially, copy the other fields.
            Data.calc_pointing()
            ra[tmp_time_ind:tmp_time_ind + this_nt] = Data.ra
            dec[tmp_time_ind:tmp_time_ind + this_nt] = Data.dec
            Data.calc_time()
            time[tmp_time_ind:tmp_time_ind + this_nt] = Data.time
            tmp_dt.append(abs(sp.mean(sp.diff(Data.time))))
            az[tmp_time_ind:tmp_time_ind + this_nt] = Data.field['CRVAL2']
            el[tmp_time_ind:tmp_time_ind + this_nt] = Data.field['CRVAL3']
            tmp_time_ind += this_nt
        mask_inds = (mask_inds_chan, mask_inds_time)
        # Make sure the bandwidths are the same for all Data Blocks.
        dt = sp.mean(tmp_dt)
        if not sp.allclose(tmp_dt, dt, rtol=1e-2):
            print tmp_dt, dt
            raise RuntimeError("Samplings not uniform.")
        self.BW = 1./2./dt
        # Now trim this down to exculd all data points that are outside the map
        # bounds.
        map = self.map
        time_stream, inds = trim_time_stream(time_stream, (ra, dec),
                (min(map.get_axis('ra')), min(map.get_axis('dec'))),
                (max(map.get_axis('ra')), max(map.get_axis('dec'))))
        ra = ra[inds]
        dec = dec[inds]
        time = time[inds]
        az = az[inds]
        el = el[inds]
        new_mask_inds = ([], [])
        for jj in range(len(mask_inds[0])):
            if mask_inds[1][jj] in inds:
                new_time_ind = sp.where(inds == mask_inds[1][jj])[0][0]
                new_mask_inds[0].append(mask_inds[0][jj])
                new_mask_inds[1].append(new_time_ind)
        mask_inds = (sp.asarray(new_mask_inds[0], dtype=int),
                     sp.asarray(new_mask_inds[1], dtype=int))
        return time_stream, ra, dec, az, el, time, mask_inds

    def get_noise_parameter(self, parameter_name):
        """Reads a desired noise parameter for the current data."""
        
        if not self.params['noise_parameter_file']:
            # If there is no data base file we can use the measured variance as
            # thermal.
            if parameter_name == "thermal":
                return self.channel_vars
            else :
                raise ValueError('Only thermal parameter can be measured '
                                 'with no noise parameter database.')
        else :
            # Get the data base entry that corresponds to the current file.
            noise_entry = (self.noise_params[self.file_middle]
                           [int(round(self.band_centres[self.band_ind]/1e6))]
                           [self.pols[self.pol_ind]])
            # In many cases the thermal noise is measured in units K**2/Hz.  To
            # get a variance we need to multiply by twice the bandwidth (factor
            # of 2 deals with negitive frequencies).
            BW = self.BW
            # Choose the noise model and get the corresponding parameters.
            noise_model = self.params["frequency_correlations"]
            if noise_model == "None":
                if parameter_name == "thermal":
                    if noise_entry.has_key("channel_var"):
                        return noise_entry["channel_var"]
                    else:
                        return self.channel_vars
            elif noise_model == "mean":
                noise_entry = noise_entry["mean_over_f"]
                if parameter_name == "thermal":
                    return noise_entry["thermal"] * BW * 2
                elif parameter_name == "mean_over_f":
                    over_f_params = (noise_entry["amplitude"],
                                     noise_entry["index"], noise_entry["f_0"])
                    return over_f_params
            elif noise_model == "measured":
                noise_entry = noise_entry["freq_modes_over_f_"
                                + repr(self.params['number_frequency_modes'])]
                if parameter_name == "thermal":
                    thermal = noise_entry["thermal"] * BW * 2
                    # Min is factor of 2 smaller than theoredical (since some
                    # noise will be in the frequency modes).
                    thermal_min = abs(T_sys**2 * 2. * BW / self.delta_freq / 2.)
                    thermal[thermal < thermal_min] = thermal_min
                    return thermal
                elif parameter_name[:12] == "over_f_mode_":
                    p = noise_entry[parameter_name]
                    # Factors of the band width applied in Noise object.
                    return (p["amplitude"], p['index'], p['f_0'], p['thermal'],
                            p['mode'])
                elif parameter_name == "all_channel":
                    # When specifying as a cross over, thermal is the
                    # amplitude.
                    thermal = sp.copy(noise_entry["thermal"])
                    thermal_min = abs(T_sys**2 / self.delta_freq / 2.)
                    thermal[thermal < thermal_min] = thermal_min
                    p = (thermal,)
                    p += (noise_entry["all_channel_index"],)
                    p += (noise_entry["all_channel_corner_f"],)
                    return p
        # If we got to this point, we didn't find a model-parameter match.
        msg = ("Invalid noise parameter name: " + parameter_name
               + " for noise model: " + noise_model)
        raise ValueError(msg)

    def get_ts_foreground_mode(self, mode_number):
        # Collect the variables we need from the class.
        foreground_modes_db = self.foreground_modes
        pol_ind = self.pol_ind
        cal_ind = 0
        file_middle = self.file_middle
        freq = self.map.get_axis('freq')
        # Get the data base key.
        key = ts_measure.get_key(file_middle)
        # From the database, read the data for this file.
        db_freq = foreground_modes_db[key + '.freq']
        db_vects = foreground_modes_db[key + '.vects']
        # Loop through the IF's in the data base and figure out which on is the
        # right one.
        n_bands_db = db_freq.shape[0]
        db_band = -1
        for ii in range(n_bands_db):
            if ((min(db_freq[ii,:]) <= min(freq))
                and (max(db_freq[ii,:]) >= max(freq))):
                db_band = ii
        if db_band == -1:
            raise RuntimeError("Did not find an overlapping band in foreground"
                               " data.")
        # Throw away the data we don't need.
        db_freq = db_freq[db_band,:]
        db_vects = db_vects[db_band,pol_ind,cal_ind,:,-1 - mode_number]
        # Now just need to find what range of frequencies correspond to the
        # band we are dealing with.
        n_chan_db = len(db_freq)
        n_chan = len(freq)
        chan_start_ind = -1
        for ii in range(n_chan_db - n_chan):
            if sp.allclose(freq, db_freq[ii:ii + n_chan]):
                chan_start_ind = ii
        if chan_start_ind == -1:
            raise RuntimeError("Could not find frequency range in foreground"
                               " data.")
        # Get that piece of the vector and normalize.
        vector = sp.copy(db_vects[chan_start_ind:chan_start_ind + n_chan])
        vector /= sp.sum(vector**2)
        return vector
        
    def execute(self, n_processes):
        """Driver method."""
        
        self.n_processes = n_processes
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix='mm_')
        if self.feedback > 0:
            print "Input root is: " + params["input_root"]
            print "Output root is: " + params["output_root"]
        # Set flag if we are doing independant channels.
        if params['frequency_correlations'] == 'None':
            self.uncorrelated_channels = True
            if self.feedback > 0:
                print "Treating frequency channels as being independant."
        else:
            self.uncorrelated_channels = False
        # Open the file that stores noise parameters.
        if params['noise_parameter_file']:
            self.noise_params = shelve.open(params['noise_parameter_file'], 'r')
        else:
            self.noise_params = None
        # If we are subtracting out foreground modes, open that file.
        if params['ts_foreground_mode_file']:
            self.foreground_modes = shelve.open(
                params['ts_foreground_mode_file'], 'r')
        else:
            self.foreground_modes = None
        # Set the map dimensioning parameters.
        self.n_ra = params['map_shape'][0]
        self.n_dec = params['map_shape'][1]
        spacing = params["pixel_spacing"]
        # Negative sign because RA increases from right to left.
        ra_spacing = -spacing/sp.cos(params['field_centre'][1]*sp.pi/180.)
        # To set some parameters, we have to read the first data file.
        first_file_name = (params["input_root"] + params["file_middles"][0]
                           + params["input_end"])
        Reader = fitsGBT.Reader(first_file_name, feedback=self.feedback)
        Blocks = Reader.read(0, (), force_tuple=True)
        self.n_chan = Blocks[0].dims[3]
        # Figure out which polarization indices to process.
        n_pols = Blocks[0].dims[1]
        pols = list(Blocks[0].field["CRVAL4"])
        self.pols = pols
        if not params['polarizations']:
            pol_inds = range(n_pols)
        else:
            pol_inds = []
            for stokes_par in params['polarizations']:
                for ii in range(n_pols):
                    if stokes_par == utils.polint2str(pols[ii]):
                        pol_inds.append(ii)
                        break
                else :
                    msg = "Polarization " + stokes_par + " not in Data files."
                    raise ce.ValueError(msg)
        # Get the centre frequncies of all the bands and figure out which bands
        # to process.
        n_bands = len(Reader.IF_set)
        band_inds = params["IFs"]
        if not band_inds:
            band_inds = range(n_bands)
        delta_freq = Blocks[0].field['CDELT1']
        self.delta_freq = delta_freq
        band_centres = []
        for Data in Blocks:
            Data.calc_freq()
            band_centres.append(Data.freq[self.n_chan//2])
        self.band_centres = band_centres
        del Blocks
        del Reader
        map_shape = (self.n_chan, self.n_ra, self.n_dec)
        for ii in pol_inds:
            for jj in band_inds:
                self.band_ind = jj
                band_centre = band_centres[jj]
                self.pol_ind = ii
                pol = pols[ii]
                if self.feedback > 1:
                    print ("Making map for polarization "
                           + utils.polint2str(pol) + " and band centred at "
                           +repr(round(band_centre/1e6)) + "MHz.")
                # Work out the file names.
                map_filename = (params["output_root"] + "dirty_map_"
                                + utils.polint2str(pol)
                                + "_" + str(int(round(band_centre/1e6)))
                                + '.npy')
                cov_filename = (params["output_root"] + "noise_inv_"
                                + utils.polint2str(pol)
                                + "_" + str(int(round((band_centre/1e6))))
                                + '.npy')
                # Initialization of the outputs.
                map = sp.zeros(map_shape, dtype=float)
                map = al.make_vect(map, axis_names=('freq', 'ra', 'dec'))
                map.set_axis_info('freq', band_centre, delta_freq)
                map.set_axis_info('ra', params['field_centre'][0], ra_spacing)
                map.set_axis_info('dec', params['field_centre'][1], 
                                   params['pixel_spacing'])
                if self.uncorrelated_channels:
                    cov_inv = al.open_memmap(cov_filename, mode='w+',
                                             shape=map_shape + map_shape[1:],
                                             dtype=float)
                    cov_inv = al.make_mat(cov_inv,
                                axis_names=('freq', 'ra', 'dec', 'ra', 'dec'),
                                row_axes=(0, 1, 2), col_axes=(0, 3, 4))
                else:
                    cov_inv = al.open_memmap(cov_filename, mode='w+',
                                             shape=map_shape + map_shape,
                                             dtype=float)
                    cov_inv = al.make_mat(cov_inv,
                                      axis_names=('freq', 'ra', 'dec',
                                                  'freq', 'ra', 'dec'),
                                      row_axes=(0, 1, 2), col_axes=(3, 4, 5))
                cov_inv.copy_axis_info(map)
                cov_inv[...] = 0
                self.map = map
                self.cov_inv = cov_inv
                # Do work.
                try:
                    self.make_map()
                except:
                    # I think yo need to do this to get a sensible error
                    # message, otherwise it will print cov_inv and fill the
                    # screen.
                    del self.cov_inv, cov_inv
                    raise
                # IO.
                # To write the noise_inverse, all we have to do is delete the
                # memeory map object.
                del self.cov_inv, cov_inv
                al.save(map_filename, map)
        if not self.noise_params is None:
            self.noise_params.close()

    def make_map(self):
        """Makes map for current polarization and band.
        
        This worker function has been split off for testing reasons.
        """
        
        params = self.params
        map = self.map
        cov_inv = self.cov_inv
        # This outer loop is over groups of files.  This is so we don't need to
        # hold the noise matrices for all the data in memory at once.
        file_middles = params['file_middles']
        n_files = len(file_middles)
        n_files_group = params['n_files_group']
        if n_files_group == 0:
            n_files_group = n_files
        for start_file_ind in range(0, n_files, n_files_group):
            this_file_middles = file_middles[start_file_ind
                                             :start_file_ind + n_files_group]
            if self.feedback > 1:
                print ("Processing data files %d to %d of %d."
                       % (start_file_ind, start_file_ind + n_files_group,
                          n_files))
            # Initialize lists for the data and the noise.
            data_list = []
            noise_list = []
            pointing_list = []
            # Loop an iterator that reads and preprocesses the data.
            # This loop can't be threaded without a lot of restructuring,
            # because iterate_data, preprocess_data and get_noise_parameter
            # all talk to each other through class attributes.
            for this_data in self.iterate_data(this_file_middles):
                # Unpack all the input data.
                time_stream, ra, dec, az, el, time, mask_inds = this_data
                n_time = time_stream.shape[1]
                if n_time < 5:
                    continue
                P = Pointing(("ra", "dec"), (ra, dec), map,
                             params['interpolation'])
                # Now build up our noise model for this piece of data.
                N = Noise(time_stream, time)
                N.add_mask(mask_inds)
                # The thermal part.
                thermal_noise = self.get_noise_parameter("thermal")
                try:
                    N.add_thermal(thermal_noise)
                except NoiseError:
                    continue
                # Frequency correlations.
                if params['frequency_correlations'] == 'mean':
                    mean_overf = self.get_noise_parameter("mean_over_f")
                    N.add_correlated_over_f(*mean_overf)
                elif params['frequency_correlations'] == 'measured':
                    # Correlated channel mode noise.
                    # The first frequency modes are completly deweighted. This
                    # was seen to be nessisary from looking at the noise
                    # spectra.
                    n_modes_discard = params['number_frequency_modes_discard']
                    try:
                        for ii in range(n_modes_discard):
                            mode_noise_params = self.get_noise_parameter(
                                    "over_f_mode_" + repr(ii))
                            # Instead of completly discarding, multiply 
                            # amplitude by 2.
                            #N.deweight_freq_mode(mode_noise_params[4])
                            mode_noise_params = ((mode_noise_params[0] * 2,)
                                                 + mode_noise_params[1:])
                            N.add_over_f_freq_mode(*mode_noise_params)
                        # The remaining frequency modes are noise weighed
                        # normally.
                        n_modes = params['number_frequency_modes']
                        for ii in range(n_modes_discard, n_modes):
                            mode_noise_params = self.get_noise_parameter(
                                "over_f_mode_" + repr(ii))
                            N.add_over_f_freq_mode(*mode_noise_params)
                        # All channel low frequency noise.
                        all_chan_noise_params = self.get_noise_parameter(
                                'all_channel')
                        # In some cases the corner frequncy will be so high that
                        # this scan will contain no information.
                        N.add_all_chan_low(*all_chan_noise_params)
                    except NoiseError:
                        if self.feedback > 1:
                            print ("Noise parameters risk numerical"
                                   " instability. Data skipped.")
                        continue
                elif params['frequency_correlations'] == 'None':
                    pass
                else:
                    raise ValueError("Invalid frequency correlations.")
                # Things to do along the time axis.  With no frequency
                # correlations, things are more stable.  Also, thermal noise
                # estimate will be high, so you need to deweight things extra.
                if params['deweight_time_mean']:
                    if params["frequency_correlations"] == 'None':
                        N.deweight_time_mean(T_large**2)
                    else:
                        N.deweight_time_mean(T_large**2)
                if params['deweight_time_slope']:
                    if params["frequency_correlations"] == 'None':
                        N.deweight_time_slope(T_large**2)
                    else:
                        N.deweight_time_slope(T_large**2)
                # If we've subtracted any modes out of the time stream,
                # deweight them.
                for mode_num in range(params['n_ts_foreground_modes']):
                    N.deweight_freq_mode(self.get_ts_foreground_mode(mode_num))
                # Store all these for later.
                pointing_list.append(P)
                noise_list.append(N)
                data_list.append(time_stream)
            # Threaded loop to finalize all the noises.  Each finalization
            # requires a large inverse.
            if self.feedback > 1:
                print "Finalizing %d noise blocks: " % len(data_list)
            # Initialize the queue of noise objects to work on.
            noise_queue = Queue()
            # Define a function that takes work off the queue and does it.
            def thread_work():
                while True:
                    thread_N = noise_queue.get()
                    # None will be the flag that there is no more work to do.
                    if thread_N is None:
                        return
                    else:
                        thread_N.finalize(frequency_correlations=(not
                                self.uncorrelated_channels),
                                preserve_matrices=False)
                        if self.feedback > 1:
                            sys.stdout.write('.')
                            sys.stdout.flush()

            # Start the worker threads.
            thread_list = []
            for ii in range(self.n_processes):
                T = threading.Thread(target=thread_work)
                T.start()
                thread_list.append(T)
            # Now put work on the queue for the threads to do.
            for N in noise_list:
                noise_queue.put(N)
            # At the end of the queue, tell the threads that they are done.
            for ii in range(self.n_processes):
                noise_queue.put(None)
            # Wait for the threads.
            for T in thread_list:
                T.join()
            if not noise_queue.empty():
                msg = "A thread had an error in Noise finialization."
                raise RuntimeError(msg)
            if self.feedback > 1:
                print
                print "Noise finalized."
            # Construct the dirty map.
            for ii in xrange(len(data_list)):
                time_stream = data_list[ii]
                N = noise_list[ii]
                P = pointing_list[ii]
                # Make the dirty map.
                weighted_time_stream = N.weight_time_stream(time_stream)
                map += P.apply_to_time_axis(weighted_time_stream)
            if self.feedback > 1:
                print "Dirty map done."
                print "Building map covariance. Frequecies finished:"
            # Now we have lists with all the data and thier noise.  
            # Accumulate it into a dirty map and its covariance.
            n_time_blocks = len(noise_list)
            # Initialize the queue of work to be done.
            index_queue = Queue()
            # I'm using a Lock to write to the memory map, not because I need
            # it but to hopefully keep simultaniouse writes from grinding
            # things to a halt.
            write_lock = threading.Lock()
            # Define a function that takes work off a queue and does it. 
            # Variable local to this function prefixed with 'thread_'.
            def thread_work():
                while True:
                    thread_inds = index_queue.get()
                    # None will be the flag that there is no more work to do.
                    if thread_inds is None:
                        return
                    if self.uncorrelated_channels:
                        thread_f_ind = thread_inds
                        thread_cov_inv_block = sp.zeros((self.n_ra, self.n_dec, 
                                    self.n_ra, self.n_dec), dtype=float)
                        for thread_kk in xrange(n_time_blocks):
                            thread_P = pointing_list[thread_kk]
                            thread_N = noise_list[thread_kk]
                            thread_P.noise_channel_to_map(thread_N, 
                                        thread_f_ind, thread_cov_inv_block)
                        cov_inv[thread_f_ind,...] += thread_cov_inv_block
                        if self.feedback > 1:
                            print thread_f_ind,
                            sys.stdout.flush()
                    else :
                        thread_f_ind = thread_inds[0]
                        thread_ra_ind = thread_inds[1]
                        thread_cov_inv_row = sp.zeros((self.n_dec, self.n_chan,
                                                       self.n_ra, self.n_dec),
                                                      dtype=float)
                        for thread_kk in xrange(n_time_blocks):
                            thread_P = pointing_list[thread_kk]
                            thread_N = noise_list[thread_kk]
                            thread_P.noise_to_map_domain(thread_N, 
                                    thread_f_ind, thread_ra_ind, 
                                    thread_cov_inv_row)
                        # Use a lock to try to stagger the processes and make
                        # them write at different times.
                        write_lock.acquire()
                        cov_inv[thread_f_ind,thread_ra_ind,...] += \
                                thread_cov_inv_row
                        write_lock.release()
                        if (self.feedback > 1
                            and thread_ra_ind == self.n_ra - 1):
                            print thread_f_ind,
                            sys.stdout.flush()

            # Start the worker threads.
            thread_list = []
            for ii in range(self.n_processes):
                T = threading.Thread(target=thread_work)
                T.start()
                thread_list.append(T)
            # Now put work on the queue for the threads to do.
            for ii in xrange(self.n_chan):
                if self.uncorrelated_channels:
                    index_queue.put(ii)
                else:
                    for jj in xrange(self.n_ra):
                        index_queue.put((ii, jj))
            # At the end of the queue, tell the threads that they are done.
            for ii in range(self.n_processes):
                index_queue.put(None)
            # Wait for the threads.
            for T in thread_list:
                T.join()
            if not index_queue.empty():
                msg = "A thread had an error while building map covariance."
                raise RuntimeError(msg)
            if self.feedback > 1:
                print
        # Now go through and make sure that the noise isn't singular by 
        # adding a bit of information to the diagonal.  This is equivalent to
        # setting a prior that all pixels are 0 += T_large. This does
        # bias the map maker slightly, but really shouldn't matter.
        if self.uncorrelated_channels:
            cov_view = cov_inv.view()
            cov_view.shape = (self.n_chan, (self.n_ra * self.n_dec)**2)
            cov_view[:,::self.n_ra * self.n_dec + 1] += 1.0 / T_large**2
        else:
            cov_inv.flat[::self.n_chan * self.n_ra * self.n_dec + 1] += \
                1.0 / T_large**2

#### Classes ####

class Pointing(object):
    """Class represents the pointing operator.

    The pointing operator converts from the map domain to the time domain in
    its native form and from the time domain to the map domain in its
    transposed form.
    
    Parameters
    ----------
    axis_names: tuple of strings
        The names of the axes in the map domain e.g. ("ra", "dec")
    coords : tuple of 1D arrays
        Tuple must be same length as `axis_names`.  The coordinates as a 
        function of time for each of the map axes.
    map : al.vect object
        The map that we will be gridding onto.  Map axes must include
        `axis_names`.  No modification to the map is made, only the axis
        information is used.
    scheme : string
        Gridding scheme to use.  Choices are 'nearest'.
    """

    def __init__(self, axis_names, coords, map, scheme='nearest'):
        # Sanity check some of the inputs.
        if len(axis_names) != len(coords):
            msg = "Got %d pointing axis names, but got %d coordinate arrays."
            raise ValueError(msg % (len(axis_names), len(coords)))
        n_pointings = coords[0].shape[0]
        for coordinate_array in coords:
            if coordinate_array.shape != (n_pointings,):
                msg = "Coordinate arrays must all be 1D and same length."
                raise ValueError(msg)
        # Save the data we need for the pointing.
        self._axis_names = axis_names
        self._coords = tuple(coords) # A tuple of 1D arrays.
        self._scheme = scheme
        # Get all coordinate information from the map.
        # Figure out which axes we are doing the pointing for, get the shape of
        # those axes, the coordinate centres and the pixel size.
        map_axis_indices = ()
        map_axes_shape = ()
        #map_axes_centre = ()
        #map_axes_delta = ()
        list_map_axis_names = list(map.axes)
        for axis_name in axis_names:
            axis_index = list_map_axis_names.index(axis_name)
            map_axis_indices += (axis_index,)
            map_axes_shape += (map.shape[axis_index],)
            #map_axes_centre += (map.info[axis_name + "_centre"],)
            #map_axes_delta += (map.info[axis_name + "_delta"],)
        self._map_axes = map_axis_indices
        self._map_shape = map_axes_shape
        
        # Store the full pointing matrix in sparse form.
        n_pointings = len(coords[0])
        n_coords = len(axis_names)
        self.dtype = np.float
        # Loop over the time stream and get the weights for each pointing.
        memory_allocated = False
        for ii in xrange(n_pointings):
            coordinate = ()
            for jj in xrange(n_coords):
                coordinate += (self._coords[jj][ii],)
            pixels, weights = map.slice_interpolate_weights(
                self._map_axes, coordinate, scheme)
            # On first iteration need to allocate memory for the sparse matrix
            # storage.
            if not memory_allocated:
                n_points_template = pixels.shape[0]
                self._pixel_inds = sp.zeros((n_pointings, n_coords,
                                              n_points_template), dtype=np.int)
                self._weights = sp.zeros((n_pointings, n_points_template),
                                         dtype=self.dtype)
                memory_allocated = True
            self._pixel_inds[ii,:,:] = pixels.transpose()
            self._weights[ii,:] = weights
    
    def get_sparse(self):
        """Return the arrays representing the pointing matrix in sparse form.

        Returns
        -------
        pixel_inds : array of ints
            Shape is (n_pointings, n_coordinates, n_pixels_per_pointing).
            These give which pixel entries are that are non zero for each
            pointing.  If a pointing is off the map, that row will be zeros.
        weights : array of floats
            Shape is (n_pointings, n_pixels_per_pointing).
            The entry of the pointing matrix corresponding with the matching
            pixel index.  If a pointing is off the map, that row will be zeros.
        """
        return self._pixel_inds, self._weights

    def apply_to_time_axis(self, time_stream, map_out=None):
        """Use this operator to convert a 'time' axis to a coordinate axis.
        
        This functions implements a fast matrix multiplication. It is roughly
        equivalent to using `algebra.partial_dot` except in axis placement.
        This function "replaces" the time axis with the map axes which is
        different from `partial_dot`'s behaviour.  This function is much more
        efficient than using `partial_dot`.
        
        For input
        `map`, the following operations should be equivalent, with the later
        much more efficient.
        
        # XXX This example is wrong and broken.
        >>> a = al.partial_dot(map)
        >>> b = self.apply_to_time_axis(map)
        >>> sp.allclose(a, b)
        True
        """
        
        if not isinstance(time_stream, al.vect):
            raise TypeError("Input data must be an algebra.vect object.")
        # Find the time axis for the input.
        for ii in range(time_stream.ndim):
            if time_stream.axes[ii] == 'time':
                time_axis = ii
                break
        else :
            raise ValueError("Input data vect doesn't have a time axis.")
        # Get some dimensions.
        n_pointings = self._pixel_inds.shape[0]
        n_pixels_template = self._pixel_inds.shape[2]
        n_axes = len(self._axis_names)
        if time_stream.shape[time_axis] != n_pointings:
            msg = ("Time stream data and pointing have different number of"
                   " time points.")
            raise ValueError(msg)
        # Get the shape and axis names of the output.
        out_shape = (time_stream.shape[:time_axis] + self._map_shape
                     + time_stream.shape[time_axis + 1:])
        out_axes = (time_stream.axes[:time_axis] + self._axis_names
                     + time_stream.axes[time_axis + 1:])
        # Allowcate output memory if not passed.
        if map_out is None:
            map_out = sp.zeros(out_shape, dtype=float)
            map_out = al.make_vect(map_out, axis_names=out_axes)
        else :
            if map_out.shape != out_shape:
                raise ValueError("Output array is the wrong shape.")
        # Initialize tuples that will index the input and the output.
        data_index = [slice(None),] * time_stream.ndim + [None]
        out_index = [slice(None),] * map_out.ndim
        # Loop over the time axis and do the dot.
        for ii in xrange(n_pointings):
            data_index[time_axis] = ii
            for kk in xrange(n_axes):
                out_index[time_axis + kk] = self._pixel_inds[ii,kk,:]
            map_out[tuple(out_index)] += (self._weights[ii,:]
                                   * time_stream[tuple(data_index)])
        return map_out

    def get_matrix(self):
        """Gets the matrix representation of the pointing operator."""

        n_pointings = self._pixel_inds.shape[0]
        n_coords = self._pixel_inds.shape[1]
        n_pixels_per_pointing = self._pixel_inds.shape[2]
        # Initialize the output matrix.
        matrix = sp.zeros((n_pointings,) + self._map_shape, dtype=self.dtype)
        matrix = al.make_mat(matrix, axis_names=("time",) + self._axis_names,
                             row_axes=(0,), col_axes=range(1, n_coords + 1))
        # Loop over the time stream and get the weights for each pointing.
        for ii in xrange(n_pointings):
            for jj in xrange( n_pixels_per_pointing):
                matrix[(ii,) + tuple(self._pixel_inds[ii,:,jj])] = \
                        self._weights[ii, jj]
        return matrix

    def noise_to_map_domain(self, Noise, f_ind, ra_ind, map_noise_inv):
        """Convert noise to map space.
        
        For performace and IO reasons this is done with a call to this function
        for each frequency row and each ra row.  All dec rows and all columns
        are handled in this function simultaniousely.

        This function is designed to be thread safe in that if it is called
        from two separate threads but with different `f_ind` or `ra_ind`, there
        should be no race conditions.
        """
        
        Noise._assert_finalized()
        if Noise._frequency_correlations:
            _mapmaker_c.update_map_noise_chan_ra_row(Noise.diagonal_inv,
                    Noise.freq_modes, Noise.time_modes, Noise.freq_mode_update,
                    Noise.time_mode_update, Noise.cross_update, 
                    self._pixel_inds, self._weights, f_ind, ra_ind,
                    map_noise_inv)
        else:
            msg = ("Noise object has no frequency correlations.  Use "
                   " `noise_channel_to_map` instead.")
            raise RuntimeError(msg)

    def noise_channel_to_map(self, Noise, f_ind, map_noise_inv):
        """Convert noise to map space.
        
        Use this function over `noise_to_map_domain` if the noise has no
        frequency correlations.
        """

        Noise._assert_finalized()
        if not Noise._frequency_correlations:
            _mapmaker_c.update_map_noise_independant_chan(Noise.diagonal_inv,
                    Noise.time_modes, Noise.time_mode_update, self._pixel_inds,
                    self._weights, f_ind, map_noise_inv)
        else:
            msg = ("Noise object has frequency correlations.  Use "
                   " `noise_to_map_domain` instead.")
            raise RuntimeError(msg)


class NoiseError(Exception):
    """Exception to raise if the there is something wrong with the noise and
    this peice of data should be ignored.
    """
    pass

class Noise(object):
    """Object that represents the noise matrix for time stream data.
    
    The noise matrix is represented as separate components each with different
    symetries.  This is so the huge matrix does not have to be stored.

    Parameters
    ----------
    time_strem_data : al.vect object
        The data for which we want to represent the noise matrix.  Only meta
        data is used, not the acctual data.
    time : 1D array
        The time axis of the data.
    """

    # Internal nomiclature: The noise matrix is divided into three parts:  The
    # 'diagonal' contains one weight for every data point and represents a
    # fully diagonal matrix.  It has contributions from thermal noise and
    # deweights masked points.
    # The 'frequency_modes' part is the only part that
    # couples frequencies.  There are only a small number of modes along the
    # frequency axes, but each mode has a full time-time covariance matrix.
    # The modes are assumed to be uncorrelated.  The number of frequency modes
    # is named 'm'
    # The 'time_modes' part deweights certain modes along the time axis.  They
    # are uncorrelated between frequencies and the same for each frequency.
    # The number time modes is named 'q'.
    # The 'update' term will generally refer to the second term in the binomial
    # inverse identity, and the 'update_modes' refers to rotation matrices in
    # the identity.

    # ---- Initialization methods. ----

    def __init__(self, time_stream_data, time):
        if len(time_stream_data.shape) != 2:
            raise ValueError("Only 2D data suported (freq, time).")
        self.n_chan = time_stream_data.shape[0]
        self.n_time = time_stream_data.shape[1]
        self.info = dict(time_stream_data.info)
        self._finalized = False
        # Some functions make this assumption.
        if not sp.alltrue(sp.diff(time) > 0):
            raise ValueError("Time must be monotonically increasing.")
        self.time = time

    def _assert_not_finalized(self):
        """Make sure the noise matrix is not finalized and can be modified."""
        if self._finalized:
            raise AssertionError("Noise model closed for modification.")
    
    def _assert_finalized(self):
        """Make sure the noise matrix is finalized and can not be modified."""
        if not self._finalized:
            raise AssertionError("Noise model still being modified.")

    def initialize_diagonal(self):
        """Create the diagonal part of the noise matrix if it doesn't exist."""
        
        self._assert_not_finalized()
        # TODO: Need to copy the axis info from self.info.
        if hasattr(self, "diagonal"):
            return
        diagonal = sp.zeros((self.n_chan, self.n_time), dtype=float)
        diagonal = al.make_mat(diagonal, axis_names=("freq", "time"), 
                               row_axes=(0, 1), col_axes=(0, 1))
        self.diagonal = diagonal

    def add_time_modes(self, n_new_modes=1):
        """Initialize time noise modes.
        """

        self._assert_not_finalized()
        # TODO: Need to copy the axis info from self.info.
        if hasattr(self, "time_modes"):
            current_q = self.time_modes.shape[0]
            new_q = current_q + n_new_modes
            old_time_modes = self.time_modes
            old_time_mode_noise = self.time_mode_noise
            time_modes = sp.zeros((new_q, self.n_time), 
                                    dtype=float)
            time_mode_noise = sp.zeros((new_q, self.n_chan, self.n_chan),
                                       dtype=float)
            time_modes[:current_q,:] = old_time_modes
            time_mode_noise[:current_q,:,:] = old_time_mode_noise
        else :
            current_q = 0
            time_modes = sp.zeros((n_new_modes, self.n_time), 
                                  dtype=float)
            time_mode_noise = sp.zeros((n_new_modes, self.n_chan, self.n_chan),
                                       dtype=float)
        time_modes = al.make_mat(time_modes, axis_names=("time_mode", "time"), 
                                 row_axes=(0,), col_axes=(1,))
        time_mode_noise = al.make_mat(time_mode_noise,
                                      axis_names=("time_mode", "freq", "freq"), 
                                      row_axes=(0,1), col_axes=(0,2))
        self.time_modes = time_modes
        self.time_mode_noise = time_mode_noise
        return current_q
    
    def add_freq_modes(self, n_new_modes=1):
        """Initialize frequency noise modes.
        """
        
        self._assert_not_finalized()
        # TODO: Need to copy the axis info from self.info.
        if hasattr(self, "freq_modes"):
            current_m = self.freq_modes.shape[0]
            new_m = current_m + n_new_modes
            old_freq_modes = self.freq_modes
            old_freq_mode_noise = self.freq_mode_noise
            freq_modes = sp.zeros((new_m, self.n_chan), 
                                    dtype=float)
            freq_mode_noise = sp.zeros((new_m, self.n_time, self.n_time),
                                       dtype=float)
            freq_modes[:current_m,:] = old_freq_modes
            freq_mode_noise[:current_m,:,:] = old_freq_mode_noise
        else :
            current_m = 0
            freq_modes = sp.zeros((n_new_modes, self.n_chan), 
                                  dtype=float)
            freq_mode_noise = sp.zeros((n_new_modes, self.n_time, self.n_time),
                                       dtype=float)
        freq_modes = al.make_mat(freq_modes, axis_names=("freq_mode", "freq"), 
                                 row_axes=(0,), col_axes=(1,))
        freq_mode_noise = al.make_mat(freq_mode_noise,
                                      axis_names=("freq_mode", "time", "time"), 
                                      row_axes=(0,1), col_axes=(0,2))
        self.freq_modes = freq_modes
        self.freq_mode_noise = freq_mode_noise
        return current_m

    # ---- Methods that build up the noise matrix. ----

    def add_thermal(self, thermal_levels):
        """Add a thermal component to the noise.
        
        This modifies the diagonal part of the noise matrix.

        Parameters
        ----------
        thermal_levels : 1D array
            The noise level of each channel (frequency).  This should be the
            thermal variance, in K**2, not K**2/Hz.
        """
        
        self.initialize_diagonal()
        if (not sp.all(sp.isfinite(thermal_levels))
            or sp.any(thermal_levels < T_small**2)):
            raise ValueError("Non finite thermal noise.")
        if isinstance(thermal_levels, sp.ndarray):
            self.diagonal += thermal_levels[:, None]
        else:
            self.diagonal += thermal_levels

    def add_mask(self, mask_inds):
        """Add a mask to the noise.
        
        This modifies the diagonal part of the noise matrix.

        Parameters
        ----------
        mask_inds : 2 element tuple of integer arrays.
            The locations of data points to be masked.
        """
        
        self.initialize_diagonal()
        self.diagonal[mask_inds] = T_infinity**2

    def deweight_time_mean(self, T=T_large**2):
        """Deweights time mean in each channel.

        This modifies the part of the noise matrix that is the same for each
        channel.
        """
        
        start = self.add_time_modes(1)
        self.time_modes[start,:] = 1.0 / sp.sqrt(self.n_time)
        # We deweight by T_medium, enough such that the mean doesn't interfere
        # with other modes in the map making, but not enough to actually
        # eliminate mean information (will have to be done explicitly) as doing
        # so would kill numerical stability.
        self.time_mode_noise[start,...] = (sp.eye(self.n_chan, dtype=float) 
                                           * T * self.n_time)

    def deweight_time_slope(self, T=T_large**2):
        """Deweights time slope in each channel.
        
        This modifies the part of the noise matrix that is the same for each
        channel.
        """

        start = self.add_time_modes(1)
        mode = self.time - sp.mean(self.time)
        mode *= 1.0 / sp.sqrt(sp.sum(mode**2)) 
        self.time_modes[start,:] = mode
        self.time_mode_noise[start,...] = (sp.eye(self.n_chan, dtype=float) 
                                           * T * self.n_time)

    def add_correlated_over_f(self, amp, index, f0):
        """Add 1/f noise that is perfectly correlated between frequencies.

        This modifies the the correlated mode part of the noise matrix. It adds
        a single mode with equal amplitude in all frequencies.

        Parameters
        ----------
        amp : float
            The amplitude of the noise at `f0`.
        index : float
            The spectral index of the spectrum (normaly near -1).
        f0 : float
            The pivot frequency.
        """
        
        mode = 1.0 / sp.sqrt(self.n_chan)
        thermal = 0.
        # Multiply by number of channels since this is the mean, not the sum
        # (Sept 15, 2011 in Kiyo's notes).
        amp *= self.n_chan
        self.add_over_f_freq_mode(amp, index, f0, thermal, mode)
        return
        
    def add_over_f_freq_mode(self, amp, index, f0, thermal, mode):
        """Add `1/f + const` noise to a given frequency mode."""
        
        # Too steep a spectra will crash finialize.
        time = self.time
        # Build the matrix.
        time_deltas = abs(time[:, None] - time)
        # Smallest time step.
        dt = sp.amin(abs(sp.diff(time)))
        # Time step for calculating the correlation function.
        # Over sample for precise interpolation.
        dt_calc = dt / 4
        n_lags = sp.amax(time_deltas) // dt_calc + 5
        # Calculate the correlation function at these lags.
        correlation_function = noise_power.calculate_overf_correlation(amp, 
            index, f0, dt_calc, n_lags)
        # If we are adding too much noise, we risk making the matrix singular.
        if sp.amax(correlation_function) > T_huge**2:
            print "Freq mode max:", sp.amax(correlation_function)
            raise NoiseError("Extremely high 1/f risks singular noise.")
        start_mode = self.add_freq_modes(1)
        self.freq_modes[start_mode,:] = mode
        corr_func_interpolator = \
            interpolate.interp1d(sp.arange(n_lags) * dt_calc,
                                 correlation_function, kind='linear')
        noise_mat = corr_func_interpolator(time_deltas)
        # Add the thermal part to the diagonal.
        # XXX This is acctually quite wrong.  The below is probably a bit
        # better.
        BW = 1. / 2. / dt
        noise_mat.flat[::self.n_time + 1] += thermal * BW * 2
        # XXX This isn't strictly being done correctly, but to do better we
        # need to measure the noise differently.
        #mode_norm = sp.sum(mode**2)  # Usually unity.
        #self.initialize_diagonal()
        #BW = 1. / 2. / dt
        #self.diagonal += mode[:,None]**2 * thermal * BW * 2 / mode_norm
        self.freq_mode_noise[start_mode,...] = noise_mat

    def deweight_freq_mode(self, mode, T=T_large**2):
        """Completly deweight a frequency mode."""
        
        n_chan = self.n_chan
        n_time = self.n_time
        start_mode = self.add_freq_modes(1)
        self.freq_modes[start_mode,:] = mode
        self.freq_mode_noise[start_mode,...]  = (sp.eye(n_time, dtype=float)
                                                 * T * n_chan)

    def add_all_chan_low(self, amps, index, f_0):
        """Deweight frequencies below, and a bit above 'f_0', 
        with 1/f like noise."""
        
        time = self.time
        n_chan = self.n_chan
        f_0 = float(f_0)
        # A bit of extra time range so the map maker doesn't think the modes
        # are periodic.
        extra_time_factor = 1.25
        time_extent = (time[-1] - time[0]) * extra_time_factor
        # Normalize the time to -pi <= t < pi
        time_normalized = time - time[0] - time_extent / 2.
        time_normalized *= 2 * sp.pi / time_extent
        df = 1. / time_extent
        # If there are no frequencies below f_0, do nothing.
        if f_0 * 1.5 <= df :
            return
        # Figure out how many modes to deweight.
        # I could concievably make the cut-off index dependant.
        frequencies = sp.arange(df, f_0 * 2.0, df)
        n_f = len(frequencies)
        n_modes = 2 * n_f
        if n_modes > self.n_time // 4:
            print f_0, n_modes, self.n_time
            raise NoiseError("To many time modes to deweight.")
        # Add allowcate mememory for the new modes.
        start_mode = self.add_time_modes(n_modes)
        # Loop through and fill the modes.
        for ii, f in enumerate(frequencies):
            this_amp_factor = (f / f_0)**index * df
            if sp.any(this_amp_factor > f_large * self.n_time):
                print "time mode amplitude:", this_amp_factor
                raise NoiseError("Deweighting amplitude too high.")
            this_amps = this_amp_factor * amps
            # The cosine (symetric) mode.
            cos_mode = sp.cos((ii+1) * time_normalized)
            norm_cos = sp.sqrt(sp.sum(cos_mode**2))
            cos_mode /= norm_cos
            cos_amp = this_amps * norm_cos**2
            self.time_modes[start_mode + 2*ii,:] = cos_mode
            cos_noise = sp.zeros((n_chan, n_chan), dtype=float)
            cos_noise.flat[::n_chan + 1] = cos_amp
            self.time_mode_noise[start_mode + 2*ii,:,:] = cos_noise
            # The sine (antisymetric) mode.
            sin_mode = sp.sin((ii+1) * time_normalized)
            norm_sin = sp.sqrt(sp.sum(sin_mode**2))
            sin_mode /= norm_sin
            sin_amp = this_amps * norm_sin**2
            self.time_modes[start_mode + 2*ii + 1,:] = sin_mode
            sin_noise = sp.zeros((n_chan, n_chan), dtype=float)
            sin_noise.flat[::n_chan + 1] = sin_amp
            self.time_mode_noise[start_mode + 2*ii + 1,:,:] = sin_noise
      
    def finalize(self, frequency_correlations=True, preserve_matrices=True):
        """Tell the class that you are done building the matrix.
        """
        
        self._assert_not_finalized()
        n_time = self.n_time
        n_chan = self.n_chan
        # This diagonal part must be set or the matrix will be singular.
        if not hasattr(self, 'diagonal'):
            raise RuntimeError("Diagonal noise component not set.")
        diagonal_inv = self.diagonal**-1
        self.diagonal_inv = al.as_alg_like(diagonal_inv, self.diagonal)
        if frequency_correlations:
            self._frequency_correlations = True
            if not hasattr(self, "freq_modes"):
                self.add_freq_modes(0)
            if not hasattr(self, "time_modes"):
                self.add_time_modes(0)
            # Calculate the inverses of all matricies.
            freq_mode_inv = al.empty_like(self.freq_mode_noise)
            for ii in xrange(self.freq_mode_noise.shape[0]):
                freq_mode_inv[ii,...] = \
                        linalg.inv(self.freq_mode_noise[ii,...])
                # Check that its positive definiate.
                A = self.freq_mode_noise[ii].view()
                A.shape = (n_time,) * 2
                e, v = linalg.eigh(A)
                e_max = max(e)
                if not sp.all(e > -1e-7*e_max) or e_max < 0:
                    print e
                    msg = ("Some freq_mode noise components not positive"
                           " definate.")
                    raise RuntimeError(msg)
                if hasattr(self, 'flag'):
                    print "Initial freq_mode condition number:", max(e)/min(e)
            time_mode_inv = al.empty_like(self.time_mode_noise)
            for ii in xrange(self.time_mode_noise.shape[0]):
                time_mode_inv[ii,...] = linalg.inv(
                    self.time_mode_noise[ii,...])
                A = self.time_mode_noise[ii].view()
                A.shape = (n_chan,) * 2
                e, v = linalg.eigh(A)
                e_max = max(e)
                if not sp.all(e > -1e-7*e_max) or e_max < 0:
                    print e
                    msg = ("Some time_mode noise components not positive"
                           " definate.")
                if hasattr(self, 'flag'):
                    print "Initial time_mode condition number:", max(e)/min(e)
            # The normal case when we are considering the full noise.
            # Calculate the term in the bracket in the matrix inversion lemma.
            # Get the size of the update term.
            # First, the rank of the correlated frequency part.
            m = self.freq_modes.shape[0]
            n_update =  m * n_time
            # Next, the rank of the all frequencies part.
            q = self.time_modes.shape[0]
            n_update += q * n_chan
            # Build the update matrix in blocks.
            freq_mode_update = sp.zeros((m, n_time, m, n_time), dtype=float)
            freq_mode_update = al.make_mat(freq_mode_update, 
                axis_names=('freq_mode', 'time', 'freq_mode', 'time'),
                row_axes=(0, 1), col_axes=(2, 3))
            cross_update = sp.zeros((m, n_time, q, n_chan), dtype=float)
            cross_update = al.make_mat(cross_update, 
                axis_names=('freq_mode', 'time', 'time_mode', 'freq'),
                row_axes=(0, 1), col_axes=(2, 3))
            time_mode_update = sp.zeros((q, n_chan, q, n_chan), dtype=float)
            time_mode_update = al.make_mat(time_mode_update, 
                axis_names=('time_mode', 'freq', 'time_mode', 'freq'),
                row_axes=(0, 1), col_axes=(2, 3))
            # Build the matrices.
            # Add the update mode noise in thier proper space.
            for ii in range(m):
                freq_mode_update[ii,:,ii,:] = freq_mode_inv[ii,:,:]
            for ii in range(q):
                time_mode_update[ii,:,ii,:] = time_mode_inv[ii,:,:]
            # Now transform the diagonal noise to this funny space and add
            # it to the update term. Do this one pair of modes at a time
            # to make things less complicated.
            for ii in xrange(m):
                for jj in xrange(m):
                    tmp_freq_update = sp.sum(self.freq_modes[ii,:,None]
                                             * self.freq_modes[jj,:,None]
                                             * diagonal_inv[:,:], 0)
                    freq_mode_update[ii,:,jj,:].flat[::n_time + 1] += \
                            tmp_freq_update
            for ii in xrange(m):
                for jj in xrange(q):
                    tmp_cross_update = (self.freq_modes[ii,None,:]
                                        * self.time_modes[jj,:,None]
                                        * diagonal_inv.transpose())
                    cross_update[ii,:,jj,:] += tmp_cross_update
            for ii in xrange(q):
                for jj in xrange(q):
                    tmp_time_update = sp.sum(self.time_modes[ii,None,:]
                                             * self.time_modes[jj,None,:]
                                             * diagonal_inv[:,:], 1)
                    time_mode_update[ii,:,jj,:].flat[::n_chan + 1] += \
                            tmp_time_update
            # Put all the update terms in one big matrix and invert it.
            update_matrix = sp.empty((n_update, n_update), dtype=float)
            # Top left.
            update_matrix[:m * n_time,:m * n_time].flat[...] = \
                freq_mode_update.flat
            # Bottom right.
            update_matrix[m * n_time:,m * n_time:].flat[...] = \
                time_mode_update.flat
            # Top right.
            update_matrix[:m * n_time,m * n_time:].flat[...] = \
                cross_update.flat
            # Bottom left.
            tmp_mat = sp.swapaxes(cross_update, 0, 2)
            tmp_mat = sp.swapaxes(tmp_mat, 1, 3)
            update_matrix[m * n_time:,:m * n_time].flat[...] = \
                tmp_mat.flat
            update_matrix_inv = linalg.inv(update_matrix)
            # TODO: Some sort of condition number check here?
            #if hasattr(self, 'flag'):
            #    e, v = linalg.eigh(update_matrix)
            #    print "Update term condition number:", max(e)/min(e)
            # Copy the update terms back to thier own matrices and store them.
            freq_mode_update.flat[...] = \
                    update_matrix_inv[:m * n_time,:m * n_time].flat
            self.freq_mode_update = freq_mode_update
            time_mode_update.flat[...] = \
                    update_matrix_inv[m * n_time:,m * n_time:].flat
            self.time_mode_update = time_mode_update
            cross_update.flat[...] = \
                    update_matrix_inv[:m * n_time,m * n_time:].flat
            self.cross_update = cross_update
            # Set flag so no more modifications to the matricies can occure.
            self._finalized = True
            # Check that the diagonal is positive to catch catastrophic
            # inversion faileurs.
            self.check_inv_pos_diagonal()

        else:
            # Ignore the channel correlations in the noise.  Ignore freq_modes.
            self._frequency_correlations = False
            # TODO: raise a warning?
            if hasattr(self, "freq_modes"):
                print "Warning, frequency mode noise ignored."
            # Calculate a time mode update term for each frequency.
            q = self.time_modes.shape[0]
            # First get the diagonal inverse of the time mode noise.
            time_mode_noise_diag = self.time_mode_noise.view()
            time_mode_noise_diag.shape = (q, self.n_chan**2)
            time_mode_noise_diag = time_mode_noise_diag[:,::self.n_chan + 1]
            time_mode_noise_inv = 1.0/time_mode_noise_diag
            # Allocate memory for the update term.
            time_mode_update = sp.zeros((self.n_chan, q, q), dtype=float)
            time_mode_update = al.make_mat(time_mode_update,
                        axis_names=('freq', 'time_mode', 'time_mode'),
                        row_axes=(0, 1), col_axes=(0,2))
            # Add in the time mode noise.
            for ii in range(q):
                for jj in range(self.n_chan):
                    time_mode_update[jj,ii,ii] += time_mode_noise_inv[ii,jj]
            # Transform the diagonal noise to time_mode space and add it in.
            for ii in range(self.n_chan):
                time_mode_update[ii,...] += sp.sum(self.diagonal_inv[ii,:] 
                                                   * self.time_modes[:,None,:] 
                                                   * self.time_modes[None,:,:],
                                                   -1)
                time_mode_update[ii,...] = linalg.inv(time_mode_update[ii,...])
            self.time_mode_update = time_mode_update
            # Set flag so no more modifications to the matricies can occure.
            self._finalized = True
        # If desired, delete the noise matrices to recover memory.  Doing
        # this means the Noise cannot be 'unfinalized'.
        if not preserve_matrices:
            if hasattr(self, "freq_mode_noise"):
                del self.freq_mode_noise
            del self.time_mode_noise

    def check_inv_pos_diagonal(self, thres=-1./T_huge**2):
        """Checks the diagonal elements of the inverse for positiveness.
        """
        
        noise_inv_diag = self.get_inverse_diagonal()
        if sp.any(noise_inv_diag < thres):
            print (sp.sum(noise_inv_diag < 0), noise_inv_diag.size)
            print (noise_inv_diag[noise_inv_diag < 0], sp.amax(noise_inv_diag))
            time_mod.sleep(300)
            raise NoiseError("Inverted noise has negitive entries on the "
                             "diagonal.")

    # ---- Methods for using the Noise Matrix. ----

    def get_inverse_diagonal(self):
        self._assert_finalized()
        return _mapmaker_c.get_noise_inv_diag(self.diagonal_inv,
                    self.freq_modes, self.time_modes, self.freq_mode_update,
                    self.time_mode_update, self.cross_update)

    def get_inverse(self):
        """Get the full noise inverse.

        This function is more for testing than accually being usefull (since in
        production we will only use part of the inverse at a time).
        """
        
        self._assert_finalized()
        n_chan = self.n_chan
        n_time = self.n_time
        if self._frequency_correlations:
            if hasattr(self, 'flag'):
                e = self.diagonal_inv.flat[:]
                print "Diagonal condition number:", max(e)/min(e)
            freq_modes = self.freq_modes
            time_modes = self.time_modes
            # Get the size of the update term.
            # First, the rank of the correlated frequency part.
            m = self.freq_modes.shape[0]
            n_update =  m * n_time
            # Next, the rank of the all frequencies part.
            q = self.time_modes.shape[0]
            n_update += q * n_chan
            # Allowcate memory.
            out = sp.zeros((n_chan, n_time, n_chan, n_time), dtype=float)
            out = al.make_mat(out, axis_names=('freq','time','freq','time'),
                              row_axes=(0, 1), col_axes=(2, 3))
            # Loop over the frequency indeces to reduce workspace memory and
            # for simplicity.
            for ii in xrange(n_chan):
                this_freq_modes1 = freq_modes.index_axis(1, ii)
                for jj in xrange(n_chan):
                    # Get only the matricies that apply to this slice.
                    this_freq_modes2 = freq_modes.index_axis(1, jj)
                    this_cross1 = self.cross_update.index_axis(3, ii)
                    this_cross2 = self.cross_update.index_axis(3, jj)
                    this_time_update = self.time_mode_update.index_axis(1, ii)
                    this_time_update = this_time_update.index_axis(2, jj)
                    # The freq_mode-freq_mode block of the update term.
                    tmp_mat = al.partial_dot(this_freq_modes1, 
                                             self.freq_mode_update)
                    tmp_mat = al.partial_dot(tmp_mat, this_freq_modes2)
                    out[ii,:,jj,:] -= tmp_mat
                    # The off diagonal blocks.
                    tmp_mat = al.partial_dot(this_cross1, time_modes)
                    tmp_mat = al.partial_dot(this_freq_modes2, tmp_mat)
                    out[ii,:,jj,:] -= tmp_mat.transpose()
                    tmp_mat = al.partial_dot(this_freq_modes1, this_cross2)
                    tmp_mat = al.partial_dot(tmp_mat, time_modes)
                    out[ii,:,jj,:] -= tmp_mat
                    # Finally the time_mode-time_mode part.
                    tmp_mat = al.partial_dot(time_modes.mat_transpose(), 
                                             this_time_update)
                    tmp_mat = al.partial_dot(tmp_mat, time_modes)
                    out[ii,:,jj,:] -= tmp_mat
                    # Multply one side by the diagonal.
                    out[ii,:,jj,:] *= self.diagonal_inv[ii,:,None]
            # Add the identity.
            out.flat[::n_chan * n_time + 1] += 1.0
            # Multiply by the thermal term.
            out[:,:,:,:] *= self.diagonal_inv[:,:]
        else:
            time_modes = self.time_modes
            time_mode_update = self.time_mode_update
            diagonal_inv = self.diagonal_inv
            q = self.time_modes.shape[0]
            # Allocate memory for the output.
            out = sp.zeros((n_chan, n_time, n_time), dtype=float)
            out = al.make_mat(out, axis_names=('freq','time','time'),
                              row_axes=(0, 1), col_axes=(0, 2))
            # First get the update term.
            tmp_update = al.partial_dot(time_modes.mat_transpose(),
                                        time_mode_update)
            out -= al.partial_dot(tmp_update, time_modes)
            # Multiply both sides by the diagonal.
            out *= diagonal_inv[:,:,None]
            out *= diagonal_inv[:,None,:]
            # Finally add the diagonal in.
            out.shape = (n_chan, n_time**2)
            out[:,::n_time + 1] += diagonal_inv
            out.shape = (n_chan, n_time, n_time)
        return out

    def weight_time_stream(self, data):
        """Noise weight a time stream data vector.
        """
        
        self._assert_finalized()
        time_modes = self.time_modes
        time_mode_update = self.time_mode_update
        diagonal_inv = self.diagonal_inv
        # Noise weight by the diagonal part of the noise.
        # These two lines replace an algebra library function.  They are much
        # faster for this case.
        diag_weighted = diagonal_inv * data
        diag_weighted = al.as_alg_like(diag_weighted, data)
        # Calculate the update term carrying the freq modes and the time modes
        # through separately.
        # Transform to the update space.
        tmp_update_term_time = al.partial_dot(time_modes, diag_weighted)
        # Multiply by the update matrix.
        update_term_time = al.partial_dot(time_mode_update,
                                          tmp_update_term_time)
        if self._frequency_correlations:
            # Calculate terms that couple frequencies.
            freq_modes = self.freq_modes
            freq_mode_update = self.freq_mode_update
            cross_update = self.cross_update
            # Transform to the update space.
            tmp_update_term_freq = al.partial_dot(freq_modes, diag_weighted)
            # Multiply by the update matrix.
            update_term_freq = (al.partial_dot(freq_mode_update, 
                                               tmp_update_term_freq)
                                + al.partial_dot(cross_update,
                                                 tmp_update_term_time))
            update_term_time += al.partial_dot(cross_update.mat_transpose(),
                                               tmp_update_term_freq)
            # Transform back.
            update_term_freq = al.partial_dot(freq_modes.mat_transpose(),
                                              update_term_freq)
            update_term_time = al.partial_dot(time_modes.mat_transpose(),
                                              update_term_time)
            # Combine.
            update_term = update_term_freq + update_term_time.transpose()
        else:
            # Transform back.
            update_term_time = al.partial_dot(time_modes.mat_transpose(),
                                              update_term_time)
            update_term = update_term_time
        # Final multiply by the diagonal component.
        update_term = al.partial_dot(diagonal_inv, update_term)
        # Combine the terms.
        out = diag_weighted - update_term
        return out 


#### Utilities ####

def trim_time_stream(data, coords, lower_bounds, upper_bounds):
    """Discards data that is out of the map bounds.

    Parameters
    ----------
    data : algebra.vect
        The time stream data.  The vect object axes must be ('freq', 'time').
    coords : tuple (length `n`) of arrays
        The pointing information for `n` map coordinates as a function of time.
    lower_bounds : tuple (length `n`) of floats.
        Lower map boundaries for each coordinate.
    upper_bounds : tuple (length `n`) of floats.
        Upper map boundaries for each coordinate.

    Returns
    -------
    trimmed_data : algebra.vect
        New data vector with ashortened time axes to exclude off map
        pointings.
    inds : 1d array of ints
        Which tim stream array indices were retained (this array can be used 
        index other arrays).
    """
    
    n_time = data.shape[1]
    # Loop through the pointing information and see which times are on the map.
    inds = []
    for ii in range(n_time):
        for jj in range(len(coords)):
            if (coords[jj][ii] <= lower_bounds[jj]
                or coords[jj][ii] >= upper_bounds[jj]):
                break
        else:
            inds.append(ii)
    inds = sp.asarray(inds, dtype=int)
    # Shorten the time axis and cast as a vect object.
    trimmed_data = sp.ascontiguousarray(data[:,inds])
    trimmed_data = al.make_vect(trimmed_data, axis_names=('freq', 'time'))
    return trimmed_data, inds
