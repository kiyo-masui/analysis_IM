#!/usr/bin/python
"""Convert psrfits with subintegration extension to SDfits.

This script does a rough conversion from psrfits to SDfits.  It was written
with the HI intensity mapping analysis pipeline in mind, so there is no
guarantee that it will conform Exactly to the SDfits standard.  The idea is
simply to arange the data in a familiar way for SDfits.

This is run using a parameter file using kiyopy parameter reading
(github.com/kiyo-masui/kiyopy).

This also needs to go into the antenna fits files to get the precise pointing
information.
"""

import datetime
import multiprocessing as mp
import sys
import os
import glob
import subprocess
import warnings

import numpy.ma as ma
import scipy as sp
import scipy.interpolate as interp
import scipy.fftpack as fft
import pyfits
import matplotlib.pyplot as plt

from time_stream import rotate_pol, cal_scale, rebin_freq
from core import fitsGBT, data_block
from kiyopy import parse_ini, utils
import kiyopy.custom_exceptions as ce
import kiyopy.pickle_method

params_init = {# Inputs.
               # Everything up to the scan number.
               "guppi_input_root" : "./",
               # Everything after the scan number
               "guppi_input_end" : ".fits",
               # The directory containing the fits scan log file (with antenna
               # and go fits subdirectories).
               "fits_log_dir" : "./",
               # Where to write converted files.
               "output_root" : "./",
               # Scans to convert.  List of scan numbers to convert.  Other
               # scans in same proceedure are automatically selected if
               # combine_map_scans = True
               "scans" : (0,),
               # Whethar to combine all the scans from a proceedure into a
               # single file.
               "combine_map_scans" : False,
               # Split cal on and cal off.
               "partition_cal" : False,
               # How many time bins to average over when resampling.  Must be a
               # power of 2 and less than 2048.
               "time_bins_to_average" : 1
               }
prefix = ''

class Converter(object) :

    def __init__(self, input_file_or_dict=None, feedback=2) :    
        # Read the parameter file.
        self.params = parse_ini.parse(input_file_or_dict, params_init, 
                                      prefix=prefix)
    
    def execute(self, nprocesses=1) :
        params = self.params
        scans = list(params["scans"])
        # Make sure that the output directory exists.
        utils.mkparents(params["output_root"])

        # Now we need to read in the Scan Log fits file.
        log_dir = params["fits_log_dir"]
        scan_log_list = pyfits.open(log_dir + "/ScanLog.fits", "readonly")
        # From the header we need the project session.
        session = scan_log_list[0].header["PROJID"].split('_')[-1]
        scan_log = scan_log_list[1].data
        self.scan_log = scan_log
        
        # Keep track of scans already processed because some scans are 
        # processed by being in the same map as another.
        finished_scans = []
        for initial_scan in scans :
            if initial_scan in finished_scans :
                continue
            self.initial_scan = initial_scan

            # Open the go fits file.
            scan_log_files = scan_log.field('FILEPATH')[
                scan_log.field('SCAN')==initial_scan]
            # Find the one that's a GO fits file.
            found_file = False
            for file_name in scan_log_files :
                if 'GO' in file_name :
                    go_file = (log_dir + "/GO" + file_name.split("GO")[1])
                    found_file = True
                    break
            if not found_file :
                raise ce.DataError("GO file not in Scan Log")
            go_hdu = pyfits.open(go_file)[0].header

            # From the go information get the source and the scan type.
            object = go_hdu["OBJECT"].strip()
            self.proceedure = go_hdu["PROCNAME"].strip().lower()

            if params["combine_map_scans"] :
                # Read the go file and figure out all the scans in the same map.
                # Check the go files for all scans make sure everything is
                # consistant.
                self.n_scans_proc = go_hdu["PROCSIZE"]
                # Which scan this is of the sequence (1 indexed).
                self.initial_scan_ind = go_hdu["PROCSEQN"]
                scans_this_file = (sp.arange(self.n_scans_proc, dtype=int) + 1 -
                                   self.initial_scan_ind + initial_scan)
                scans_this_file = list(scans_this_file)

            else :
                scans_this_file = [initial_scan]
            finished_scans += scans_this_file

            # Initialize a list to store all the data that will be saved to a
            # single fits file (generally 8 scans).
            Block_list = []
            # Loop over the scans to process for this output file.
            for scan in scans_this_file :
                # The acctual reading of the guppi fits file needs to be split
                # off in a different process due to a memory leak in pyfits.
                
                # Make a pipe over which we will receive out data back.
                P_here, P_far = mp.Pipe()
                # Start the forked process.
                p = mp.Process(target=self.processfile, args=(scan, P_far))
                p.start()
                Data = P_here.recv()
                p.join()
                
                # Store our processed data.
                if Data == -1 :
                    # Scan proceedure aborted.
                    raise ce.DataError("Scan proceedures do not agree."
                            " Perhase a scan was aborted. Scans: " + str(scan)
                            + ", " + str(initial_scan) + " in directory: "
                            + log_dir)
                elif Data is None :
                    warnings.warn("Missing psrfits file: " + guppi_file)
                else :
                    Block_list.append(Data)

            # End loop over scans (input files).
            # Now we can write our list of scans to disk.
            if len(scans_this_file) > 1 :
                str_scan_range = (str(scans_this_file[0]) + '-' +
                                  str(scans_this_file[-1]))
            else :
                str_scan_range = str(scans_this_file[0])
            out_file = (params["output_root"] + session + '_' + object + '_' +
                        self.proceedure + '_' + str_scan_range + '.fits')
            
            # Output data is pretty large so we'd better protect the pyfits part
            # in a process lest memory leaks kill us.
            if len(Block_list) > 0 :
                p = mp.Process(target=out_write, args=(Block_list, out_file))
                p.start()
                del Block_list
                p.join()
        # End loop over maps (output files or input 'scans' parameter).
            
    def processfile(self, scan, Pipe) :
        params = self.params
        scan_log = self.scan_log
        log_dir = params["fits_log_dir"]

        # Whethar we will fold the data to look for the cal or not.
        get_cal = params["partition_cal"]
        if get_cal :
            ncal = 2
        else :
            ncal = 1

        # First find the GO fits file.
        scan_log_files = scan_log.field('FILEPATH')[
            scan_log.field('SCAN')==scan]
        # Find the one that's a GO fits file.
        found_file = False
        for file_name in scan_log_files :
            if 'GO' in file_name :
                go_file = (log_dir + "/GO" + file_name.split("GO")[1])
                found_file = True
                break
        if not found_file :
            raise ce.DataError("GO file not in Scan Log")
        # If we are combining multiple scans from one proceedure, make sure
        # that they all came from the same one (ie that the proceedure
        # wasn't aborted).
        if params["combine_map_scans"] :
            go_hdu = pyfits.open(go_file)[0].header
            if (self.n_scans_proc != go_hdu["PROCSIZE"] or go_hdu["PROCSEQN"] 
                - self.initial_scan_ind != scan - self.initial_scan) :
                # Rather than crashing here, send an error code and crash the
                # main program.
                Pipe.send(-1)
                return

        # Now find the Antenna fits file.
        found_file = False
        for file_name in scan_log_files :
            if 'Antenna' in file_name :
                antenna_file = (log_dir + "/Antenna" +
                                file_name.split("Antenna")[1])
                found_file = True
                break
        if not found_file :
            raise ce.DataError("Antenna file not in Scan Log")
        # Open the antenna fits file.
        antenna_hdu_list = pyfits.open(antenna_file, 'readonly')
        ant_main = antenna_hdu_list[0].header
        ant_header = antenna_hdu_list[2].header
        ant_data = antenna_hdu_list[2].data
        # Get the time of the scan midpoint in seconds since 00:00 UTC.
        ant_scan_mid = ant_main["SDMJD"]%1.0 * 24.0 * 3600.0
        n_ant = ant_header["NAXIS2"]
        # Antenna sample speed is 10 Hz always.
        ant_periode = 0.100
        ant_times = ant_scan_mid + ant_periode*(sp.arange(n_ant) -
                                                n_ant/2.0)
        # Get the az and el out of the antenna.
        ant_az = ant_data.field("OBSC_AZ")
        ant_el = ant_data.field("OBSC_EL")
        az_interp = interp.interp1d(ant_times, ant_az,
                                          kind="nearest")
        el_interp = interp.interp1d(ant_times, ant_el,
                                          kind="nearest")
        
        # Get the guppi file name.
        guppi_file = (params["guppi_input_root"] + "%04d"%scan +
                      params["guppi_input_end"])
        # Sometimes guppi files are just missing.  This shouldn't crash the
        # program.
        if not os.path.isfile(guppi_file) :
            Pipe.send(None)
            return
        print "Converting file: " + guppi_file,
        if params['combine_map_scans'] :
           print (" (scan " + str(go_hdu["PROCSEQN"]) + " of " +
                  str(self.n_scans_proc) + ')')
        else :
            print
        psrhdu_list = pyfits.open(guppi_file, 'readonly')
        # A record array with each row holding about 2 seconds of data.
        psrdata = psrhdu_list[1].data
        psrheader = psrhdu_list[1].header
        psrmain = psrhdu_list[0].header
        # Some numbers we need to pull from the header.
        if not psrheader["TTYPE17"] == "DATA" :
            raise ce.DataError("Expected to find DATA as the 17th field")
        data_shape = eval(psrheader["TDIM17"])
        if data_shape[0] != 1:
            raise ce.DataError("Data not shaped as expected.")
        # Change data shape to C ordering.
        data_shape = (list(data_shape[1:]))
        data_shape.reverse()
        data_shape = tuple(data_shape)
        # Number of times in a record (only a fraction of the total
        # number of times in the fits file).
        ntime = data_shape[0]
        npol = data_shape[1]
        nfreq = data_shape[2]
        
        # Number of times to average to together when rebinning.
        n_bins_ave = params["time_bins_to_average"]
        if ntime%n_bins_ave != 0 :
            raise ValueError("Number of time bins to average must divide "
                             "the number of time bins in a sub integration.")
        # Number of time bins AFTER rebinning
        n_bins = ntime//n_bins_ave
        
        # Figure out the file start time
        start_string = psrmain["DATE-OBS"] # UTC string including date.
        # Since 00h00 UTC.
        start_seconds = psrmain["STT_SMJD"] + psrmain["STT_OFFS"]
        # String that will be converted to the time at each sample.
        time_string = start_string.split('T')[0] + 'T'
        time_string = time_string + '%02d:%02d:%06.3f'
        # Make a time string for the next day in the off chance that the scan
        # strattles 00:00:00 UTC. This calculation is very rarely needed and
        # suficiently complex that I won't write the whole algorithm more
        # generally, lest I make a mistake.
        time_obj = datetime.datetime.strptime(start_string.split('T')[0],
                                              "%Y-%M-%d")
        time_obj = time_obj + datetime.timedelta(days=1)
        time_string_next_day = (time_obj.strftime("%Y-%M-%d") +
                                'T%02d:%02d:%06.3f')
        # Figure out the sample period after rebinning.
        sample_time = psrheader["TBIN"]
        resampled_time = sample_time*n_bins_ave
        
        # In current scan startegy, Zenith angle is approximatly 
        # constant over a file.  We will verify that this matches the antenna 
        # fits file as a good (but scan strategy specific) check.
        if self.proceedure == 'ralongmap' :
            zenith_angle = psrdata[0]["TEL_ZEN"]
            if not sp.allclose(90.0 - zenith_angle, ant_el, atol=0.1) :
                raise ce.DataError("Antenna el disagrees with guppi el.")

        # Allowcate memory for this whole scan.
        # final_shape is ordered (ntime, npol, ncal, nfreq) where ntime
        # is after rebining and combing the records.
        final_shape = (n_bins*len(psrdata), npol, ncal,
                              nfreq)
        Data = data_block.DataBlock(sp.empty(final_shape, dtype=float))
        # Allowcate memory for the the fields that are time dependant.
        Data.set_field('CRVAL2', sp.empty(final_shape[0]),
                       ('time',), '1D')
        Data.set_field('CRVAL3', sp.empty(final_shape[0]),
                       ('time',), '1D')
        Data.set_field('DATE-OBS', sp.empty(final_shape[0], dtype='S22'),
                       ('time',), '22A')

        # Copy as much field information as we can without looping.
        crpix1 = nfreq//2 + 1
        crval1 = psrdata.field("DAT_FREQ")[0, crpix1 - 1]*1e6
        Data.set_field('CRVAL1', crval1, (), '1D')
        Data.set_field('CRPIX1', crpix1, (), '1I')
        Data.set_field('SCAN', scan, (), '1I')
        Data.set_field('BANDWID', psrmain['OBSBW']*1e6, (), '1D')
        Data.set_field('CDELT1', psrheader['CHAN_BW']*1e6, (), '1D')
        Data.set_field('OBJECT', psrmain['SRC_NAME'], (), '32A')
        Data.set_field('EXPOSURE', resampled_time/ncal, (), '1D')
        if psrheader["POL_TYPE"].strip() == 'IQUV' :
            Data.set_field('CRVAL4', [1,2,3,4], ('pol',), '1I')
        else :
            raise ce.DataError("Unrecognized polarizations.")
        Data.set_field('OBSERVER', psrmain['OBSERVER'], (), '32A')
        if get_cal :
            Data.set_field('CAL', ['T', 'F'], ('cal',), '1A')
        else :
            Data.set_field('CAL', ['F'], ('cal',), '1A')
        Data.verify()

        # Loop over the records and modify the data.
        print "Record: ",
        for jj, record in enumerate(psrdata) :
            # Print the record we are currently on and force a flush to
            # standard out (so we can watch the progress update).
            print jj,
            sys.stdout.flush()
            # Get the data field.
            data = record["DATA"]
            # Convert the Data to the proper shape and type.
            data = format_data(data, ntime, npol, nfreq)
            # Now get the scale and offsets for the data and apply them.
            scls = sp.array(record["DAT_SCL"], dtype=sp.float32)
            scls.shape = (1, npol, nfreq)
            offs = sp.array(record["DAT_OFFS"], dtype=sp.float32)
            offs.shape = (1, npol, nfreq)

            if get_cal :
                # Figure out the number of time bins in a cal period.
                try :
                    cal_period = 1.0/psrmain["CAL_FREQ"]
                except ZeroDivisionError :
                    # Default value that we usually use.
                    cal_period = 0.065536
                n_bins_cal = cal_period/sample_time
                if n_bins_cal%1.0 > 0.00001 :
                    raise NotImplementedError("Need an integer number of "
                                              "samples per cal periode.")
                n_bins_cal = int(round(n_bins_cal))
                if n_bins_ave%n_bins_cal != 0 :
                    raise ValueError("You should rebin on a multiple of "
                                     "the cal period.  Change parameter "
                                     "time_bins_to_average.")
                # We can speed up the separating of the cal by rebinning in
                # time first.
                tmp_n_bins_ave = n_bins_ave//n_bins_cal
                data.shape = (n_bins, tmp_n_bins_ave, n_bins_cal, npol,
                              nfreq)
                # Use mean not median due to discritization.
                data = scls*sp.mean(data, 1) + offs
                data.shape = (n_bins*n_bins_cal, npol, nfreq)
                # Finally devided into cal on and cal off.
                data = separate_cal(data, n_bins_cal)
            else :
                # Just rebin in time and add a length 1 cal index.
                data.shape = (n_bins, n_bins_ave, npol, nfreq)
                data = sp.mean(data, 1)
                data = scls*data + offs
                data.shape = (n_bins, npol, ncal, nfreq)
            
            # Now figure out the timing information.
            time = start_seconds + record["OFFS_SUB"]
            time = time + (sp.arange(-n_bins/2.0 + 0.5, n_bins/2.0 + 0.5)
                           * resampled_time)
            # Make sure that our pointing data covers the guppi data time
            # range within the 0.1 second accuracy.
            if ((max(time) - max(ant_times) > 0.1) or 
                (min(time) - min(ant_times) < 0.1)) :
                raise ce.DataError("Guppi data time range not covered by "
                                   "antenna.")
            az = az_interp(time)
            el = el_interp(time)
            
            # Convert the time to hours, minutes and seconds.
            hours = time//3600
            minutes = (time%3600)//60
            seconds = time%60

            # Put all the data into the DataBlock for writing out.
            Data.data[jj*n_bins:(jj+1)*n_bins,:,:,:] = data
            Data.field['CRVAL2'][jj*n_bins:(jj+1)*n_bins] = az
            Data.field['CRVAL3'][jj*n_bins:(jj+1)*n_bins] = el
            # If the scan strattles midnight UTC, this is extra
            # complecated.
            for mm, kk in enumerate(xrange(jj*n_bins, (jj+1)*n_bins)) :
                if hours[mm] >= 24 :
                    tmp_time = (time_string_next_day
                                % (hours[mm]-24, minutes[mm], seconds[mm]))
                else :
                    tmp_time = time_string%(hours[mm], minutes[mm],
                                            seconds[mm])
                Data.field["DATE-OBS"][kk] = tmp_time
        Data.add_history('Converted from a guppi sub-integration fits file.',
                         (utils.abbreviate_file_path(guppi_file),
                          utils.abbreviate_file_path(antenna_file),
                          utils.abbreviate_file_path(go_file)))
        # End loop over records.
        print
        # Send Data back on the pipe.
        Pipe.send(Data)

# Functions for the most error prone parts of the above script.  Unit testing
# applied to these.

def out_write(Block_list, file_name) :
    """Fits file writer protected in a function so it can be multiprocessed, 
    thus avoiding memory leaks in pyfits."""

    Writer = fitsGBT.Writer(Block_list)
    Writer.write(file_name)

def format_data(data, ntime, npol, nfreq) :
    """Does the first reformating of the guppi data.  Reshapes it and converts
    the data type to float.
    
    The input data matrix is not preserved in this function.
    """
    if (not data.dtype == sp.uint8) or (len(data.shape) != 1) :
        raise TypeError("Expected flat uint8 data.")
    out_data = sp.empty((ntime, npol, nfreq), dtype=sp.float32)
    # Reshape the data making the last index the one that is
    # averaged over when we rebin.
    data.shape = (ntime, npol, nfreq)
    # Q, U, V are signed characters, not signed ones.  We need to
    # convert them.
    out_data[:,0,:] = data[:,0,:]
    data.dtype = sp.int8
    out_data[:,1:,:] = data[:,1:,:]
    data.dtype = sp.uint8

    return out_data

def get_cal_mask(data, n_bins_cal) :
    """Figures out where the cal is turning on and off.
    
    This function has been separated from separate_cal for testing purposes."""
    
    ntime = data.shape[0]
    
    if ntime%n_bins_cal != 0 :
        raise ValueError("Number of bins in a cal period must devide the number"
                        " of time bins in a subintegration (generally number"
                        " of bins in a cal period should be a power of 2).")

    # Fold the Stokes I data on the cal period to figure out the phase.
    # Sum over the frequencies.
    folded_data = sp.sum(data[:,0,:], -1, dtype=float)
    # Fold the data at the cal period.
    folded_data.shape = ((ntime//n_bins_cal, n_bins_cal) +
                         folded_data.shape[1:])
    folded_data = sp.mean(folded_data, 0)
    # Split the data into cal on and cal offs.
    base = sp.mean(folded_data)
    diff = sp.std(folded_data)
    transition = sp.logical_and(folded_data < base + 0.90*diff,
                                folded_data > base - 0.90*diff)
    if sp.sum(transition) > 2:
        profile_str = repr((folded_data - sp.mean(folded_data))
                           / sp.std(folded_data))
        raise ce.DataError("Cal profile has too many "
                           "transitions.  Profile array: "
                           + profile_str)
    # Find the last bin off before on.
    last_off_ind = None
    for kk in xrange(n_bins_cal) :
        if ((folded_data[kk] < base-0.80*diff) and not
            (folded_data[(kk+1)%n_bins_cal] < base-0.80*diff)):
            last_off_ind = kk
            break
    if last_off_ind is None :
        profile_str = repr((folded_data - sp.mean(folded_data))
                           / sp.std(folded_data))
        raise ce.DataError("Cal profile not lining up "
                           "correctly.  Profile array: "
                           + profile_str)
    # If last_off_ind+1 is a transitional bin, then we only need
    # to cut 2 bins.  If it's in the cal-on category, we should
    # cut 4 bins.
    # Lims are indicies (first index included, second
    # excluded).
    first_on = (last_off_ind + 2)%n_bins_cal
    if folded_data[(last_off_ind+1)%n_bins_cal] < base+0.80*diff :
        n_cal_state = n_bins_cal//2 - 1
        n_blank = 1
    else :
        n_cal_state = n_bins_cal//2 - 2
        n_blank = 2

    long_profile = sp.concatenate((folded_data, folded_data), 0)
    first_off = first_on + n_bins_cal//2
    if (sp.any(long_profile[first_on:first_on+n_cal_state] < base + 0.90*diff)
        or sp.any(long_profile[first_off:first_off+n_cal_state] > base -
                  0.90*diff)) :
        profile_str = repr((folded_data - sp.mean(folded_data))
                           / sp.std(folded_data))
        raise ce.DataError("Cal profile not lining up "
                           "correctly.  Profile array: "
                           + profile_str)

    return first_on, n_blank

def separate_cal(data, n_bins_cal) :
    """Function separates data into cal_on and cal off.
    
    No Guarantee that data argument remains unchanged."""
    
    # Allowcate memeory for output    
    ntime, npol, nfreq = data.shape
    n_bins_after_cal = ntime//n_bins_cal
    out_data = sp.empty((n_bins_after_cal, npol, 2, nfreq), dtype=sp.float32)
    
    # Get the phase offset of the cal.
    try :
        first_on, n_blank = get_cal_mask(data, n_bins_cal)
    except ce.DataError :
        print ": Discarded record due to bad profile. ",
        out_data[:] = float('nan')
    else :
        # How many samples for each cal state.
        n_cal_state = n_bins_cal//2 - n_blank
        first_off = (first_on + n_bins_cal//2) % n_bins_cal

        # Reshape data to add an index to average over.
        data.shape = (n_bins_after_cal, n_bins_cal) + data.shape[1:]

        # Get the masks for the on and off data.
        inds = sp.arange(n_bins_cal)
        if first_on == min((sp.arange(n_cal_state) +
                        first_on)% n_bins_cal) :
            on_mask = sp.logical_and(inds >= first_on, inds < 
                                     first_on+n_cal_state)
        else :
            on_mask = sp.logical_or(inds >= first_on, inds < 
                                (first_on + n_cal_state) % n_bins_cal)
        if first_off == min((sp.arange(n_cal_state) +
                        first_off)% n_bins_cal) :
            off_mask = sp.logical_and(inds >= first_off, inds < 
                                  first_off + n_cal_state)
        else :
            off_mask = sp.logical_or(inds >= first_off, inds < 
                                 (first_off + n_cal_state) % n_bins_cal)

        # Find cal on and cal off averages.  Always use mean not median due to
        # discretization noise.
        out_data[:,:,0,:] = sp.mean(data[:,on_mask,:,:], 1)
        out_data[:,:,1,:] = sp.mean(data[:,off_mask,:,:], 1)

    return out_data

params_init_checker = {
                      "input_root" : "",
                      "file_middles" : ("",),
                      "output_root" : "",
                      "output_end" : ".pdf"
                      }

class DataChecker(object) :

    def __init__(self, parameter_file_or_dict, feedback=2) :    
        self.params = parse_ini.parse(parameter_file_or_dict, 
                params_init_checker, prefix="")

    def execute(self) :
        params = self.params
        # Loop over files, make one set of plots per file.
        for middle in params["file_middles"] :
            p = mp.Process(target=self.processfile, args=(middle,))
            p.start()
            p.join()

    def process_file(self, middle) :
        """Split off to fix pyfits memory leak."""
        params = self.params
        # Construct the file name and read in all scans.
        file_name = params["input_root"] + middle + ".fits"
        Reader = fitsGBT.Reader(file_name)
        Blocks = Reader.read((), ())
        # Initialize a few variables.
        counts = 0
        cal_sum_unscaled = 0
        cal_sum = 0
        cal_time = ma.zeros((0, 4))
        sys_time = ma.zeros((0, 4))
        cal_noise_spec = 0
        # Get the number of times in the first block and shorten to a
        # number that should be smaller than all blocks.
        nt = int(Blocks[0].dims[0]*.9)
        # Get the frequency axis.  Must be before loop because the data is
        # rebined in the loop.
        Blocks[0].calc_freq()
        f = Blocks[0].freq
        for Data in Blocks :
            # Rotate to XX, YY etc.
            rotate_pol.rotate(Data, (-5, -7, -8, -6))
            this_count = ma.count(Data.data[:,:,0,:] 
                                  + Data.data[:,:,1,:], 0)
            cal_sum_unscaled += ma.sum(Data.data[:,:,0,:] +
                    Data.data[:,:,1,:], 0)
            # Time series of the cal temperture.
            cal_time = sp.concatenate((cal_time, ma.mean(Data.data[:,:,0,:]
                - Data.data[:,:,1,:], -1).filled(-1)), 0)
            # Everything else done in cal units.
            cal_scale.scale_by_cal(Data)
            # Time serise of the system temperture.
            sys_time = sp.concatenate((sys_time, ma.mean(Data.data[:,:,0,:]
                + Data.data[:,:,1,:], -1).filled(-5)), 0)
            # Accumulate variouse sums.
            counts += this_count
            cal_sum += ma.sum(Data.data[:,:,0,:] + Data.data[:,:,1,:], 0)
            # Take power spectrum of on-off/on+off.
            rebin_freq.rebin(Data, 512, mean=True, by_nbins=True)
            cal_diff = ((Data.data[:,[0,-1],0,:] 
                         - Data.data[:,[0,-1],1,:])
                        / (Data.data[:,[0,-1],0,:] 
                           + Data.data[:,[0,-1],1,:]))
            cal_diff -= ma.mean(cal_diff, 0)
            cal_diff = cal_diff.filled(0)[0:nt,...]
            power = abs(fft.fft(cal_diff, axis=0)[range(nt//2+1)])
            power = power**2/nt
            cal_noise_spec += power
        # Normalize.
        cal_sum_unscaled /= 2*counts
        cal_sum /= 2*counts
        # Get time steps and frequency wdith for noise power normalization.
        Data = Blocks[0]
        Data.calc_time()
        dt = abs(sp.mean(sp.diff(Data.time)))
        # Note that Data was rebined in the loop.
        dnu = abs(Data.field["CDELT1"])
        cal_noise_spec *= dt*dnu/len(Blocks)
        # Power spectrum independant axis.
        ps_freqs = sp.arange(nt//2 + 1, dtype=float)
        ps_freqs /= (nt//2 + 1)*dt*2
        # Long time axis.
        t_total = sp.arange(cal_time.shape[0])*dt
        # Make plots.
        h = plt.figure(figsize=(10,10))
        # Unscaled temperature spectrum.
        plt.subplot(3, 2, 1)
        plt.plot(f/1e6, sp.rollaxis(cal_sum_unscaled, -1))
        plt.xlim((7e2, 9e2))
        plt.xlabel("frequency (MHz)")
        plt.title("Temperture spectrum")
        # Temperture spectrum in terms of noise cal. 4 Polarizations.
        plt.subplot(3, 2, 2)
        plt.plot(f/1e6, sp.rollaxis(cal_sum, -1))
        plt.ylim((-10, 40))
        plt.xlim((7e2, 9e2))
        plt.xlabel("frequency (MHz)")
        plt.title("Temperture spectrum in cal units")
        # Time serise of cal T.
        plt.subplot(3, 2, 3)
        plt.plot(t_total, cal_time)
        plt.xlim((0,dt*3500))
        plt.xlabel("time (s)")
        plt.title("Cal time series")
        # Time series of system T.
        plt.subplot(3, 2, 4)
        plt.plot(t_total, sys_time)
        plt.xlim((0,dt*3500))
        plt.xlabel("time (s)")
        plt.ylim((-4, 35))
        plt.title("System time series in cal units")
        # XX cal PS.
        plt.subplot(3, 2, 5)
        plt.loglog(ps_freqs, cal_noise_spec[:,0,:])
        plt.xlim((1.0/60, 1/(2*dt)))
        plt.ylim((1e-1, 1e3))
        plt.xlabel("frequency (Hz)")
        plt.title("XX cal power spectrum")
        # YY cal PS.
        plt.subplot(3, 2, 6)
        plt.loglog(ps_freqs, cal_noise_spec[:,1,:])
        plt.xlim((1.0/60, 1/(2*dt)))
        plt.ylim((1e-1, 1e3))
        plt.xlabel("frequency (Hz)")
        plt.title("YY cal power spectrum")
        # Adjust spacing.
        plt.subplots_adjust(hspace=.4)
        # Save the figure.
        plt.savefig(params['output_root'] + middle
                + params['output_end'])


# This manager scripts together all the operations that need to be performed to
# the data at GBT.

params_init_manager = {
                      "default_guppi_input_root" : "",
                      "fits_log_root" : "",
                      "output_root" : "",
                      "archive_root" : "",
                      "quality_check_root" : "",
                      "sessions_to_archive" : (),
                      "sessions_to_force" : (),
                      "sessions" : [],
                      "log_file" : "",
                      "error_file" : "",
                      "dry_run" : False
                      }

class DataManager(object) :

    def __init__(self, data_log, feedback=2) :    
        # Read the data log.
        self.params = parse_ini.parse(data_log, params_init_manager, 
                                      prefix="")

    def execute(self, nprocessors=1) :
        params = self.params
        # Redirect the standart in and out if desired.
        if params["log_file"] :
            self.f_log = open(params["log_file"], 'a')
            old_out = sys.stdout
            sys.stdout = self.f_log
        else :
            self.f_log = None
        if params["error_file"] :
            self.f_err = open(params["error_file"], 'a')
            old_err = sys.stderr
            sys.sdterr = self.f_err
        else :
            self.f_err = None
        # Loop over the sessions and process them.
        try :
            for session in params["sessions"] :
                self.process_session(session)
        finally :
            if params["log_file"] :
                sys.stdout = old_out
                self.f_log.close()
            if params["error_file"] :
                sys.stderr = old_err
                self.f_err.close()

    def process_session(self, session) :
        params = self.params
        # Get all the files names and such.
        number = session["number"]
        try :
            guppi_root = session["guppi_input_root"]
        except KeyError :
            guppi_root = params["default_guppi_input_root"]
        guppi_dir = guppi_root + session["guppi_dir"] + "/"
        fits_root = params["fits_log_root"] + str(number) + "/"
        outroot = params["output_root"]
        print ("Processing sesson " + str(number) + ", in guppi directory "
                + guppi_dir)
        # Check that all the directories are present.
        if (not os.path.isdir(guppi_dir) or not os.path.isdir(fits_root)) :
            print "Skipped do to missing data."
            return
        # First thing we do is start to rsync to the archive.
        if number in params['sessions_to_archive'] and not params["dry_run"]:
            SyncProc = subprocess.Popen("rsync -av -essh " + guppi_dir + " " 
                    + params["archive_root"] + session["guppi_dir"] + "/",
                    bufsize=4096, shell=True, stdout=self.f_log, 
                    stderr=self.f_err)
        # Check if this session is in the list of force even if output file
        # exist.
        if number in params['sessions_to_force'] :
            force_session = True
        else :
            force_session = False
        # Now loop over the fields for this session.
        for source in session["sources"] :
            field = source[0]
            # scans is a list of integers containing at least one scan number 
            # from each process.
            scans = source[1]
            # File name pattern that output files should match.
            converted_pattern  = (params["output_root"] + str(number)
                            + '_' + field + '*.fits')
            if force_session :
                scans_to_convert = scans
            else :
                scans_to_convert = []
                for scan in scans :
                    # Check if the scan has already been convereted to SDfits
                    # by looking for the output file.
                    matched_file = check_file_matching_scan(scan, 
                                                        converted_pattern)
                    if not matched_file is None :
                        print ("Skipped because scan already converted "
                               "in file: "
                                + utils.abbreviate_file_path(matched_file))
                    else :
                        scans_to_convert.append(scan)
            if len(scans_to_convert) > 0 :
                # Set up the converter for the scans that havn't been converted
                # yet.
                # Probably shouldn't hard code these.
                converter_params = {"guppi_input_end" : "_0001.fits",
                                    "combine_map_scans" : True,
                                    "partition_cal" : True,
                                    "time_bins_to_average" : 128
                                    }
                # The parameters that are acctually assigned dynamically.
                converter_params["scans"] = tuple(scans_to_convert)
                converter_params["guppi_input_root"] = (guppi_dir + "guppi_"
                        + session["guppi_session"] + "_" + field + "_")
                converter_params["fits_log_dir"] = fits_root
                converter_params["output_root"] = outroot
                # Convert the files from this session.
                if not params["dry_run"] :
                    Converter(converter_params).execute()
            # Check that the new files are in fact present.
            # Make a list of these files for the next section of the pipeline.
            if params["dry_run"] :
                # The files we need to set up the rest of the proceedure will
                # be missing, just return.
                return
            file_middles = []
            for scan in scans :
                converted_file = check_file_matching_scan(scan,
                        converted_pattern)
                if converted_file is None :
                    raise RuntimeError("Scan not converted: " + str(scan))
                # Get the middle (unique) part of the file name.
                converted_file = converted_file.split('/')[-1]
                converted_file = converted_file.split('.')[0]
                file_middles.append(converted_file)
            # Now figure out which scans have had data quality checks
            # already run.
            if force_session :
                files_to_check = file_middles
            else :
                files_to_check = []
                for file_middle in file_middles :
                    # Figure out which files have already been checked by
                    # looking for the output file.
                    checked_out_fname = (params["quality_check_root"] +
                            file_middle + ".pdf")
                    if os.path.isfile(checked_out_fname) :
                        print ("Skipped because file already checked "
                            "in file: " + utils.abbreviate_file_path(
                            checked_out_fname))
                    else :
                        files_to_check.append(file_middle)
            if len(files_to_check) != 0 :
                # Build up input parameters for DataChecker.
                checker_params = {
                                  "output_end" : ".pdf",
                                  "file_middles" : tuple(files_to_check),
                                  "input_root" : params["output_root"],
                                  "output_root" : params["quality_check_root"]
                                  }
                # Execute the data checker.
                DataChecker(checker_params).execute()
            # Wait for the rsync to terminate:
            if number in params['sessions_to_archive']:
                SyncProc.wait()
                if SyncProc.returncode :
                    raise RuntimeError("rsync failed with exit code: " 
                                       + str(SyncProc.returncode))

def check_file_matching_scan(scan, match_str) :
    """
    Looks for a file that matches a provided pattern and contains given scan.

    Parameters
    ----------
    scan : integer
        Scan numbr to look for.
    match_str : string
        Look for file names that match this string.  Assumed to be of the form:
        "*_<lo_scan>-<hi_scan>.extension", where <lo_scan> and <hi_scan> are
        integers and the file name matchs if <lo_scan> <= `scan` and 
        <hi_scan> >= `scan`.

    Returns
    -------
    file_name : string
        The file name of the matching file.  None if no file found.
    """

    file_list = glob.glob(match_str)

    for file_name in file_list :
        # Remove the '.fits'.
        s = file_name.split('.')[-2]
        # Get the scan range part.
        s = s.split('_')[-1]
        lo, hi = s.split('-')
        lo = int(lo)
        hi = int(hi)
        # See if our scan is in that range.
        if scan >= lo and scan <= hi :
            return file_name
    return None


        
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        Converter(str(sys.argv[1])).execute()
    elif len(sys.argv) == 3:
        if sys.argv[1] == str("auto") :
            DataManager(str(sys.argv[2])).execute()
    elif len(sys.argv) == 1:
        Converter().execute()
    else :
        print ("Usage : python psrfits_to_sdfits.py [input file] or"
               " python psrfits_to_sdfits.py auto [data log file]")


