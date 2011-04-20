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

import time
import multiprocessing as mp
import sys
import os.path

import scipy as sp
import scipy.interpolate as interp
import pyfits

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
                p = mp.Process(target=self.ProcessFile, args=(scan, P_far))
                p.start()
                Data = P_here.recv()
                p.join()
                
                # Store our processed data.
                if not Data is None :
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
            p = mp.Process(target=out_write, args=(Block_list, out_file))
            p.start()
            del Block_list
            p.join()
        # End loop over maps (output files or input 'scans' parameter).
            
    def ProcessFile(self, scan, Pipe) :
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
                print self.initial_scan, scan
                print go_hdu["PROCSEQN"], self.initial_scan_ind
                raise ce.DataError("Scans don't agree on proceedure "
                                   "sequence. Aborted scan?  Scans: "
                                   + str(scan) + " and " + 
                                   str(self.initial_scan))

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
        # Sometimes guppi files are just missing.
        if not os.path.isfile(guppi_file) :
            print "Missing psrfits file: " + guppi_file
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
        # Figure out the sample periode after rebinning.
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
                cal_period = 1.0/psrmain["CAL_FREQ"]
                n_bins_cal = cal_period/sample_time
                if n_bins_cal%1.0 > 0.00001 :
                    raise NotImplementedError("Need an integer number of "
                                              "samples per cal periode.")
                n_bins_cal = int(round(n_bins_cal))
                if n_bins_ave%n_bins_cal != 0 :
                    raise ValueError("You should rebin on a multiple of "
                                     "the cal period.  Change parameter "
                                     "time_bins _to_average.")
                # We can speed up the separating of the cal by rebinning in
                # time first.
                tmp_n_bins_ave = n_bins_ave//n_bins_cal
                data.shape = (n_bins, tmp_n_bins_ave, n_bins_cal, npol,
                              nfreq)
                # Mean here instead of median for speed.
                data = scls*sp.mean(data, 1) + offs
                data.shape = (n_bins*n_bins_cal, npol, nfreq)
                # Finally devided into cal on and cal off.
                data = separate_cal(data, n_bins_cal)
            else :
                # Just rebin in time and add a length 1 cal index.
                data.shape = (n_bins, n_bins_ave, npol, nfreq)
                data = sp.median(data, 1)
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
            if sp.any(hours >= 24) :
                # TODO
                raise NotImplementedError("Date rollover in scan.")
            else :
                # sp.char.mod not included in scipy 0.7.
                #time = char.mod(time_string, (hours, minutes, seconds))
                for mm, kk in enumerate(xrange(jj*n_bins, (jj+1)*n_bins)) :
                    tmp_time = time_string%(hours[mm], minutes[mm],
                                            seconds[mm])
                    Data.field["DATE-OBS"][kk] = tmp_time
        Data.add_history('Converted from a guppi sub-integration fits'
                         + ' file.',
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
    folded_data = sp.median(folded_data, 0)
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

        # Find cal on and cal off averages.
        out_data[:,:,0,:] = sp.median(data[:,on_mask,:,:], 1)
        out_data[:,:,1,:] = sp.median(data[:,off_mask,:,:], 1)

    return out_data

        
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        Converter(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        Converter().execute()

