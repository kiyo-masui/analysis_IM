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
import time as tm

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
from utils import misc

params_init = {# Inputs.
               # Everything up to the scan number.
               "guppi_input_roots" : ["./"],
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
               "time_bins_to_average" : 1,
               # Scans to not convert.
               "blacklist" : [],
               # Approximate amount of data to fold to find the cal phase.
               # Each section of data this long will be folded and the cal
               # phase independanly solved for.  In seconds.
               # This should shorter than the time it takes for guppi and the
               # noise cal to change relative phase significantly.  Empirically
               # this seems to be about a minute for guppi samples of 0.001024s
               # and cal period of 64 times this.
               "cal_fold_time" : 30
               }
prefix = ''

class ScanSet(object):

    def __init__(self, fits_data_dir, initial_scan, get_set=True, blacklist=[],
                 feedback=2):
        """Gets the scan log file and figures out what scans belong in this
        set."""

        self.feedback = feedback
        self.fits_data_dir = fits_data_dir
        # Get the fits log file.
        scan_log_hdulist = pyfits.open(fits_data_dir + "/ScanLog.fits", "readonly")
        # From the header we need the project session.
        self.session = scan_log_hdulist[0].header["PROJID"].split('_')[-1]
        # Most of the information we need is in the data extension.
        self.scan_log = scan_log_hdulist[1].data
        scan_log = self.scan_log
        
        # Open the go fits file.
        scan_log_files = scan_log.field('FILEPATH')[
            scan_log.field('SCAN')==initial_scan]
        go_file = fits_data_dir + get_filename_from_key(scan_log_files, "GO")
        go_hdu = pyfits.open(go_file)[0].header
        # From the go information get the source and the scan type.
        self.object = go_hdu["OBJECT"].strip()
        self.proceedure = go_hdu["PROCNAME"].strip().lower()
        if get_set:
            # Read the go file and figure out all the scans in this set.
            # Check the go files for all scans make sure everything is
            # consistant.
            n_scans_proc = go_hdu["PROCSIZE"]
            # Which scan this is of the sequence (1 indexed).
            initial_scan_ind = go_hdu["PROCSEQN"]
            scan_set = (sp.arange(n_scans_proc, dtype=int) + 1 
                               - initial_scan_ind + initial_scan)
            self.scan_set = list(scan_set)
        else :
            self.scan_set = [initial_scan]
        self.apply_blacklist(blacklist)
        self.scan_set.sort()

    def apply_blacklist(self, blacklist):
        """Discards any scan in the blacklist from this set."""

        for scan in list(self.scan_set):
            if scan in blacklist:
                self.scan_set.remove(scan)

    def prepare_scans(self, guppi_data_roots, guppi_data_extension=''):
        """Creates a Scan object and finds file names for each scan in set."""

        self.Scan_objects = []
        proc_info = None
        if self.feedback > 1:
            print ("Searching for guppi files with scan numbers: "
                    + repr(self.scan_set) + ", in directories: "
                    + repr(guppi_data_roots))
        for scan in self.scan_set:
            # Create the scan objects that will process the data.
            # Missing guppi file should not crash program, but should
            # print a warning.
            try:
                Scan_object = Scan(self.fits_data_dir, guppi_data_roots, scan,
                                   guppi_data_extension, feedback=self.feedback)
            except ce.DataError:
                print "Warning, did not find a file for scan: " + repr(scan)
                continue
            # Check some basic info that should be the same for all the scans.
            if not proc_info:
                proc_info = Scan_object.get_proc_info()
            else:
                if proc_info != Scan_object.get_proc_info():
                    msg = ("Scans in same proceedure have conflicting"
                           " parameters.  Perhapse a schduling block was"
                           " canceled.")
                    print proc_info
                    print Scan_object.get_proc_info()
                    raise ce.DataError(msg)
            self.Scan_objects.append(Scan_object)

    def convert_scans(self, time_bins_to_average, partition_cal=False,
                cal_fold_time=0):
        """Converts scans to SDfits."""
        
        self.DataBlocks = []
        for Scan in self.Scan_objects:
            self.DataBlocks.append(Scan.convert(time_bins_to_average,
                                   partition_cal, cal_fold_time))

    def write_out(self, froot):
        """Writes converted data to disk."""
    
        # Figure out the scan range.
        if len(self.scan_set) > 1 :
            str_scan_range = (str(self.scan_set[0]) + '-' +
                              str(self.scan_set[-1]))
        else :
            str_scan_range = str(self.scan_set[0])
        fname = (froot + self.session + '_' + self.object + '_' +
                    self.proceedure + '_' + str_scan_range + '.fits')
        Writer = fitsGBT.Writer(self.DataBlocks, feedback=self.feedback)
        Writer.write(fname)



class Scan(object):

    def __init__(self, fits_data_dir, guppi_data_roots, scan_no,
                 guppi_data_extension='', feedback=2):
        """Constructor just figures out the paths of the required files."""
        
        self.scan_no = scan_no
        self.feedback = feedback
        scan_log_hdulist = pyfits.open(fits_data_dir + "/ScanLog.fits", "readonly")
        # From the header we need the project session.
        self.session = scan_log_hdulist[0].header["PROJID"].split('_')[-1]
        # Most of the information we need is in the data extension.
        scan_log = scan_log_hdulist[1].data
        # Find the GO fits file.
        scan_log_files = scan_log.field('FILEPATH')[
            scan_log.field('SCAN') == scan_no]
        self.go_file = fits_data_dir + get_filename_from_key(scan_log_files, "GO")
        # Get some parameters from the go fits file that we will use to
        # identify correct guppi files.
        go_hdu = pyfits.open(self.go_file)[0].header
        self.object = go_hdu["OBJECT"].strip()
        self.proceedure = go_hdu["PROCNAME"].strip().lower()
        self.proc_size = go_hdu['PROCSIZE']
        self.proc_scan_no = go_hdu['PROCSEQN']
        self.scan_start_time = misc.time2float(go_hdu['DATE-OBS'])
        # Find the Antenna fits file.
        self.antenna_file = (fits_data_dir + get_filename_from_key(scan_log_files, 
                                                        "/Antenna"))
        # Find all possible guppi files. `guppi_data_roots` is a list of
        # directories for which we can search for the proper files. Not all
        # matching files will acctually be from this scan so we will have to
        # check the headers for the session number.  In addition
        candidate_guppi_files = []
        for input_root in guppi_data_roots:
            candidate_guppi_files += glob.glob(input_root + '*' + ("%04d" % scan_no)
                                     + '*' + guppi_data_extension)
        # Now open the headers of these candidate files to see which acctually
        # belong to this scan.
        guppi_files = []
        guppi_file_start_times = []
        for guppi_file in candidate_guppi_files:
            # Get the main header.
            # Occational file with raise an IOError.  Protect against
            # this.
            try:
                hdr = pyfits.getheader(guppi_file, 0)
            except IOError:
                print "Warning, failed to open file: " + guppi_file
                continue
            # Check if the guppi file matchs up with the fits files.
            session = hdr['PROJID'].split('_')[-1]
            object = hdr['SRC_NAME']
            scan_start_time = misc.time2float(hdr['DATE-OBS'])
            file_start_time = misc.time2float(hdr['DATE'])
            if session != self.session:
                # Guppi and GO don't agree on session.
                pass
            elif object != self.object:
                # Guppi and GO don't agree on source.
                pass
            elif abs(scan_start_time - self.scan_start_time) > 5:
                # Guppi and GO don't agree on the scan start time to within 5s.
                pass
            else:
                # This is a matching file.
                guppi_files.append(guppi_file)
                guppi_file_start_times.append(file_start_time)
        # Make sure there are some files left.
        if not guppi_files:
            raise ce.DataError("No guppi files to convert.")
        # If there are multiple guppi files in this scan, sort them by file
        # creation time.
        times_and_files = zip(guppi_file_start_times, guppi_files)
        times_and_files.sort()
        sorted_times, self.guppi_files = zip(*times_and_files)
        if self.feedback >= 1:
            msg = ("Using antenna file: " + self.antenna_file 
                   + " and guppi files: " + repr(self.guppi_files))
            print msg

    def get_proc_info(self):

        proc_info = {}
        proc_info['session'] = self.session,
        proc_info['object'] = self.object
        proc_info['proceedure'] = self.proceedure
        scans_this_file = (sp.arange(self.proc_size, dtype=int) + 1 
                           - self.proc_scan_no + self.scan_no)
        scans_this_file = list(scans_this_file)
        proc_info['proc_scans'] = scans_this_file

        return proc_info
        
    def get_antenna_data(self):
        """Returns an interpoation object containing all the pointing information
        for this file."""
        
        if hasattr(self, 'antenna_data'):
            pass
        else:
            # Open the antenna fits file.
            antenna_hdu_list = pyfits.open(self.antenna_file, 'readonly')
            ant_main = antenna_hdu_list[0].header
            ant_header = antenna_hdu_list[2].header
            ant_data = antenna_hdu_list[2].data
            # Get the time of the scan midpoint in seconds since 00:00 UTC.
            # Times in Julian Days.
            ant_times = sp.array(ant_data.field("DMJD"), dtype=float)
            # Times in days since midnight UTC on day of scan start. Time in this
            # array become greater than one if scan spans midnight UTC.
            ant_times -= (ant_times[0] - ant_times[0]%1.0)
            # Seconds.
            ant_times *= 24.0 * 3600.0
            # Get the az and el out of the antenna.
            ant_az = ant_data.field("OBSC_AZ")
            ant_el = ant_data.field("OBSC_EL")
            az_interp = interp.interp1d(ant_times, ant_az)
            el_interp = interp.interp1d(ant_times, ant_el)
            self.antenna_data =  (az_interp, el_interp, (min(ant_times),
                                  max(ant_times)))
        return self.antenna_data


    def convert(self, time_bins_to_average, partition_cal=False,
                cal_fold_time=0):
        """Converts guppi scan to SDfits."""
        
        Blocks = []
        for guppi_file in self.guppi_files:
            Data = self.process_guppi_file(guppi_file, time_bins_to_average,
                                           partition_cal, cal_fold_time)
            Blocks.append(Data)
        Blocks = tuple(Blocks)
        if len(Blocks) == 1:
            return Blocks[0]
        else:
            return self.merge_blocks(Blocks)

    def merge_blocks(self, Blocks):
        """Merge multiple DataBlock objects from the same guppi scan."""

        # Need to stitch the DataBlocks together.
        # On the first pass, check consistancy and get time axis sizes.
        first_Block = True
        n_time_blocks = []
        for Data in Blocks:
            # Get the number of time entries in this block.
            n_time_blocks.append(Data.dims[0])
            Data.calc_time()
            time_axis = Data.time
            time_delta = abs(sp.mean(sp.diff(time_axis)))
            if first_Block:
                # Get a bunch of stuff that should be the same in all
                # blocks.
                back_shape = Data.dims[1:]
                scan = Data.field['SCAN']
                crval4 = Data.field['CRVAL4']
                crval1 = Data.field['CRVAL1']
                crpix1 = Data.field['CRPIX1']
                cdelt1 = Data.field['CDELT1']
                bandwid = Data.field['BANDWID']
                cal = Data.field['CAL']
                object = Data.field['OBJECT']
                exposure = Data.field['EXPOSURE']
                first_Block = False
            else:
                # Check that this block begins one sample after the last
                # block ended.
                if abs(next_time - time_axis[0]) > 0.003: # ms tolerance.
                    msg = "Time axis not perfectly contiguouse."
                    raise ce.DataError(msg)
                # Check a bunch of stuff that should be identical to the
                # first Block.
                if Data.dims[1:] != back_shape:
                    msg = "Shapes incompatible."
                    raise ce.DataError(msg)
                if Data.field['SCAN'] != scan:
                    msg = "Scan number not compatible"
                    raise ce.DataError(msg)
                if not sp.all(Data.field['CRVAL4'] == crval4):
                    msg = "Polarizations incompatible."
                    raise ce.DataError(msg)
                if Data.field['CRVAL1'] != crval1:
                    msg = "Frequency centre not compatible"
                    raise ce.DataError(msg)
                if Data.field['CRPIX1'] != crpix1:
                    msg = "Frequency centre pixel not compatible"
                    raise ce.DataError(msg)
                if Data.field['CDELT1'] != cdelt1:
                    msg = "Frequency pixel width not compatible"
                    raise ce.DataError(msg)
                if Data.field['BANDWID'] != bandwid:
                    msg = "Bandwidth width not compatible"
                    raise ce.DataError(msg)
                if not sp.all(Data.field['CAL'] == cal):
                    msg = "Cal not compatible."
                    raise ce.DataError(msg)
                if Data.field['OBJECT'] != object:
                    msg = "Object not compatible."
                    raise ce.DataError(msg)
                if Data.field['EXPOSURE'] != exposure:
                    msg = "Exposure not compatible."
                    raise ce.DataError(msg)
            # To make sure the data files are perfectly contiguouse,
            # calculate the expected next time bin.
            next_time = time_axis[-1] + time_delta
        # Allowcate memory for this whole scan.
        # final_shape is ordered (ntime, npol, ncal, nfreq) where ntime
        # is after rebining and combing the records.
        n_time = sp.sum(n_time_blocks)
        final_shape = (n_time,) + back_shape
        d = ma.empty(final_shape, dtype=float)
        Data = data_block.DataBlock(d, copy=False)
        # Allowcate memory for the the fields that are time dependant.
        Data.set_field('CRVAL2', sp.empty(final_shape[0]),
                       ('time',), '1D')
        Data.set_field('CRVAL3', sp.empty(final_shape[0]),
                       ('time',), '1D')
        Data.set_field('DATE-OBS', sp.empty(final_shape[0], dtype='S23'),
                       ('time',), '23A')
        # Copy all the other fields over.
        Data.set_field('CRVAL1', crval1, (), '1D')
        Data.set_field('CRPIX1', crpix1, (), '1I')
        Data.set_field('SCAN', scan, (), '1I')
        Data.set_field('BANDWID', bandwid, (), '1D')
        Data.set_field('CDELT1', cdelt1, (), '1D')
        Data.set_field('OBJECT', object, (), '32A')
        Data.set_field('EXPOSURE', exposure, (), '1D')
        Data.set_field('CRVAL4', crval4, ('pol',), '1I')
        Data.set_field('CAL', cal, ('cal',), '1A')
        Data.verify()
        # Loop though and copy the time dependant feilds.
        this_time_start = 0
        for ii, ThisData in enumerate(Blocks):
            this_n_time = n_time_blocks[ii]
            this_time_slice = slice(this_time_start, this_time_start
                                    + this_n_time)
            Data.data[this_time_slice,...] = ThisData.data
            Data.field['CRVAL2'][this_time_slice] = ThisData.field['CRVAL2']
            Data.field['CRVAL3'][this_time_slice] = ThisData.field['CRVAL3']
            Data.field['DATE-OBS'][this_time_slice] = ThisData.field['DATE-OBS']
            this_time_start += this_n_time
        Data.history = data_block.merge_histories(*Blocks)
        return Data

    def process_guppi_file(self, guppi_file, time_bins_to_average,
                           partition_cal=False, cal_fold_time=0):
        """Converts a guppi_file to DataBlock object."""

        # Whethar we will fold the data to look for the cal or not.
        if partition_cal:
            ncal = 2
        else :
            ncal = 1
        # Open up the fits file.  Be sure to memory map it in case it doesn't
        # fit in memory.
        if self.feedback >= 1:
            print "Converting guppi file: " + guppi_file
        psrhdu_list = pyfits.open(guppi_file, 'readonly', memmap=True)
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
        n_bins_ave = time_bins_to_average
        if ntime%n_bins_ave != 0 :
            raise ValueError("Number of time bins to average must divide "
                             "the number of time bins in a sub integration.")
        # Number of time bins AFTER rebinning
        n_bins = ntime//n_bins_ave
        # Figure out the file start time
        start_string = psrmain["DATE-OBS"] # UTC string including date.
        # Python date_time object for the beginning of the day (00h00 UTC).
        day_start = datetime.datetime.strptime(start_string.split('T')[0],
                                               "%Y-%m-%d")
        # Scan start time since 00h00 UTC.
        start_seconds = psrmain["STT_SMJD"] + psrmain["STT_OFFS"]
        # Figure out the sample period after rebinning.
        sample_time = psrheader["TBIN"]
        resampled_time = sample_time*n_bins_ave
        # Now get the pointing data from the antenna fits file.
        az_interp, el_interp, ant_time_range = self.get_antenna_data()
        
        # Occationally, the guppi integrates longer than the scan and the
        # antenna fits file won't cover the first or last few records.
        # Check for this and trunkate the rocords appropriately.
        n_records = len(psrdata)
        start_record = 0
        record_offsets = psrdata.field("OFFS_SUB")
        for ii in range(3) :
            earliest = (start_seconds + record_offsets[start_record]
                        - sample_time*ntime/2)
            latest = (start_seconds + record_offsets[n_records-1]
                      + sample_time*ntime/2)
            if earliest < ant_time_range[0] :
                start_record -= 1
                print "Discarded first record due to antenna time mismatch."
            if latest > ant_time_range[1] :
                n_records -= 1
                print "Discarded last record due to antenna time mismatch."

        # Allowcate memory for this whole scan.
        # final_shape is ordered (ntime, npol, ncal, nfreq) where ntime
        # is after rebining and combing the records.
        final_shape = (n_bins*(n_records - start_record), npol, ncal,
                              nfreq)
        d = ma.empty(final_shape, dtype=float)
        Data = data_block.DataBlock(d, copy=False)
        # Allowcate memory for the the fields that are time dependant.
        Data.set_field('CRVAL2', sp.empty(final_shape[0]),
                       ('time',), '1D')
        Data.set_field('CRVAL3', sp.empty(final_shape[0]),
                       ('time',), '1D')
        Data.set_field('DATE-OBS', sp.empty(final_shape[0], dtype='S23'),
                       ('time',), '23A')

        # Copy as much field information as we can without looping.
        crpix1 = nfreq//2 + 1
        crval1 = psrdata.field("DAT_FREQ")[0, crpix1 - 1]*1e6
        Data.set_field('CRVAL1', crval1, (), '1D')
        Data.set_field('CRPIX1', crpix1, (), '1I')
        Data.set_field('SCAN', self.scan_no, (), '1I')
        Data.set_field('BANDWID', psrmain['OBSBW']*1e6, (), '1D')
        Data.set_field('CDELT1', psrheader['CHAN_BW']*1e6, (), '1D')
        Data.set_field('OBJECT', psrmain['SRC_NAME'], (), '32A')
        Data.set_field('EXPOSURE', resampled_time/ncal, (), '1D')
        if psrheader["POL_TYPE"].strip() == 'IQUV' :
            Data.set_field('CRVAL4', [1,2,3,4], ('pol',), '1I')
        else :
            raise ce.DataError("Unsupported polarizations.")
        if partition_cal :
            Data.set_field('CAL', ['T', 'F'], ('cal',), '1A')
        else :
            Data.set_field('CAL', ['F'], ('cal',), '1A')
        Data.verify()
        
        # We will devide the data into cal_fold_time chunks and
        # solve for the cal phase in each one.  Figure out how many records go
        # in each chunk.
        # Figure out the number of time bins in a cal period.
        if partition_cal:
            try :
                cal_period = 1.0/psrmain["CAL_FREQ"]
            except ZeroDivisionError :
                # Default value that we usually use.
                cal_period = 0.065536
            n_bins_cal = cal_period/sample_time
            if n_bins_cal%1.0 > 0.00001 :
                raise NotImplementedError("Need an integer number of "
                                          "samples per cal period.")
            n_bins_cal = int(round(n_bins_cal))
            if n_bins_ave%n_bins_cal != 0 :
                raise ValueError("You should rebin on a multiple of "
                                 "the cal period.  Change parameter "
                                         "time_bins_to_average.")
            total_time = (n_records - start_record) * ntime * sample_time
            if cal_fold_time:
                n_divisions = total_time / cal_fold_time
                n_divisions = int(round(n_divisions))
                n_records_div = float(n_records - start_record)/n_divisions
                division_starts = [int(round(start_record + n_records_div*ii))
                                   for ii in range(n_divisions+1)]
            else:
                division_starts = range(start_record, n_records+1)

        # Loop over the records and modify the data.
        if self.feedback > 3 :
            print "Record: ",
        for jj in range(start_record, n_records) :
            record = psrdata[jj]
            # Print the record we are currently on and force a flush to
            # standard out (so we can watch the progress update).
            if self.feedback > 3 :
                print jj,
                sys.stdout.flush()
            # Do one pass through the current division of the data to find the
            # folding offsets.
            if partition_cal and jj in division_starts :
                div_num = division_starts.index(jj)
                # Loop over the data and accumulate it into one array profile.
                record_profile = sp.zeros((ntime, nfreq), dtype=float)
                for kk in range(division_starts[div_num],
                                division_starts[div_num+1]):
                    record_tmp = psrdata[kk]
                    # Data comes from telescope as wierd unsigned and signed 
                    # ints.  Handle carefully. Only interested in stokes I.
                    tmp_data = record_tmp["DATA"]
                    tmp_data.shape = (ntime, npol, nfreq)
                    # The 'I' polarization should already be correctly formated
                    # as an uint8.
                    tmp_data = tmp_data[:,0,:]
                    # There are actually time intependant a frequency dependant
                    # scalings that are ignored here.  They shouldn't affect
                    # the shape of the profile much.
                    record_profile += tmp_data
                # `get_cal_mask` expects 3D data.
                record_profile.shape = (ntime, 1, nfreq)
                # Try to get the profile of this division.  If we fail, set the
                # Phase info to None, which cause the phase to be solved for on
                # a record by record basis.
                try :
                    cal_phase_info = get_cal_mask(record_profile, n_bins_cal)
                except ce.DataError :
                    cal_phase_info = None
            # Get the data field.
            data = record["DATA"]
            # Convert the Data to the proper shape and type.
            data = format_data(data, ntime, npol, nfreq)

            if partition_cal :
                data = separate_cal(data, n_bins_cal, cal_phase_info)
                # Down sample from the cal period to the final number of time
                # bins.
                data.shape = (n_bins, ntime//n_bins_cal//n_bins, npol,
                              ncal, nfreq)
                # Use mean not median due to discritization.
                data = sp.mean(data, 1)
            else :
                # Just rebin in time and add a length 1 cal index.
                data.shape = (n_bins, n_bins_ave, npol, ncal, nfreq)
                data = sp.mean(data, 1)
            # Now get the scale and offsets for the data and apply them.
            scls = sp.array(record["DAT_SCL"], dtype=sp.float32)
            scls.shape = (1, npol, 1, nfreq)
            offs = sp.array(record["DAT_OFFS"], dtype=sp.float32)
            offs.shape = (1, npol, 1, nfreq)
            data = scls*data + offs
            
            # Now figure out the timing information.
            time = start_seconds + record["OFFS_SUB"]
            time = time + (sp.arange(-n_bins/2.0 + 0.5, n_bins/2.0 + 0.5)
                           * resampled_time)
            # Make sure that there are no time gaps between records.
            # This is nessisary for the cal separation to be valid.
            if jj > start_record:
                record_time = record["OFFS_SUB"] - last_offs_sub
                if not sp.allclose(record_time, n_bins*resampled_time):
                    msg = "Time gaps between records"
                    raise ce.DataError(msg)
            last_offs_sub = record["OFFS_SUB"] 
            # Get the pointing.
            az = az_interp(time)
            el = el_interp(time)
            # Put all the data into the DataBlock for writing out.
            Data.data[jj*n_bins:(jj+1)*n_bins,:,:,:] = data
            Data.field['CRVAL2'][jj*n_bins:(jj+1)*n_bins] = az
            Data.field['CRVAL3'][jj*n_bins:(jj+1)*n_bins] = el
            # Get the time stamp for each bin.
            for mm, kk in enumerate(xrange(jj*n_bins, (jj+1)*n_bins)) :
                # Format the time and add it to the date.
                this_time_obj = datetime.timedelta(seconds=time[mm])
                this_date_time = day_start + this_time_obj
                # Deal will partial seconds explicitly since there is not
                # builtin function.
                milliseconds = "%06d" % this_date_time.microsecond
                if len(milliseconds) != 6:
                    raise RuntimeError("Huh?")
                # Truncation better than rounding (corner case errors).
                milliseconds = milliseconds[:3]
                this_time_str = (this_date_time.strftime("%Y-%m-%dT%H:%M:%S.") 
                                 + milliseconds)
                Data.field["DATE-OBS"][kk] = this_time_str
        Data.add_history('Converted from a guppi sub-integration fits file.',
                         (utils.abbreviate_file_path(guppi_file),
                          utils.abbreviate_file_path(self.antenna_file),
                          utils.abbreviate_file_path(self.go_file)))
        # End loop over records.
        if self.feedback > 3 :
            print
        return Data

class Converter(object):

    def __init__(self, input_file_or_dict=None, feedback=2) :    
        self.feedback = feedback
        if self.feedback > 1:
            print "Starting Data Converter."
        # Read the parameter file.
        self.params = parse_ini.parse(input_file_or_dict, params_init, 
                                      prefix=prefix)
    
    def execute(self, nprocesses=1) :
        params = self.params
        initial_scans = list(params["scans"])
        # Make sure that the output directory exists.
        utils.mkparents(params["output_root"])
        
        if False: # Alwasys multiprocess to avoid memory leaks.
            # Single process case.
            # Loop though all the listed scans and convert them.
            for scan in initial_scans:
                self.execute_set(scan)
        else:
            # Multi process case.
            n_new = nprocesses
            proc_list = range(n_new)
            n_tasks = len(initial_scans)
            for ii in range(n_new + n_tasks):
                # End the task number ii - n_new.
                if ii >= n_new:
                    proc_list[ii % n_new].join()
                    if proc_list[ii % n_new].exitcode != 0:
                        raise RuntimeError("A thread failed with exit code: "
                                          + str(proc_list[ii%n_new].exitcode))
                # Start task number ii.
                if ii < n_tasks:
                    proc_list[ii % n_new] = mp.Process(
                            target=self.execute_set, args=(initial_scans[ii],))
                    proc_list[ii % n_new].start()

    def execute_set(self, scan):
        params = self.params
        Set = ScanSet(params["fits_log_dir"], scan,
                      get_set=params['combine_map_scans'],
                      blacklist=params['blacklist'], feedback=self.feedback)
        Set.prepare_scans(params['guppi_input_roots'],
                          params['guppi_input_end'])
        Set.convert_scans(params['time_bins_to_average'],
                          params['partition_cal'], params['cal_fold_time'])
        Set.write_out(params['output_root'])


# Functions for the most error prone parts of the above script.  Unit testing
# applied to these.

def get_filename_from_key(file_list, key) :
    found_file = False
    for file_name in file_list :
        if key in file_name :
            file = ("/" + key + file_name.split(key)[1])
            found_file = True
            break
    if not found_file :
        raise ce.DataError(key + " file not in scan log.")
    else :
        return file

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
        msg = ("Number of bins in a cal period must devide the number"
               " of time bins in a subintegration (generally number"
               " of bins in a cal period should be a power of 2).")
        raise ValueError(msg)
    
    # For solving for the cal phase, throw away high varience channels (RFI).
    I_data = data[:,0,:]
    I_vars = sp.var(I_data, 0)
    mean_var = sp.mean(I_vars)
    good_chans = I_vars < 5 * mean_var
    I_data = I_data * good_chans
    # Fold the Stokes I data on the cal period to figure out the phase.
    # Sum over the frequencies.
    folded_data = sp.sum(I_data, -1, dtype=float)
    # Fold the data at the cal period.
    folded_data.shape = ((ntime//n_bins_cal, n_bins_cal) +
                         folded_data.shape[1:])
    folded_data = sp.mean(folded_data, 0)
    # Split the data into cal on and cal offs.
    base = sp.mean(folded_data)
    diff = sp.std(folded_data)
    transition = sp.logical_and(folded_data < base + 0.80*diff,
                                folded_data > base - 0.80*diff)
    if sp.sum(transition) > 2:
        profile_str = repr((folded_data - sp.mean(folded_data))
                           / sp.std(folded_data))
        raise ce.DataError("Cal profile has too many "
                           "transitions.  Profile array: "
                           + profile_str)
    # Find the last bin off before on.
    last_off_ind = None
    for kk in xrange(n_bins_cal) :
        if ((folded_data[kk] < base-0.70*diff) and not
            (folded_data[(kk+1)%n_bins_cal] < base-0.70*diff)):
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
    if folded_data[(last_off_ind+1)%n_bins_cal] < base+0.70*diff :
        n_cal_state = n_bins_cal//2 - 1
        n_blank = 1
    else :
        n_cal_state = n_bins_cal//2 - 2
        n_blank = 2

    long_profile = sp.concatenate((folded_data, folded_data), 0)
    first_off = first_on + n_bins_cal//2
    if (sp.any(long_profile[first_on:first_on+n_cal_state] < base + 0.80*diff)
        or sp.any(long_profile[first_off:first_off+n_cal_state] > base -
                  0.80*diff)) :
        profile_str = repr((folded_data - sp.mean(folded_data))
                           / sp.std(folded_data))
        raise ce.DataError("Cal profile not lining up "
                           "correctly.  Profile array: "
                           + profile_str)

    return first_on, n_blank

def separate_cal(data, n_bins_cal, cal_mask=None) :
    """Function separates data into cal_on and cal off.
    
    No Guarantee that data argument remains unchanged."""
    
    # Allowcate memeory for output    
    ntime, npol, nfreq = data.shape
    n_bins_after_cal = ntime//n_bins_cal
    out_data = sp.zeros((n_bins_after_cal, npol, 2, nfreq), dtype=sp.float32)
    
    # Get the phase offset of the cal.
    try :
        if cal_mask is None:
            first_on, n_blank = get_cal_mask(data, n_bins_cal)
        else :
            first_on, n_blank = cal_mask
    except ce.DataError :
        print "Discarded record due to bad profile. "
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
        # This loop is much faster than the built in numpy mean() for some
        # reason.
        for ii in range(n_bins_cal) :
            if on_mask[ii]:
                out_data[:,:,0,:] += data[:,ii,:,:]
            elif off_mask[ii]:
                out_data[:,:,1,:] += data[:,ii,:,:]
        out_data[:,:,0,:] /= sp.sum(on_mask)
        out_data[:,:,1,:] /= sp.sum(off_mask)

    return out_data


# -------- Data Checker ---------

# Makes a bunch of plots that give an indication of how the data looks.

params_init_checker = {
                      "input_root" : "",
                      "file_middles" : ("",),
                      "output_root" : "",
                      "output_end" : ".pdf"
                      }

class DataChecker(object) :

    def __init__(self, parameter_file_or_dict, feedback=2):
        self.feedback = feedback
        if self.feedback > 1:
            print "Starting Data Checker."
        self.params = parse_ini.parse(parameter_file_or_dict, 
                params_init_checker, prefix="")

    def execute(self, nprocesses=1) :
        params = self.params
        n = len(params["file_middles"])
        np = nprocesses
        procs = [None]*np
        # Loop over files, make one set of plots per file.
        for ii in range(n+np) :
            if ii >= np :
                procs[ii%np].join()
            if ii < n :
                p = mp.Process(target=self.process_file,
                               args=(params["file_middles"][ii],))
                p.start()
                procs[ii%np] = p

    def process_file(self, middle) :
        """Split off to fix pyfits memory leak."""
        params = self.params
        # Construct the file name and read in all scans.
        file_name = params["input_root"] + middle + ".fits"
        Reader = fitsGBT.Reader(file_name, feedback=self.feedback)
        Blocks = Reader.read((), (), force_tuple=True)
        # Plotting limits need to be adjusted for on-off scans.
        if file_name.find("onoff") != -1 :
            onoff=True
        else :
            onoff=False
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
        plt.title("System temperature - mean over time")
        # Temperture spectrum in terms of noise cal. 4 Polarizations.
        plt.subplot(3, 2, 2)
        plt.plot(f/1e6, sp.rollaxis(cal_sum, -1))
        if onoff :
            plt.ylim((-1, 60))
        else :
            plt.ylim((-10, 40))
        plt.xlim((7e2, 9e2))
        plt.xlabel("frequency (MHz)")
        plt.title("System temperature in cal units")
        # Time serise of cal T.
        plt.subplot(3, 2, 3)
        plt.plot(t_total, cal_time)
        if onoff :
            plt.xlim((0,dt*900))
        else :
            plt.xlim((0,dt*3500))
        plt.xlabel("time (s)")
        plt.title("Noise cal temperature - mean over frequency")
        # Time series of system T.
        plt.subplot(3, 2, 4)
        plt.plot(t_total, sys_time)
        plt.xlabel("time (s)")
        if onoff :
            plt.ylim((-4, 90))
            plt.xlim((0,dt*900))
        else :
            plt.ylim((-4, 35))
            plt.xlim((0,dt*3500))
        plt.title("System temperature in cal units")
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


# -------- Data Manager --------

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
                      "rsync_file" : "",
                      "nprocesses": 1,
                      "dry_run" : False
                      }

class DataManager(object) :

    def __init__(self, data_log, feedback=2):
        self.feedback = feedback
        # Read the data log.
        if self.feedback > 1:
            print "Starting Data Manager."
        self.params = parse_ini.parse(data_log, params_init_manager, 
                                      prefix="")

    def execute(self, nprocesses=1) :
        params = self.params
        # Redirect the standart in and out if desired.
        if params["log_file"] :
            self.f_log = open(params["log_file"], 'w')
            old_out = sys.stdout
            sys.stdout = self.f_log
        else :
            self.f_log = None
        if params["error_file"] :
            self.f_err = open(params["error_file"], 'w')
            old_err = sys.stderr
            sys.stderr = self.f_err
        else :
            self.f_err = None
        if params["rsync_file"] :
            self.r_out = open(params["rsync_file"], 'w')
        else :
            self.r_out = None
        # Loop over the sessions and process them.
        sessions = params["sessions"]
        sessions.reverse()
        for session in sessions :
            self.process_session(session)

    def process_session(self, session) :
        params = self.params
        # Get all the files names and such.
        number = session["number"]
        try :
            guppi_root = session["guppi_input_root"]
        except KeyError :
            guppi_root = params["default_guppi_input_root"]
        # Some observations are spread over multiple directories, in which case
        # session["guppi_dir"] will be a list.
        if isinstance(session["guppi_dir"], list) :
            session_dirs = session["guppi_dir"]
        else :
            session_dirs = [session["guppi_dir"]]
        # Version of session_dirs with absolute path.
        guppi_dirs = [guppi_root + dir + "/" 
                      for dir in session_dirs]
        fits_root = params["fits_log_root"] + "%02d"%number + "/"
        outroot = params["output_root"]
        print ("Processing sesson " + str(number) + ", in guppi directories "
                + str(guppi_dirs))
        # Check that all the directories are present.
        for dir in guppi_dirs :
            if (not os.path.isdir(dir) or not os.path.isdir(fits_root)) :
                print "Skipped do to missing data."
                return
        # First thing we do is start to rsync to the archive.
        if number in params['sessions_to_archive'] and not params["dry_run"]:
            # Construct the multiple line command that will tarnsfer all
            # guppi directories to teh archive.
            program = "rsync -av -essh "
            command = ""
            for ii in range(len(session_dirs)) :
                command += program
                command += guppi_dirs[ii] + " "
                command += (params["archive_root"] + str(number) + "_" 
                            + session_dirs[ii] + "/")
                command += ";"
            print command
            SyncProc = subprocess.Popen(command, shell=True,
                           stdout=self.r_out, stderr=subprocess.STDOUT)
        # Check if this session is in the list to force even if output file
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
            converted_pattern  = (params["output_root"] + "%02d"%number
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
                converter_params = {"guppi_input_end" : ".fits",
                                    "combine_map_scans" : True,
                                    "partition_cal" : True,
                                    "time_bins_to_average" : 128,
                                    "cal_fold_time" : 10
                                    }
                # The parameters that are acctually assigned dynamically.
                converter_params["scans"] = tuple(scans_to_convert)
                converter_params["fits_log_dir"] = fits_root
                converter_params["output_root"] = outroot
                try :
                    blacklist = session["blacklist"]
                except KeyError :
                    blacklist = []
                converter_params["blacklist"] = blacklist
                # Guppi input root should be a list of all possible roots for
                # this session.
                converter_params["guppi_input_roots"] = []
                for dir in guppi_dirs :
                    converter_params["guppi_input_roots"].append(dir)
                # Convert the files from this session.
                if not params["dry_run"] :
                    C = Converter(converter_params, feedback=self.feedback)
                    C.execute(params['nprocesses'])
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
                DC = DataChecker(checker_params, feedback=self.feedback)
                DC.execute(params["nprocesses"])
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
        intstrs = s.split('-')
        if len(intstrs) == 2:
            # File contains a range of scans.
            lo = int(intstrs[0])
            hi = int(intstrs[1])
            # See if our scan is in that range.
            if scan >= lo and scan <= hi :
                return file_name
        elif len(intstrs) == 1:
            # File contains one scan.
            if scan == int(intstrs[0]):
                return file_name
    return None


        
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        Converter(str(sys.argv[1])).execute()
    elif len(sys.argv) == 3 and sys.argv[1] == str("auto") :
        DataManager(str(sys.argv[2])).execute()
    elif len(sys.argv) == 1:
        Converter().execute()
    else :
        print ("Usage : python psrfits_to_sdfits.py [input file] or"
               " python psrfits_to_sdfits.py auto [data log file]")

