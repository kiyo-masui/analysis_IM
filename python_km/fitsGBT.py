"""Read a GBT spectrometer FITS file (SDfits) and do some preliminary
processing of the time stream.

The FITS file is assumed to be in a certain format ocrresponding to regular
GBT spectrometer data.
"""

import scipy as sp
import numpy.ma as ma
import pyfits

import kiyopy.custom_exceptions as ce

class Processor() :
    """Class that opens a GBT Spectrometer Fits file, reads data from it and
    does time stream processing.

    This class opens the a GBT Fits File upon initialization and closes it upon
    deletion.  It contains routines for reading individual scans and IFs from
    the file.  This class reads data but does not store data.

    Processing:  When data is read, this class does any processing that can be
    done at a local time.  Hanning smoothing, cal isolation, RFI flagging.

    Arguments:
        fname: Required intialization argument.  FITS file name to be read.  
            The file is assumed to have a certain entries and be arranged a
            certain way corresponding to the GBT spectrometer data.
    """

    # Note to programmer: assignments of a subset of one array to another
    # are often by reference, not by value.  To be safe, any array you are
    # going to modify should be forced to assign by value with the
    # sp.array(an_array) function.

    # These are flags for performing variouse checks.  They will eventually
    # become inputs of some sort.
    verify_ordering = 1;
    verify_keynames = 1;
    verify_lengths = 1;
    verify_bands = 1;
    feedback = 2;

    # Some parameters of GBT spectrometer data.  These are assumed to be
    # correct with only minimal checking.
    npol = 4 # number of polarizations
    ncal = 2 # number of noise cal states (on, off)

    # Flags for internal use:


    # All objects that are of length # have names ending in #:
    # number of records in the fits file --- _all.
    # number of records in one IF and one scan --- _sif.
    # number of records in sif with one polarization/cal --- _tsc
    #   (also number of times in a scan)

    def __init__(self, fname) :
        """
        See class docstring (for now).
        """

        self.fname = fname
        if self.feedback > 0 :
            print "Opened GBT fits file: \n ", fname

        # The passed file name is assumed to be a GBT spectrometer fits file.
        self.hdulist = pyfits.open(self.fname, 'readonly')
        # Separate in to the useful sub objects.  These assignments are all
        # done by reference, so this is efficient.
        self.fitsdata = self.hdulist[1].data
        # The records in fitsdata are not guaranteed to be in proper order.
        # Mostly the IFs are all out of whack.  However, once you isolate an 
        # IF everything should be well ordered.

        # Get the scans and IF of all records.  This is later used to isolate a
        # single IF and scan.  Also get the set of unique IFs and scans, so we
        # know what is in the file.
        self._scans_all = self.fitsdata.field('SCAN')
        self.scan_set = sp.unique(self._scans_all)
        self._IFs_all = self.fitsdata.field('CRVAL1')/1E6 # MHz
        # Round the frequencies as we only need to tell the difference between
        # one IF and the other.
        # TODO: Figure out what the difference frequencies acctually are.
        self._IFs_all = self._IFs_all.round(0) 
        self.IF_set = sp.unique(self._IFs_all)
        # TODO: This is inefficient.  sp.unique has this for numpy > 1.3.
        IF_unique_inds = []
        for tempIF in self.IF_set :
            inds, = sp.where(self._IFs_all == tempIF)
            IF_unique_inds.append(inds[0])

        # Get IF bandwidths.  Require them to be the same for each IF for now.
        IF_bandwidths = self.fitsdata[(IF_unique_inds,)]['BANDWID']/1E6
        IF_bandwidths = IF_bandwidths.round(0)
        self.IF_bandwidth = IF_bandwidths[0]
        if self.verify_bands :
            for tband in IF_bandwidths :
                if tband != self.IF_bandwidth :
                    raise ce.DataError("IF bandwidths not all the same. "
                                    "This needs to be implemented.")
        self.nfreq_IF = len(self.fitsdata[0]['DATA']) # number of frequencies
        # in an IF.  There are overlapping frequencies between IFs, but for
        # now process all.
        # TODO get frequencies.

    def get_scan_IF_inds(self, scan_ind, IF_ind) :
        """Gets the record indices of the fits file that correspond to the
        given scan and IF.

        Note that the scans are numbered with 0 corresponding to the first scan
        in the file i.e., it is not the session scan number."""

        # Should check valid scan IF, and raise value errors as apropriate
        thescan = self.scan_set[scan_ind]
        theIF = self.IF_set[IF_ind]
        
        # Find all the records that correspond to this IF and this scan.
        # These indicies *should now be ordered in time, cal (on off)
        # and in polarization, once the IF is isolated.
        (inds_sif,) = sp.where(sp.logical_and(self._IFs_all==theIF, 
                                        self._scans_all==thescan))
        
        # Reform to organize by pol, cal, etc.
        ntimes = len(inds_sif)//self.npol//self.ncal
        inds_sif = sp.reshape(inds_sif, (ntimes, self.npol, self.ncal))

        if self.verify_ordering > 0:
            # We expect noise cal to be on for every second record.
            for thecal in range(self.ncal) :
                tmp = sp.unique(self.fitsdata.field('CAL')[inds_sif[:,:,thecal]])
                if len(tmp) > 1 :
                    raise ce.DataError("Calibration (ON/OFF) not in "
                                    "perfect order in file: "+self.fname)
            # Polarization should cycle through 4 modes (-5,-7,-8,-6)
            for thepol in range(self.npol) :
                tmp = sp.unique(self.fitsdata.field('CRVAL4')
                            [inds_sif[:,thepol,:]])
                if len(tmp) > 1 :
                    raise ce.DataError("Polarizations not in perfect order in "
                                    "file: "+self.fname)
            # We expect the entries to be sorted in time
            tmp = (self.fitsdata.field('LST')
                       [inds_sif[:,0,0]])
            lastLST = 0
            for thisLST in tmp :
                if lastLST > thisLST :
                    # Rollover at 86400 is normal.
                    if lastLST > 86395 and thisLST < 5 :
                        raise ce.DataError("LST roll over. Normal but not"
                                        " supported yet.") 
                    else :
                        raise ce.DataError("Times not in perfect order.")
                lastLST=thisLST
            del tmp

        return inds_sif


    
    def read_process(self, scans=None, IFs=None) :
        """ Read and process data from the fits file.

        This method reads data from the fits file and does basic time stream
        processing.

        Arguments:
            fname: FITS file name to be read.  The file is assumed to have a
                certain entries and be arranged a certain way corresponding 
                to the GBT spectrometer data.
            scans: Which scans in the file to be processed.  A list of 
                integers, with 0 corresponding to the lowest numbered scan.
                Default is all of them.
            IFs: Which intermediate frequencies (also called frequency windows)
                to process.  A list of integers with 0 coorsponding to the 
                lowest frequency present. Default is all of them.
                TODO : Overlapping frequency windows stiched together somehow.

        Returns: A dictionary filled with data.
            TODO: Fill this in with detailed keys names and such.
        """
        # The default is to read all scans and all IFs.
        if not scans :
            scans = range(len(self.scan_set))
        if not IFs :
            IFs = range(len(self.IF_set))
        if self.feedback > 0 and self.feedback < 3 :
            print "Precessing scans :", scans, " and IFs ", IFs
        
        # TODO
        # First determine which where the IF bands overlap and the indicies
        # that we will keep for each IF.
        # sortedIFs = IFs.sort()
        # acctually, do this later.  We can calibrate to figure this all out.
        nfreq_total = self.nfreq_IF * len(IFs) # temporary
        records_per_time = self.npol * self.ncal

        # Initialize output arrays.  Will use concatenation, so time index is
        # initialized zero lengthed.
        # ma = masked_array.  For flagging out bad data.
        P_total = ma.zeros((0, self.npol, self.ncal, nfreq_total), order='C',
                              dtype=float)
        LST_total = sp.zeros((0,))

        # Now loop over all the scans and IFs and deal with each individually
        for thescan in self.scan_set[scans]:
            first_iteration_IF = True
            for theIF in self.IF_set[IFs] :
               
                # Now we are ready to get a IF/scan set of data and work on it.
                # Masked array (ma) combines flags and data.
                P_IF = ma.empty((ntimes_scan, self.npol, self.ncal, 
                                 self.nfreq_IF), order='C',dtype=float)
                for ii in range(ntimes_scan) :
                    for jj in range(self.npol) :
                        for kk in range(self.ncal) :
                            # XXX A good place for a bug ... test me.
                            index = ii*self.npol*self.ncal + jj*self.ncal + kk
                            P_IF[ii,jj,kk,:] = ( self.fitsdata.field('DATA') 
                                                          [inds_sif[index]] )
                # Pointing and time data only needed once per time.
                if first_iteration_IF :
                    azimuth_scan = ( sp.array(self.fitsdata.field('AZIMUTH') 
                                         [inds_sif[0::records_per_time]]) )
                    elevation_scan = ( sp.array(self.fitsdata.field('ELEVATIO') 
                                         [inds_sif[0::records_per_time]]) )
                    # Local Sidereal Time.
                    LST_scan = sp.array(self.fitsdata.field('LST')
                                   [inds_sif[0::records_per_time]])
                # Make Flag matrix to store locaiton of bad data
                P_IF[sp.logical_not(sp.isfinite(P_IF))] = ma.masked
                if self.feedback > 2 :
                    n_masked = ma.count_masked(P_IF)
                    print "NANs this scan and IF :", n_masked, ", ",
                    print float(n_masked)/float(sp.size(P_IF)), "%."

                # Perform Hanning smoothing frequency space convolution.
                # Masked array automatically flags data adjacent to masked
                # data.
                P_IF[:,:,:,1:-1] = ( 0.25*P_IF[:,:,:,:-2] + 
                                     0.50*P_IF[:,:,:,1:-1] + 
                                     0.25*P_IF[:,:,:,2:] )
                # End points where not smoothed.
                P_IF[:,:,:,0] = ma.masked
                P_IF[:,:,:,-1] = ma.masked
                
                # TODO : RFI flagging.
                
                # Report number of flags.
                if self.feedback > 2 :
                    n_masked = ma.count_masked(P_IF)
                    print "Flagged this scan and IF :", n_masked, ", ",
                    print float(n_masked)/float(sp.size(P_IF)), "%."    
                # For now this is all, just build final array.
                if first_iteration_IF :
                    P_scan = P_IF
                else :
                    P_scan = ma.concatenate((P_scan, P_IF), 3)
                
                first_iteration_IF = False
                # End theIF loop.

            P_total = ma.concatenate((P_total, P_scan), 0)
            LST_total = sp.concatenate((LST_total, LST_scan))
            # End thescan loop.
        
        if self.feedback > 1 :
            n_masked = ma.count_masked(P_total)
            print "Flagged total :", n_masked, ", ",
            print float(n_masked)/float(sp.size(P_total)), "%." 

        return {"P":P_total, "LST":LST_total}


    def __del__(self) :
        self.hdulist.close()
        if self.feedback > 3 :
            print "Closed ", self.fname


