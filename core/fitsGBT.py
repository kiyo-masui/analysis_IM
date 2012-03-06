"""Classes that read and write GBT spectrometer FITS files (SDfits).

The FITS file is assumed to be in a certain format ocrresponding to regular
GBT spectrometer data.
"""

import scipy as sp
import numpy.ma as ma
import pyfits

import kiyopy.custom_exceptions as ce
import kiyopy.utils as ku
import base_fits as bf
import data_block as db

# Here we list the fields that we want to track.  The fields listed below will
# be read and stored every time a fits files is read and will then be written
# when data is written back to fits.

# If I didn't initially think that the field you need should be copied, you can
# still get that data from the initial data file, which should be listed in the
# file history.

# The only field that is allowed to vary over the 'freq' axis right now is
# 'DATA' (not listed below since it is always read).

# If you find that there are other fields you need that are not on this list,
# feel free to add them.  This should not create an problems with other
# people's code.

# For a description of the fields:
#    https://safe.nrao.edu/wiki/bin/view/Main/SdfitsDetails
fields_and_axes = { 
                   'SCAN' : (),
                   'CRVAL1' : (), # Band centre frequency
                   'BANDWID' : (),
                   'CRPIX1' : (), # Center frequency channel
                   'CDELT1' : (), # Frequency Channel width
                   'OBJECT' : (),
                   'TIMESTAMP' : (),
                   'OBSERVER' : (),
                   'RESTFREQ' : (),
                   'DURATION' : (),
                   'EXPOSURE' : ('time', 'cal'),
                   'CTYPE2' : (), # Type of longetudinal axis (probably AZ)
                   'CTYPE3' : (),
                   'DATE-OBS' : ('time', ),
                   'LST' : ('time', ),
                   # These pointings refer to the structural telescope, not
                   # the beam (slightly different).
                   'ELEVATIO' : ('time', ),
                   'AZIMUTH' : ('time', ),
                   'OBSFREQ' : ('time', ),
                   # These pointings are corrected for pointing calibration
                   # and I think refraction.
                   'CRVAL2' : ('time', ), # Azimuth
                   'CRVAL3' : ('time', ), # Elevation
                   # Polarization indicies
                   'CRVAL4' : ('pol', ),
                   'CAL' : ('cal', )
                   }

# These globals are the cards for (our custom) history entries in a fits header
card_hist = 'DB-HIST'
card_detail = 'DB-DET'


class Reader(object) :
    """Class that opens a GBT Spectrometer Fits file and reads data.

    This class opens the a GBT Fits File upon initialization and closes it upon
    deletion.  It contains routines for reading individual scans and IFs from
    the file.  This class reads data but does not store data.  Data is stored
    in DataBlock objects, which are returned by the 'read(scans, IFs)' method.

    Arguments:
        fname: Required intialization argument.  FITS file name to be read.  
            The file is assumed to have a certain entries and be arranged a
            certain way corresponding to the GBT spectrometer data.
    """

    # Note to programmer: assignments of a subset of one array to another
    # are often by reference, not by value.  To be safe, any array you are
    # going to modify should be forced to assign by value with the
    # sp.array(an_array) function.

    def __init__(self, fname, feedback=2, checking=1) :
        """Init script for the fitsGBT Reader class.

        The reader is initialised with the fits file name to be read.
        Optionally, one can supply up to two integers: feedback (default 2)
        indicating how much junk to print, and checking (default 1) specifying
        to what level the input file is checked for sensible data.
        """

        self.feedback = feedback
        self.checking = checking
        if checking > 0 :
            self.verify_ordering = 1
        else :
            self.verify_ordering = 0

        self.fname = fname

        # The passed file name is assumed to be a GBT spectrometer fits file.
        self.hdulist = pyfits.open(self.fname, 'readonly')
        if len(self.hdulist) < 2 :
            raise ce.DataError("File missing data extension")
        if self.feedback > 0 :
            print "Opened GBT fits file: ", ku.abbreviate_file_path(fname)
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
        # Sort scans for reapeatable ordering.
        self.scan_set.sort()
        self._IFs_all = self.fitsdata.field('CRVAL1')/1E6 # MHz
        # Round the frequencies as we only need to tell the difference between
        # one IF and the other.
        self._IFs_all = self._IFs_all.round(0) 
        self.IF_set = sp.unique(self._IFs_all)
        self.IF_set.sort()

    def get_scan_IF_inds(self, scan_ind, IF_ind) :
        """Gets the record indices of the fits file that correspond to the
        given scan and IF.

        Note that the scans are numbered with 0 corresponding to the first scan
        in the file i.e., it is not the session scan number."""

        # TODO: Should check valid scan IF, and raise value errors as apropriate
        thescan = self.scan_set[scan_ind]
        theIF = self.IF_set[IF_ind]
        
        # Find all the records that correspond to this IF and this scan.
        # These indicies *should now be ordered in time, cal (on off)
        # and in polarization, once the IF is isolated.
        (inds_sif,) = sp.where(sp.logical_and(self._IFs_all==theIF, 
                                        self._scans_all==thescan))
        ncal = len(sp.unique(self.fitsdata.field('CAL')[inds_sif]))
        npol = len(sp.unique(self.fitsdata.field('CRVAL4')[inds_sif]))
        
        # Reform to organize by pol, cal, etc.
        ntimes = len(inds_sif)//npol//ncal
        inds_sif = sp.reshape(inds_sif, (ntimes, npol, ncal))

        if self.verify_ordering > 0:
            # We expect noise cal to be on for every second record.
            for thecal in range(ncal) :
                tmp = sp.unique(self.fitsdata.field('CAL')[inds_sif[:,:,thecal]])
                if len(tmp) > 1 :
                    raise ce.DataError("Calibration (ON/OFF) not in "
                                    "perfect order in file: "+self.fname)
            # Polarization should cycle through 4 modes (-5,-7,-8,-6)
            for thepol in range(npol) :
                tmp = sp.unique(self.fitsdata.field('CRVAL4')
                            [inds_sif[:,thepol,:]])
                if len(tmp) > 1 :
                    raise ce.DataError("Polarizations not in perfect order in "
                                    "file: "+self.fname)
            # We expect the entries to be sorted in time and for time to not
            # change across pol and cal.
            lastLST = 0
            for ii in range(ntimes) :
                # Sometimes won't have the LST.
                try :
                    thisLST = self.fitsdata.field('LST')[inds_sif[ii,0,0]]
                # If 'LST' is missing raises a KeyError in later versions of
                # pyfits, and a NameError in earlier ones.
                except (KeyError, NameError) :
                    break
                if not (sp.allclose(self.fitsdata.field('LST')
                        [inds_sif[ii,:,:]] - thisLST, 0)) :
                    raise ce.DataError("LST change across cal or pol in "
                                       "file: " + self.fname)

        return inds_sif

    def set_history(self, Block) :
        """Reads the file history and sets the corresponding Block history."""

        prihdr = self.hdulist[0].header
        # If there is no history, return.
        try :
            ii = prihdr.ascardlist().index_of(card_hist)
        except KeyError :
            return
        n_cards = len(prihdr.ascardlist().keys())
        while ii < n_cards :
            if prihdr.ascardlist().keys()[ii] == card_hist :
                hist_entry = prihdr[ii]
                details = []
            elif prihdr.ascardlist().keys()[ii] == card_detail :
                details.append(prihdr[ii])
            ii = ii + 1
            if ii == n_cards or prihdr.ascardlist().keys()[ii] == card_hist :
                Block.add_history(hist_entry, details)

    def read(self, scans=None, bands=None, force_tuple=False, IFs=None) :
        """Read in data from the fits file.

        This method reads data from the fits file including the files history
        and basically every peice of data that could be needed.  It is done,
        one scan and one IF at a time.  Each scan and IF is returned in an
        instance of the DataBlock class (defined in another module of this
        package).

        Parameters
        ----------
        scans : tuple of integers
            Which scans in the file to be processed.  A list of 
            integers, with 0 corresponding to the lowest numbered scan.
            Default is all of them.
        bands : tuple of integers
            Which intermediate frequencies (also called frequency windows)
            to process.  A list of integers with 0 coorsponding to the 
            lowest frequency present. Default is all of them.
            TODO : Overlapping frequency windows stiched together somehow.
        force_tuple: By default, if there is only a single output Data
            Block, it is returned not wraped in a tuple, but if we want to
            loop over the output we can force the output to be a tuple,
            even if it only has one element.
        IFs : tuple of integers
            Depricated, use `bands`.

        Returns
        -------
        Instance of the DataBlock class, or a tuple of these instances
        (if asked for multiple scans and IFs).
        """
        
        # `bands` and `IFs` is are two names for the same parameter.
        if not bands is None and IFs is None:
            IFs = bands
        # We want scans and IFs to be a sequence of indicies.
        if scans is None :
            scans = range(len(self.scan_set))
        elif not hasattr(scans, '__iter__') :
            scans = (scans, )
        elif len(scans) == 0 :
            scans = range(len(self.scan_set))
        if IFs is None :
            IFs = range(len(self.IF_set))
        elif not hasattr(IFs, '__iter__') :
            IFs = (IFs, )
        elif len(IFs) == 0 :
            IFs = range(len(self.IF_set))
        
        # Sequence of output DataBlock objects.
        output = ()
        for scan_ind in scans :
            for IF_ind in IFs :
                # Choose the appropriate records from the file, get that data.
                inds_sif = self.get_scan_IF_inds(scan_ind, IF_ind)
                Data_sif = db.DataBlock(self.fitsdata.field('DATA')[inds_sif])
                # Masked data is stored in FITS files as float('nan')
                Data_sif.data[sp.logical_not(sp.isfinite(
                                   Data_sif.data))] = ma.masked
                # Now iterate over the fields and add them
                # to the data block.
                for field, axis in fields_and_axes.iteritems() :
                    # See if this fits file has the key we are looking for.
                    try:
                        names = self.fitsdata.names
                    except AttributeError:
                        names = self.fitsdata._names
                    if not field in names :
                        continue
                    # First get the 'FITS' format string.
                    field_format = self.hdulist[1].columns.formats[
                                self.hdulist[1].columns.names.index(field)]
                    if axis :
                        # From the indices in inds_sif, we only need a
                        # subset: which_data will subscript inds_sif.
                        temp_data = self.fitsdata.field(field)[inds_sif]
                        # For reshaping at the end.
                        field_shape = []
                        for ii, single_axis in enumerate(Data_sif.axes[0:-1]) :
                            # For each axis, slice out all the data except the
                            # stuff we need.
                            which_data = [slice(None)] * 3
                            if single_axis in axis :
                                field_shape.append(Data_sif.dims[ii])
                            else :
                                which_data[ii] = [0]
                            temp_data = temp_data[tuple(which_data)]
                        temp_data.shape = tuple(field_shape)
                        Data_sif.set_field(field, temp_data, axis, field_format)
                    else :
                        Data_sif.set_field(field, self.fitsdata.field(field)
                            [inds_sif[0,0,0]], axis, field_format)
                if hasattr(self, 'history') :
                    Data_sif.history = db.History(self.history)
                else :
                    self.history =bf.get_history_header(self.hdulist[0].header)
                    #self.set_history(Data_sif)
                    fname_abbr = ku.abbreviate_file_path(self.fname)
                    self.history.add('Read from file.', ('File name: ' + 
                                         fname_abbr, ))
                    Data_sif.history = db.History(self.history)
                Data_sif.verify()
                output = output + (Data_sif, )
        if self.feedback > 2 :
            print 'Read finished.'
        if len(output) == 1 and not force_tuple :
            return output[0]
        else :
            return output

    def __del__(self) :
        self.hdulist.close()
        if self.feedback > 3 :
            print "Closed ", self.fname


class Writer() :
    """Class that writes data back to fits files.

    This class acculumates data stored in DataBlock objects using the
    'add_data(DataBlock)' method.  Once the user has added all the data, 
    she can then call the 'write(file_name)' method to write it to file.
    """
    
    def __init__(self, Blocks=None, feedback=2) :
        """Init script for the fitsGBT Writer.

        The writer can optionally be initialized with a sequence of DataBlock
        objects to be written to file, but these can also be added later with
        the add data method.  It also can take an optional feedback parameter
        (integer default 2) which specifies the amount of junk to print.
        """
        
        self.feedback = feedback

        self.first_block_added = True
        self.field = {}
        self.formats = {}
        if not Blocks is None:
            self.add_data(Blocks)

    def add_data(self, Blocks) :
        """Interface for adding DataBlock objects to the Writter.
        
        This method can be passed either a single DataBlock object or any
        sequence of DataBlock objects.  They will all be added to the Writer's
        data which can eventually be written as a fits file.
        """

        if not hasattr(Blocks, '__iter__') :
            self._add_single_block(Blocks)
        else :
            for Block in Blocks :
                self._add_single_block(Block)

    def _add_single_block(self, Block) :
        """Adds all the data in a DataBlock Object to the Writer such that it
        can be written to a fits file eventually."""
        
        Block.verify()
        # Merge the histories
        if self.first_block_added :
            self.history = db.History(Block.history)
        else :
            self.history = db.merge_histories(self.history, Block)
        # Some dimensioning and such
        dims = tuple(Block.dims)
        n_records = dims[0]*dims[1]*dims[2]
        block_shape = dims[0:-1]
        # For now automatically determine the format for the data field.
        data_format = str(dims[-1]) + 'E'
        if self.first_block_added :
            self.data_format = data_format
        elif self.data_format != data_format :
            raise ce.DataError('Data shape miss match: freq axis must be same'
                               ' length for all DataBlocks added to Wirter.')

        # Copy the reshaped data from the DataBlock
        data = sp.array(ma.filled(Block.data, float('nan')))
        if self.first_block_added :
            self.data = data.reshape((n_records, dims[3]))
        else :
            self.data = sp.concatenate((self.data, data.reshape((
                                        n_records, dims[3]))), axis=0)

        # Now get all stored fields for writing out.
        for field, axes in Block.field_axes.iteritems() :
            # Need to expand the field data to the full ntimes x npol x ncal
            # length (with lots of repitition).  We will use np broadcasting.
            broadcast_shape = [1,1,1]
            for axis in axes :
                axis_ind = list(Block.axes).index(axis)
                broadcast_shape[axis_ind] = dims[axis_ind]
            # Allowcate memory for the new full field.
            data_type = Block.field[field].dtype
            field_data = sp.empty(block_shape, dtype=data_type)
            # Copy data with the entries, expanding dummy axes.
            field_data[:,:,:] = sp.reshape(Block.field[field],
                                                 broadcast_shape)
            if self.first_block_added :
                self.field[field] = field_data.reshape(n_records)
                self.formats[field] = Block.field_formats[field]
            else :
                self.field[field] = sp.concatenate((self.field[field],
                                        field_data.reshape(n_records)), axis=0)
                if self.formats[field] != Block.field_formats[field] :
                    raise ce.DataError('Format miss match in added data blocks'
                                       ' and field: ' + field)
        self.first_block_added = False

    def write(self, file_name) :
        """Write stored data to file.
        
        Take all the data stored in the Writer (from added DataBlocks) and
        write it to a fits file with the passed file name.
        """

        # Add the data
        Col = pyfits.Column(name='DATA', format=self.data_format, 
                            array=self.data)
        columns = [Col,]
        
        # Add all the other stored fields.
        for field_name in self.field.iterkeys() :
            Col = pyfits.Column(name=field_name,
                                format=self.formats[field_name],
                                array=self.field[field_name])
            columns.append(Col)
        coldefs = pyfits.ColDefs(columns)
        # Creat fits header data units, one for the table and the mandatory
        # primary.
        tbhdu = pyfits.new_table(coldefs)
        prihdu = pyfits.PrimaryHDU()
        # Add the write history.
        fname_abbr = ku.abbreviate_file_path(file_name)
        self.history.add('Written to file.', ('File name: ' + fname_abbr,))
        # Add the history to the header.
        bf.write_history_header(prihdu.header, self.history)

        # Combine the HDUs and write to file.
        hdulist = pyfits.HDUList([prihdu, tbhdu])
        hdulist.writeto(file_name, clobber=True)
        if self.feedback > 0 :
            print 'Wrote data to file: ' + fname_abbr

