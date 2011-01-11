"""This module contains the class that holds an IF and scan of GBT data"""

import scipy as sp
import matplotlib as plt

import utils
import base_data


class DataBlock(base_data.BaseData) :
    """Class that holds an single IF and scan of GBT data.

    This is the main vessel for storing an transporting GBT data.  This class
    can be used with the fitsGBT.py module to be read and written as a properly
    formatted fits file.  The rax data is accessed and updated through the
    'data' and 'field' attributes of this class and associated hleper
    functions.

    Please remember that when working with the 'data' attribute, that it is a
    numpy MaskedArray class, not a normal numpy array.  Take care to use the
    masked versions of any numpy functions to preserve the mask.  This is
    especially usefull for flagging bad data and RFI.

    Most of the functionality here is inherited from
    data_container.DataContainer.
    """
    
    # These are the valid axes that a data field can vary over.  Any other
    # field can vary over only the first three of these.
    axes = ('time', 'pol', 'cal', 'freq')

    # The following methods calculate useful quantities, but assume certain
    # fields exist.  They should be valid if the DataBlock was read from a GBT
    # fits file.
    
    def calc_pointing(self) :
        """Calculates the telescope pointing.

        At every time the Ra and Dec of the telescope time is calculated.
        These are stored as attributes (not fields) named ra and dec.  This
        requires the fields 'CRVAL3', 'CRVAL2' and 'DATE-OBS' to be set.
        """
        self.ra = sp.zeros(self.dims[0])
        self.dec = sp.zeros(self.dims[0])
        for ii in range(self.dims[0]) :
            self.ra[ii], self.dec[ii] = utils.elaz2radecGBT(
                                            self.field['CRVAL3'][ii],
                                            self.field['CRVAL2'][ii],
                                            self.field['DATE-OBS'][ii])

    def calc_freq(self) :
        """Calculates the frequency axis.
        
        The frequency axis is stored as a attribute (not a field) named freq.
        This requires the fileds 'CRVAL1', 'CPIX1' and 'CDELT1' to be set.
        """
        self.freq = ((sp.arange(self.dims[-1], dtype=float) + 1.0 - 
                     self.field['CRPIX1'])*self.field['CDELT1'] + 
                     self.field['CRVAL1'])

    def plot_spectra(self, times=(), pols=(), cals=(), time_average=False) :
        """Make a plot along frequency axis."""
        #

    def calc_time(self) :
        self.time = sp.empty(self.dims[0])
        for ii in range(self.dims[0]) :
            self.time[ii] = utils.time2float(self.field['DATE-OBS'][ii])


# Clone some extra functions:

merge_histories = base_data.merge_histories
History = base_data.History
