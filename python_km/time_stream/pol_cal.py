"""Module that performs polarization calibration.
"""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single
import map.tools
from core import fits_map

# base_single.BaseSingle is a base class that knows how to read an input file
# loop over files, scans and IFs.  It does all the input and output.  The only
# thing it doesn't know how to do is the science, which is what will be added
# here.
class Calibrate(base_single.BaseSingle) :
    """Pipeline module that corrects for polarization leakage.
    """
    
    # Here we define a bunch of stuff that BaseSingle needs to know to do its
    # thing.
    # prefix is a few letters that are added to all parameter names that are
    # read from a configureation file.
    prefix = 'pc_'
    # These are the parameters that should be read from file.  These are in
    # addition to the ones defined at the top of base_single.py.
    params_init = {
                   'mueler_file' : 'default_fname'
                   }

    # The base single initialization method does a bunch of stuff, but we want
    # to add one thing.  We want to read a mueler matrix from file.
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
        base_single.BaseSingle.__init__(self, parameter_file_or_dict,
                                        feedback)
        # Read in the mueler matrix file.
        mueler_file_name = self.params['mueler_file']
        self.mueler = read_mueler(mueler_file_name)
    
    # This function tells BaseSingle what science to do.  Data is a
    # core.data_block.DataBlock object.  It holds all the data for a single
    # scan and a single IF.  BaseSingle knows how to loop over all of these.
    # More on DataBlock objects in the calibrate function below.
    def action(self, Data) :
        calibrate_pol(Data, self.mueler)
        Data.add_history('Corrected for polarization leakage.', 
                         ('Mueler matrix file: ' + self.params['mueler_file'],))

        return Data


# This fuction will need to be written.
def read_mueler(file_name) :
    """Reads a Mueler matrix from file."""
    return None


def calibrate_pol(Data, mueler) :
    """Subtracts a Map out of Data."""
        
    # Data is a DataBlock object.  It holds everything you need to know about
    # the data in a single scan and IF.  You should get to know them very well.
    # Data.data is a numpy masked array (see numpy documentation) and holds the
    # acctual data.  It is a 4 dimensional array.  The demensions are (in
    # order): (time, pol, cal, freq).  Each dimension can be any length which
    # you can figure out by looking at Data.dims = sp.shape(Data.data).
    # Data.field is a python dictionary that holds all the other data that you
    # might care about from the origional fits file.  For instance,
    # Data.field['CAL'] is an array with length dims[2].  It normally has
    # values ['T', 'F']. Data.field['CRVAL4'] tells you about the polarization
    # axis of Data.data.  By SDfits convension each polarization is represented
    # by an integer: 1=I, 2=Q, 3=U, 4=V, -5=XX, -6=YY, -7=XY, -8=YX.

    # Also this depends on having the polarizations rotated correctly to IQUV.
    # Some code to do this has been hacked together in the rotate_pol module,
    # but I don't trust it yet.

    # Some dimension checks.
    # We expect 4 polarizations.
    if not Data.dims[1] == 4 :
        raise ce.DataError('Require 4 polarizations.')
    # We expect polarizations to be in order IQUV.
    if (Data.field['CRVAL4'][0] != 1 or Data.field['CRVAL4'][1] != 2 or
        Data.field['CRVAL4'][2] != 3 or Data.field['CRVAL4'][3] != 4) :
        raise ce.DataError('Expected the polarization basis to be IQUV.')

    # A useful function that might need:
    Data.calc_freq()
    # Now data has an atribute Data.freq which is an array that gives the
    # frequency along the last axis.


