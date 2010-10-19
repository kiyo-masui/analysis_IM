#!/usr/bin/python
"""This module performs hanning smoothing.
"""

import numpy.ma as ma

import base_single
import kiyopy.custom_exceptions as ce


class Hanning(base_single.BaseSingle) :
    """Pipeline module that performs hanning smmothing on the data."""

    prefix = 'ha_'

    def action(self, Data):
        do_hanning(Data)
        Data.add_history('Hanning smoothed.')
        return Data

def do_hanning(Data) :
    """Perform Hanning smoothing.

    This function accepts a DataBlock class and returns nothing.  It changes
    the data attribute of the passed Data.
    """
    
    # Perform Hanning smoothing frequency space convolution.
    # Masked array automatically flags data adjacent to masked
    # data.
    if len(Data.axes) != 4 or Data.axes[-1] != 'freq' :
        raise ce.DataError('Hanning smoothing expects data be 4D with frequency '
                        'being the last axis.')
    # Channels adjacent to masked data automatically masked.
    Data.data[:,:,:,1:-1] = ( 0.25*Data.data[:,:,:,:-2] + 
                         0.50*Data.data[:,:,:,1:-1] + 
                         0.25*Data.data[:,:,:,2:] )
    # End points where not smoothed, mask them as bad data.
    Data.data[:,:,:,0] = ma.masked
    Data.data[:,:,:,-1] = ma.masked


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Hanning(str(sys.argv[1])).execute()

