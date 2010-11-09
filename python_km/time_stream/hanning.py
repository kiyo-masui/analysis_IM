#!/usr/bin/python
"""This module performs hanning smoothing.
"""

import numpy.ma as ma
from scipy import weave
from scipy.weave import converters

import base_single
import kiyopy.custom_exceptions as ce


class Hanning(base_single.BaseSingle) :
    """Pipeline module that performs hanning smmothing on the data."""

    prefix = 'ha_'

    def action(self, Data) :
        hanning_smooth(Data)
        Data.add_history('Hanning smoothed.')
        return Data

def hanning_smooth(Data) :
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
    
    # Below copiled c code replaces the above numpy code.  Speed up of ~20.
    # Not working.
    #data = Data.data.data
    #mask = Data.data.mask
    #nt, npol, ncal, nf = data.shape

    #code = """
    #       #line 43 "hanning.py"
    #       double this_data, last_data;
    #       int this_mask, last_mask, ind;
    #
    #       for (int tii=0; tii<nt; tii++) {
    #           ind = tii*npol*ncal*nf;
    #           for (int pjj=0; pjj<npol; pjj++) {
    #               ind = ind + pjj*ncal*nf;
    #               for (int ckk=0; ckk<ncal; ckk++) {
    #                   ind = ind + ckk*nf;
    #                   last_data = data[ind];
    #                   last_mask = mask[ind];
    #                   for (int fmm=1; fmm<nf-1; fmm++) {
    #                       this_data = data[ind+fmm];
    #                       this_mask = mask[ind+fmm];
    #                       if (last_mask || mask[ind+fmm+1])
    #                           mask[ind+fmm] = 1;
    #                       else if (!this_mask)
    #                           data[ind+fmm] = 0.25*last_data + 0.5*this_data 
    #                                           + 0.25*data[ind+fmm+1];
    #                       last_data = this_data;
    #                       last_mask = this_mask;
    #                       }
    #                   }
    #               }
    #           }
    #       """
    #weave.inline(code, ['data', 'mask', 'nt', 'npol', 'ncal', 'nf'])
    
    # End points where not smoothed, mask them as bad data.
    Data.data[:,:,:,0] = ma.masked
    Data.data[:,:,:,-1] = ma.masked


# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Hanning(str(sys.argv[1])).execute()

