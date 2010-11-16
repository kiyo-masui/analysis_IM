#! /usr/bin/python
"""This module stitches the frequency windows together.
"""

import copy

import numpy.ma as ma
import scipy as sp

import kiyopy.custom_exceptions as ce
import base_single

class Stitch(base_single.BaseSingle) :
    """Pipeline module stitches frequency windows together.
    """

    prefix = 'sw_'
    params_init = {
                   }
    
    def scan_action(self, scan_blocks) :
        Stiched = stitch(scan_blocks)
        Stiched.add_history('Stiched Frequency windows together.')
        return Stiched

def stitch(blocks) :
    """Stitches to gether data blocks from different frequency windows.

    Accepts a tuple of Data Block objects to be stitched together.
    """
    
    # Sort data blocks by frequeny.
    try :
        blocks = sorted(blocks, key=lambda Data : -1*Data.field['CRVAL1'])
    except KeyError :
        raise ce.DataError('All blocks must have frequency axis'
                           ' information.')
    # Stitched data starts life as a copy of one of the old data blocks.
    OutData = copy.deepcopy(blocks[0])
    
    # First make sure all the data is compatible.
    for Data in blocks :
        # Make sure the data we need is here.
        if not (Data.field.has_key('CRVAL1') 
                and Data.field.has_key('CDELT1')
                and Data.field.has_key('CRPIX1')) :
            raise ce.DataError('All blocks must have frequency axis'
                               ' information.')
        # Make sure all the data not involved in the stitching is the same.
        # For now enforce that CDELT1 be the same for all IFs.
        for key, field_data in OutData.field.iteritems() :
            if not key in ('CRVAL1', 'CRPIX1', 'OBSFREQ', 'RESTFREQ') :
                if not Data.field.has_key(key) :
                    raise ce.DataError('All blocks must have the same data '
                                       'fields.')
                # Treat strings differently.
                if OutData.field_formats[key][-1] != 'A' :
                    if not sp.allclose(field_data, Data.field[key]) :
                        raise ce.DataError('All blocks to be stitched must '
                                           'have matching data fields execpt '
                                           'for frequency axis information.')
                else :
                    if not sp.alltrue(field_data == Data.field[key]) :
                        raise ce.DataError('All blocks to be stitched must '
                                           'have matching data fields execpt '
                                           'for frequency axis information.')
    # For now assume that the frequencies are reversed ordered.
    if OutData.field['CDELT1'] >= 0 :
        raise NotImplementedError('Expected frequency steps to be negitive.')
    
    delt = abs(OutData.field['CDELT1']) 
    # Loop over data and stitch.
    for Data in blocks[1:] :
        # Get the freqeuncy axes
        OutData.calc_freq()
        Data.calc_freq()

        n_over = list(Data.freq >= OutData.freq[-1]).count(True)
        if n_over == 0 :
            raise ce.DataError('Frequency windows do not overlap.')
        # Use mean, not sum, to normalize in case of flagged data.
        factor = ma.mean(OutData.data[:,:,:,-2*n_over//3:-n_over//3], axis=3)
        factor /= ma.mean(Data.data[:,:,:,n_over//3:2*n_over//3], axis=3)
        Data.data *= factor[:,:,:,sp.newaxis]
        OutData.set_data(ma.concatenate((OutData.data[:,:,:,:-n_over//2], 
                                         Data.data[:,:,:,n_over//2:]), axis=3))
    OutData.calc_freq()
    OutData.set_field('BANDWIDTH', abs(OutData.freq[0] - OutData.freq[-1]),
                      (), format='D')

    return OutData



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Stitch(str(sys.argv[1])).execute()

