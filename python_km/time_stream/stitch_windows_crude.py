#! /usr/bin/python
"""This module stitches the frequency windows together.
"""

import copy

import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single

class Stitch(base_single.BaseSingle) :
    """Pipeline module that flags rfi and other forms of bad data.

    For lots of information look at the doc-string for flag_data.apply_cuts.
    """

    prefix = 'sw_'
    params_init = {
                   }
    
    def scan_action(self, scan_blocks) :
        return scan_blocks

def stitch(blocks) :
    """Stitches to gether data blocks from different frequency windows.

    Accepts a tuple of Data Block objects to be stitched together.
    """
    
    # Stitched data starts life as a copy of one of the old data blocks.
    NewData = copy.deepcopy(blocks[0])
    
    for Data in blocks :
        # Make sure the data we need is here.
        if not (Data.field.has_key('CRVAL1') 
                and Data.field.has_key('CDELT1')
                and Data.field.has_key('CRPIX1')) :
            raise ce.DataError('All blocks must have frequency axis'
                               ' information.')
        # Make sure all the data not involved in the stitching is the same.
        for key, value in NewData.field.iteritems() :
            pass



# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Stitch(str(sys.argv[1])).execute()

