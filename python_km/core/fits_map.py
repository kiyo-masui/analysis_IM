"""Classes and utilities to read and write maps to an from fits.

This will eventually be as comprehensive as fitsGBT.py, but for now there are
just a few helper functions.
"""

import scipy as sp
import numpy.ma as ma
import pyfits

card_hist = 'DB-HIST'
card_detail = 'DB-DET'

def write(Map, file_name, feedback=2) :
    """Write a map to fits file.

    Map should be a map_data.MapData object.
    """

    # First create a primary with the history and such:
    prihdu = pyfits.PrimaryHDU()
    # Add history to the primary.
    history_keys  = Map.history.keys()
    history_keys.sort()
    for hist in history_keys :
        details = Map.history[hist]
        # Chop off the number, since they are already sorted.
        hcard = pyfits.Card(card_hist, hist[5:])
        prihdu.header.ascardlist().append(hcard)
        for detail in details :
            dcard = pyfits.Card(card_detail, detail)
            prihdu.header.ascardlist().append(dcard)
    hcard = pyfits.Card(card_hist, 'Written to file.')
    prihdu.header.ascardlist().append(hcard)
    dcard = pyfits.Card(card_detail, 'File name: ' + file_name)
    prihdu.header.ascardlist().append(dcard)
    
    # Creat an image HDU.
    map = sp.array(ma.filled(Map.data, float('nan')))
    imhdu = pyfits.ImageHDU(sp.swapaxes(map, 0, 2), name='MAP')

    # Creat the HDU list and write to file.
    hdulist = pyfits.HDUList([prihdu, imhdu])
    hdulist.writeto(file_name, clobber=True)
    if feedback > 0 :
        print 'Wrote data to file: ' + file_name


