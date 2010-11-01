"""Classes and utilities to read and write maps to an from fits.

This will eventually be as comprehensive as fitsGBT.py, but for now there are
just a few helper functions.
"""

import scipy as sp
import numpy.ma as ma
import pyfits

import kiyopy.custom_exceptions as ce
import kiyopy.utils as ku

card_hist = 'DB-HIST'
card_detail = 'DB-DET'

# List of fields that should be read an written to fits files.  Since this is
# an image, they are always scalars, never arrays.  The are stored int he
# header, not with the image.
fields = (
          'BANDWID',
          'OBJECT',
          'CTYPE1', # Axis type (Ra)  
          'CRVAL1', # Ra axis centre
          'CRPIX1', # Centre pixel
          'CDELT1', # Pixel width
          'CTYPE2', # Axis type (Dec)  
          'CRVAL2',
          'CRPIX2',
          'CDELT2',
          'CTYPE3', # Axis type (Freq)  
          'CRVAL3', 
          'CRPIX3',
          'CDELT3'
          )

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
    dcard = pyfits.Card(card_detail, 'File name: ' + 
                        ku.abbreviate_file_path(file_name))
    prihdu.header.ascardlist().append(dcard)
    
    # Creat an image HDU.
    map = sp.array(ma.filled(Map.data, float('nan')))
    imhdu = pyfits.ImageHDU(sp.swapaxes(map, 0, 2), name='MAP')

    # Add extra data to the HDU
    for key in fields :
        if Map.field_axes[key] != () :
            raise ce.DataError('Only 0D data can be written to a Fits Map '
                               'Header.')
        card = pyfits.Card(key, Map.field[key].item())
        imhdu.header.ascardlist().append(card)

    # Creat the HDU list and write to file.
    hdulist = pyfits.HDUList([prihdu, imhdu])
    hdulist.writeto(file_name, clobber=True)
    if feedback > 0 :
        print 'Wrote data to file: ' + ku.abbreviate_file_path(file_name)


