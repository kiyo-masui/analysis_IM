"""Classes and utilities to read and write maps to an from fits.

This will eventually be as comprehensive as fitsGBT.py, but for now there are
just a few helper functions.
"""

import scipy as sp
import numpy.ma as ma
import pyfits

import data_map
import base_fits as bf
import base_data
import kiyopy.custom_exceptions as ce
import kiyopy.utils as ku

card_hist = 'DB-HIST'
card_detail = 'DB-DET'

# List of fields that should be read and written to fits files.  Since this is
# an image, they are always scalars, never arrays.  They are stored in the
# header, not with the image.
fields = (
          'POL',
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

def write(Maps, file_name, feedback=2) :
    """Write a map to fits file.

    Map should be a map_data.MapData object.
    """

    # If only a single Map was passed, make it iterable.
    if not hasattr(Maps, '__iter__') :
        Maps = (Maps,)
    # First create a primary with the history and such:
    prihdu = pyfits.PrimaryHDU()
    # Add history to the primary.
    fname_abbr = ku.abbreviate_file_path(file_name)
    # Add final history entry and store in the primary header.
    history = base_data.merge_histories(*Maps)
    history.add('Written to file.', 'File name: ' + fname_abbr)
    bf.write_history_header(prihdu.header, history)
    list_of_hdus = [prihdu]
    
    for ii, Map in enumerate(Maps) :
        # Creat an image HDU.
        map = sp.array(ma.filled(Map.data, float('nan')))
        imhdu = pyfits.ImageHDU(sp.swapaxes(map, 0, 2), name='MAP%d'%ii)

        # Add extra data to the HDU
        for key in Map.field.iterkeys() :
            if Map.field_axes[key] != () :
                raise ce.DataError('Only 0D data can be written to a Fits Map '
                                   'Header.')
            card = pyfits.Card(key, Map.field[key].item())
            imhdu.header.ascardlist().append(card)
        list_of_hdus.append(imhdu)

    # Creat the HDU list and write to file.
    hdulist = pyfits.HDUList(list_of_hdus)
    hdulist.writeto(file_name, clobber=True)
    if feedback > 0 :
        print 'Wrote data to file: ' + fname_abbr


def read(file_name, map_inds=None, feedback=2) :
    """Read a map from an image fits file.
    """

    fname_abbr = ku.abbreviate_file_path(file_name)
    if feedback > 0 :
        print 'Opening file: ' + fname_abbr
    # Open the fits file.
    hdulist = pyfits.open(file_name)
    history = bf.get_history_header(hdulist[0].header)
    history.add('Read from file.', 'File name: ' + fname_abbr)
    map_list = []

    if map_inds is None :
        map_inds = range(1, len(hdulist))
    elif not hasattr(map_inds, '__iter__') :
        map_inds = (map_inds, )
    elif len(scans) == 0 :
        map_inds = range(1, len(hdulist))
    
    for ii in map_inds :
        data = hdulist[ii].data

        # Set the data attriute
        Map = data_map.DataMap()
        Map.set_data(sp.swapaxes(data,0,2))
        # Masked data is stored in FITS files as float('nan')
        Map.data[sp.logical_not(sp.isfinite(Map.data))] = ma.masked
        Map.history = history

        # Set the other fields.
        for field_name in fields :
            if not field_name in hdulist[1].header.keys() :
                continue
            value = hdulist[1].header[field_name]
            Map.set_field(field_name, value)
        map_list.append(Map)
    if len(map_list) == 1:
        map_list = map_list[0]

    return map_list


