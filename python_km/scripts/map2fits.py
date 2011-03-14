"""Converts a map stored as an algebra vect object (.npy) into a fits file."""

import sys

import scipy as sp

from core import algebra, data_map, fits_map, hist

def convert(map_file, history_file=None) :
    """Main function."""
    
    map = algebra.load(map_file)
    map = algebra.make_vect(map)

    if map.axes != ('freq', 'ra', 'dec') :
        raise NotImplementedError("Exepected input map to be organized "
                                  "('freq', 'ra', 'dec').")
    
    new_shape = map.shape[1:] + (map.shape[0],)
    
    # Make the out file name assuming the input file end in .npy.  This is a
    # hack and someone should fix it sometime.
    out_fname = map_file.split('/')[-1][:-4] + '.fits'

    Map_fits = data_map.DataMap(sp.rollaxis(map, 0, 3))
    # Set axis names
    Map_fits.set_field('CTYPE3', 'FREQ--HZ', (), '32A')
    Map_fits.set_field('CTYPE1', 'RA---DEG', (), '32A')
    Map_fits.set_field('CTYPE2', 'DEC--DEG', (), '32A')
    # Copy frequency axis (now the third axis not the first).
    Map_fits.set_field('CRVAL3', map.info['freq_centre'], (), 'D')
    Map_fits.set_field('CRPIX3', new_shape[2]//2+1, (), 'D')
    Map_fits.set_field('CDELT3', map.info['freq_delta'], (), 'D')
    # Set the other two axes.
    Map_fits.set_field('CRVAL1', map.info['ra_centre'], (), 'D')
    Map_fits.set_field('CRPIX1', new_shape[0]//2+1, (), 'D')
    Map_fits.set_field('CDELT1', map.info['ra_delta'], (), 'D')
    Map_fits.set_field('CRVAL2', map.info['dec_centre'], (), 'D')
    Map_fits.set_field('CRPIX2', new_shape[1]//2+1, (), 'D')
    Map_fits.set_field('CDELT2', map.info['dec_delta'], (), 'D')
    
    # Copy the file history if provided.
    if not history_file is None :
        history = hist.read(history_file)
        history.add("Converted map to fits.", ("File name: " + out_fname,))
        Map_fits.history = history
    
    # Verify contents and write out.
    Map_fits.verify()
    fits_map.write(Map_fits, out_fname)


if __name__ == '__main__' :
    if len(sys.argv) == 3 :
        convert(str(sys.argv[1]), sys.argv[2])
    elif len(sys.argv) == 2 :
        convert(str(sys.argv[1]))
    else :
        print ("Usage: 'map2fits.py mapfile.npy' or 'map2fits.py mapfile.npy "
               "historyfile.hist'")

