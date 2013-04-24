"""Make a dirty map (and noise matrix) from data.

Loops over polarizations and only consideres only the 0th cal state (if you
want something else, change it in the time stream).

TODO: finish this?
"""

import copy

import scipy as sp
import numpy.ma as ma
import numpy.linalg as linalg

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from core import utils, data_block, fitsGBT, algebra, hist
import tools

def build_hitmap(file_middles, input_end, output_root, scans, IFs,
                 field_centre, map_shape, pixel_spacing):
    # Rename some commonly used parameters.
    ra_spacing = -spacing/sp.cos(field_centre[1]*sp.pi/180.)
    if len(IFs) != 1:
        raise ce.FileParameterTypeError('Can only process a single IF.')

    all_file_names = []

    # Flag for the first block processed (will allowcate memory on the
    # first iteration).
    first_block = True
    # Loop over the files to process.
    try:
        for file_middle in file_middles:
            input_fname = (input_root + file_middle +
                            input_end)
            # Read in the data, and loop over data blocks.
            Reader = fitsGBT.Reader(input_fname, feedback=feedback)
            Blocks = Reader.read(scans, IFs)

            # Calculate the time varience at each frequency.  This will
            # be used as weights in most algorithms.
            for Data in Blocks:
                dims = Data.dims
                # On first pass set up the map parameters.
                if first_block:
                    shape = map_shape + (dims[-1],)
                    Data.calc_freq()
                    centre_freq = Data.freq[dims[-1]//2]
                    delta_freq = Data.field['CDELT1']
                    # Allocate memory for the map.
                    map_data = sp.zeros(shape, dtype=float)
                    map_data = algebra.make_vect(map_data,
                                    axis_names=('ra', 'dec', 'freq'))
                    first_block=False
                # Figure out the pointing pixel index and the frequency
                # indicies.
                Data.calc_pointing()
                ra_inds = tools.calc_inds(Data.ra,
                            field_centre[0], shape[0],
                            ra_spacing)
                dec_inds = tools.calc_inds(Data.dec,
                                        field_centre[1],
                                        shape[1], pixel_spacing)
                data = Data.data[:,pol_ind,0,:]
                pixel_hits[:] = 0
                pixel_list = pixel_counts(data, ra_inds, dec_inds,
                                            pixel_hits, map_shape=shape[0:2])
                # End Blocks for loop.

            # Free up all that memory and flush memory maps to file.
            del map_data

            # Save the file names for the history.
            all_file_names.append(kiyopy.utils.abbreviate_file_path(
                map_file_name))
    except ce.NextIteration:
        pass
        # End polarization for loop.
    history.add("Made dirty map.", all_file_names)
    h_file_name = (output_root + 'history.hist')
    history.write(h_file_name)

def add_data_2_map(data, ra_inds, dec_inds, map, noise_i=None, weight=1):
    """Add a data masked array to a map.

    This function also adds the weight to the noise matrix for diagonal noise.
    """

    ntime = len(ra_inds)
    shape = sp.shape(map)
    if len(dec_inds) != ntime or len(data[:,0]) != ntime:
        raise ValueError('Time axis of data, ra_inds and dec_inds must be'
                         ' same length.')
    if not noise_i is None and map.shape != noise_i.shape:
        raise ValueError('Inverse noise array must be the same size as the map'
                         ' or None.')

    for time_ind in range(ntime):
        if (ra_inds[time_ind] >= 0 and ra_inds[time_ind] < shape[0] and
            dec_inds[time_ind] >= 0 and dec_inds[time_ind] < shape[1]):
            # Get unmasked
            unmasked_inds = sp.logical_not(ma.getmaskarray(data[time_ind,:]))
            ind_map = (ra_inds[time_ind], dec_inds[time_ind], unmasked_inds)
            map[ind_map] += (weight*data)[time_ind, unmasked_inds]
            if not noise_i is None:
                if not hasattr(weight, '__iter__'):
                    noise_i[ind_map] += weight
                else:
                    noise_i[ind_map] += weight[unmasked_inds]

def pixel_counts(data, ra_inds, dec_inds, pixel_hits, map_shape=(-1,-1)):
    """Counts the hits on each unique pixel.

    Returns pix_list, a list of tuples, each tuple is a (ra,dec) index on a
    map pixel hit on this scan.  The list only contains unique entries.  The
    array pixel_hits (preallocated for performance), is
    filled with the number of hits on each of these pixels as a function of
    frequency index. Only the entries pixel_hits[:len(pix_list),:]
    are meaningful.
    """

    if ra_inds.shape != dec_inds.shape or ra_inds.ndim != 1:
        raise ValueError('Ra and Dec arrays not properly shaped.')
    if (pixel_hits.shape[-1] != data.shape[-1] or pixel_hits.shape[0] <
        len(ra_inds)):
        raise ValueError('counts not allowcated to right shape.')

    pix_list = []
    for ii in range(len(ra_inds)):
        pix = (ra_inds[ii], dec_inds[ii])
        if ((map_shape[0] > -1 and pix[0] >= map_shape[0]) or
            (map_shape[1] > -1 and pix[1] >= map_shape[1]) or
            pix[0] < 0 or pix[1] < 0):
            continue
        elif not pix in pix_list:
            pix_list.append(pix)
        unmasked_freqs = sp.logical_not(ma.getmaskarray(data)[ii,:])
        pixel_hits[pix_list.index(pix), unmasked_freqs] += 1

    return pix_list


# For running this module from the command line
if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        DirtyMap(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else:
        DirtyMap().execute()

