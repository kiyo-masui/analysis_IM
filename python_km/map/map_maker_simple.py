"""First crack at a map maker.

This map maker assumes uncorrelated data, i.e. it is the same as griding the
data, averageing over pointings that land in the same grid.
"""

import scipy as sp
import numpy.ma as ma

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from core import utils, data_block, fitsGBT, data_map, fits_map

# Parameters prefixed with 'mm_' when read from file.
params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./testoutput",
               'output_end' : ".map.fits",
               # What data to process within each file.
               'scans' : (),
               'IFs' : (0,),
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               'map_shape' : (40, 40),
               'pixel_spacing' : 0.125, # degrees
               # What times streams to include in map.
               'cal_weights' : (1.0, 1.0),
               'pol_weights' : (1.0, 0.0, 0.0, 1.0)
               }

# Fields that should be copied from times stream data.
fields_to_copy = (
                  'BANDWID',
                  'OBJECT'
                  )

class MapMaker(object) :
    """Simple data gridder."""
   
    def __init__(self, parameter_file_or_dict=None) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix='mm_')

    def execute(self) :
        params = self.params
        last_slash = params['output_root'].rfind('/')
        if last_slash > 0 : # Protect against last_slash == -1 or 0.
            outdir = params['output_root'][:last_slash]
            kiyopy.utils.mkdir_p(outdir)
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix='mm_')
        # Rename some commonly used parameters.
        shape = params['map_shape']
        spacing = params['pixel_spacing']
        ra_spacing = -spacing/sp.cos(params['field_centre'][1]*2.0*sp.pi/180.)

        if len(params['IFs']) != 1 :
            raise ce.FileParameterTypeError('Can only process a single IF.')

        # Flag for the first block processed (will allowcate memory on the 
        # first iteration).
        first_block = True

        # Generate bins for ra, dec
        ra_bins = calc_bins(params['field_centre'][0], shape[0], 
                            ra_spacing, 'middle')
        dec_bins = calc_bins(params['field_centre'][1], shape[1], 
                             params['pixel_spacing'], 'middle')

        # Loop over the files to process.
        for file_middle in params['file_middles'] :
            input_fname = (params['input_root'] + file_middle +
                           params['input_end'])
            output_fname = (params['output_root']
                            + file_middle + params['output_end'])
            
            # Read in the data, and loop over data blocks.
            Reader = fitsGBT.Reader(input_fname)
            Blocks = Reader.read(params['scans'], params['IFs'])
            for Data in Blocks :
                dims = Data.dims
                Data.calc_freq()
                if first_block :
                    shape = shape + (dims[-1],)
                    freq_bins = Data.freq/1.0e6

                    # Allowcate memory for the map and pointing counts.
                    counts = sp.zeros(shape, dtype=int)
                    Map = data_map.DataMap(ma.zeros(shape, dtype=float))
                    # Copy some data over that all the data_blocks should have
                    # in common.
                    # Do not copy RA and DEC to the map, since they cannot be
                    # stored in the fits file.
                    #Map.set_field('RA', ra_bins, ('long',), '1E') 
                    #Map.set_field('DEC', dec_bins, ('lat',), '1E')
                    Map.history = Data.history
                    for key in fields_to_copy :
                        Map.set_field(key, Data.field[key], 
                           Data.field_axes[key], Data.field_formats[key])
                    Map.set_field('CTYPE3', 'FREQ--HZ', (), '32A')
                    Map.set_field('CTYPE1', 'RA---DEG', (), '32A')
                    Map.set_field('CTYPE2', 'DEC--DEG', (), '32A')
                    # Copy frequency axis (now the third axis not the first).
                    Map.set_field('CRVAL3', Data.field['CRVAL1'], (), 'D')
                    Map.set_field('CRPIX3', Data.field['CRPIX1'], (), 'D')
                    Map.set_field('CDELT3', Data.field['CDELT1'], (), 'D')
                    # Set the other two axes.
                    Map.set_field('CRPIX1', shape[0]//2, (), 'D')
                    Map.set_field('CRVAL1', ra_bins[shape[0]//2], (), 'D')
                    Map.set_field('CDELT1', ra_spacing, (), 'D')
                    Map.set_field('CRPIX2', shape[1]//2, (), 'D')
                    Map.set_field('CRVAL2', dec_bins[shape[1]//2], (), 'D')
                    Map.set_field('CDELT2', spacing, (), 'D')

                else :
                    Map.history = data_map.merge_histories(Map, Data)

                # Figure out the pointing pixel index and the frequency 
                # indicies.
                Data.calc_pointing()
                ra_inds = calc_inds(Data.ra, params['field_centre'][0], 
                                    shape[0], ra_spacing)
                dec_inds = calc_inds(Data.dec, params['field_centre'][1], 
                                     shape[1], params['pixel_spacing'])
                # Create a linear combination of the cals and pols that we 
                # want to map.
                data = ma.zeros((dims[0],dims[3]))
                for ii_pol in range(dims[1]) :
                    for jj_cal in range(dims[2]) :
                        data += (Data.data[:,ii_pol,jj_cal,:]
                                 * params['pol_weights'][ii_pol]
                                 * params['cal_weights'][jj_cal])

                add_data_2_map(data, ra_inds, dec_inds, Map.data, counts)
                first_block = False
        uncovered_inds = sp.where(counts==0)
        Map.data[uncovered_inds] = ma.masked
        Map.data /= counts
        Map.add_history('Gridded data with map_maker_simple.')
        
        Map.verify()
        fits_map.write(Map, params['output_root'] + params['output_end'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini')





def add_data_2_map(data, ra_inds, dec_inds, map, counts) :
    """Add a data masked array to a map."""

    ntime = len(ra_inds)
    shape = sp.shape(map)
    if len(dec_inds) != ntime or len(data[:,0]) != ntime :
        raise ValueError('Time axis of data, ra_inds and dec_inds must be'
                         ' same length.')
    if sp.shape(counts) != sp.shape(map) or len(map[0,0,:]) != len(data[0,:]) :
        raise ValueError('Map-counts shape mismatch or data frequency axis '
                         'length mismatch.')
    for time_ind in range(ntime) :
        if (ra_inds[time_ind] >= 0 and ra_inds[time_ind] < shape[0] and
            dec_inds[time_ind] >= 0 and dec_inds[time_ind] < shape[1]) :
            # Get unmasked
            unmasked_inds = sp.logical_not(ma.getmaskarray(data[time_ind,:]))
            map[ra_inds[time_ind], dec_inds[time_ind], unmasked_inds] += \
                      data[time_ind, unmasked_inds]
            counts[ra_inds[time_ind], dec_inds[time_ind], unmasked_inds] += 1


def calc_inds(pointing, centre, shape, spacing=1.) :
    """Calculates a 1D map index corresponding to a location.

    Given a position, this funciton returns the 1D index on map of that
    position.  The map is specified by its centre, number of pixels across and
    the pixel width.

    Function is vectorized for input pointing.
    
    This fuction might be sufficiently general that it could be moved to
    core.utils.  We will see.
    """
    
    if type(pointing) is list or type(pointing) is tuple :
        pointing = sp.array(pointing, dtype=float)
    shape = int(shape)

    inds = sp.array((pointing - centre) / spacing + shape/2.0, dtype=int)
    
    return inds

def calc_bins(centre, shape, spacing=1., edge='lower') :
    """Calculates bin edges.

    Calculates a 1D array of regularly spaced bins.  User specifies whethar
    lower, middle, or upper side bins are to be calculated.
    """
    
    if edge == 'left' :
        offset = 0.
    elif edge == 'middle' :
        offset = spacing/2.0
    elif edge == 'right' :
        offset = float(spacing)
    else :
        raise ValueError("Argument edge must be either 'left', 'middle' or "
                         "'right'")
    shape = int(shape)
    
    bins = (spacing * sp.arange(-shape/2.0, shape/2.0, dtype=float)
            + centre + offset)
    
    return bins


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        MapMaker(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        MapMaker().execute()
                


