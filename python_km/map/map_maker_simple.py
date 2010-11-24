"""First crack at a map maker.

This map maker assumes uncorrelated data, i.e. it is the same as griding the
data, averageing over pointings that land in the same grid.
"""

import copy

import scipy as sp
import numpy.ma as ma

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from core import utils, data_block, fitsGBT, data_map, fits_map
import tools

# Parameters prefixed with 'mm_' when read from file.
params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./testoutput",
               'output_end' : ".fits",
               # What data to process within each file.
               'scans' : (),
               'IFs' : (0,),
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               'map_shape' : (40, 40),
               'pixel_spacing' : 0.125, # degrees
               # What time streams to include in map. Should sum to 1.
               'cal_weights' : (0.5, 0.5),
               'pol_weights' : (0.5, 0.0, 0.0, 0.5),
               # Weight each file by its varience.
               'variance_weight' : False
               }

class MapMaker(object) :
    """Simple data gridder."""
   
    def __init__(self, parameter_file_or_dict=None) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix='mm_')

    def execute(self, nprocesses=1) :
        """Function that acctually does the work.

        The nprocesses parameter does not do anything yet.  It is just there
        for compatibility with the pipeline manager.
        """
        params = self.params
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix='mm_')
        # Rename some commonly used parameters.
        shape = params['map_shape']
        spacing = params['pixel_spacing']
        ra_spacing = -spacing/sp.cos(params['field_centre'][1]*sp.pi/180.)

        if len(params['IFs']) != 1 :
            raise ce.FileParameterTypeError('Can only process a single IF.')

        # Flag for the first block processed (will allowcate memory on the 
        # first iteration).
        first_block = True

        # Generate bins for ra, dec
        ra_bins = tools.calc_bins(params['field_centre'][0], shape[0], 
                            ra_spacing, 'middle')
        dec_bins = tools.calc_bins(params['field_centre'][1], shape[1], 
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
            if params['variance_weight'] :
                weight = 0.0
                file_counts = 0.0
                for Data in Blocks :
                    for ii_pol in range(len(params['pol_weights'])) :
                        for jj_cal in range(len(params['cal_weights'])) :
                            weight += (ma.sum(Data.data[:,ii_pol,jj_cal,:]**2,
                                              0) * params['pol_weights'][ii_pol]
                                       * params['cal_weights'][jj_cal])
                            file_counts += (Data.dims[0] - ma.count_masked(
                                            Data.data[:,ii_pol,jj_cal,:], 0)
                                            * params['pol_weights'][ii_pol]
                                            * params['cal_weights'][jj_cal])
                weight = file_counts/weight.filled(99999.)
            else :
                weight = 1.0

            for Data in Blocks :
                dims = Data.dims
                Data.calc_freq()
                if first_block :
                    freq_bins = Data.freq/1.0e6
                    shape = shape + (dims[-1],)

                    # Allowcate memory for the map and pointing counts.
                    counts = ma.zeros(shape, dtype=float)
                    Map = tools.set_up_map(Data, params['field_centre'], 
                                           shape[0:2], (ra_spacing, spacing))
                else :
                    Map.history = data_map.merge_histories(Map, Data)

                # Figure out the pointing pixel index and the frequency 
                # indicies.
                Data.calc_pointing()
                ra_inds = tools.calc_inds(Data.ra, params['field_centre'][0], 
                                    shape[0], ra_spacing)
                dec_inds = tools.calc_inds(Data.dec, params['field_centre'][1], 
                                     shape[1], params['pixel_spacing'])
                # Create a linear combination of the cals and pols that we 
                # want to map.
                data = ma.zeros((dims[0],dims[3]))
                for ii_pol in range(dims[1]) :
                    for jj_cal in range(dims[2]) :
                        data += (Data.data[:,ii_pol,jj_cal,:]
                                 * params['pol_weights'][ii_pol]
                                 * params['cal_weights'][jj_cal])

                add_data_2_map(data, ra_inds, dec_inds, Map.data,
                               counts, weight)
                first_block = False
        uncovered_inds = sp.where(counts==0)
        Map.data[uncovered_inds] = ma.masked
        counts[uncovered_inds] = ma.masked
        Map.data /= counts
        noise = 1/counts
        Noise = copy.deepcopy(Map)
        Map.add_history('Gridded data with map_maker_simple.')
        Noise.data = noise
        
        Map.verify()
        Noise.verify()
        fits_map.write((Map, Noise), 
                       params['output_root'] + params['output_end'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini')


def add_data_2_map(data, ra_inds, dec_inds, map, counts, weight=1) :
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
                      (weight*data)[time_ind, unmasked_inds]
            if not hasattr(weight, '__iter__') :
                counts[ra_inds[time_ind], dec_inds[time_ind], unmasked_inds] +=\
                          weight
            else :
                counts[ra_inds[time_ind], dec_inds[time_ind], unmasked_inds] +=\
                          weight[unmasked_inds]


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        MapMaker(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        MapMaker().execute()
                


