"""First crack at a map maker.

This map maker assumes uncorrelated data, i.e. it is the same as griding the
data, averageing over pointings that land in the same grid.
"""

import scipy as sp
import numpy.ma as ma

from kiyopy import parse_ini
import kiyopy.custom_exceptions as ce
from core import utils, data_block, fitsGBT

params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ["testfile_GBTfits"],
               'input_end' : ".fits",
               'output_root' : "./",
               'output_end' : ".map.fits",
               # What data to process within each file.
               'scans' : [],
               'IF' : 0,
               # Map parameters (Ra (deg), Dec (deg), f (MHz)).
               'resample_freqs' : False,
               # centre[2], shape[2] and freq_spacing are ignored unless
               # resample_freqs is True.
               'centre' : (315.0, 0.0, 800.0),
               'shape' : (40, 40, 200),
               'spacing' : 0.125, # degrees
               'freq_spacing' : 1.0 # MHz
               }


def mk_map_simple(parameter_file_or_dict=None) :
    """Main function for this simple map maker."""
    
    # Read in the parameters.
    params = parse_ini.parse(parameter_file_or_dict, params_init)
    shape = params['shape']

    if not type(params['IF']) is int :
        raise ce.FileParameterTypeError('Can only process a single IF.')

    # Flag for the first block processed (will allowcate memory on the first
    # iteration).
    first_block = True

    # Generate bins for ra, dec
    ra_bins = calc_bins(params['centre'][0], 
                        params['spacing']/sp.cos(params['centre'][1]),
                        shape[0], 'middle')
    dec_bins = calc_bins(params['centre'][1], params['spacing'],
                        shape[1], 'middle')

    # Loop over the files to process.
    for file_middle in params['file_middles'] :
        input_fname = (params['input_root'] + file_middle +
                       params['input_end'])
        output_fname = (params['output_root']
                        + file_middle + params['output_end'])
        Writer = fitsGBT.Writer()
        
        # Read in the data, and loop over data blocks.
        Reader = fitsGBT.Reader(input_fname)
        Blocks = Reader.read(params['scans'], params['IF'])
        for Data in Blocks :
            dims = Data.dims
            Data.calc_freq()
            if first_block :
                if not params['resample_freqs'] :
                    shape[2] = dims[-1]
                    freq_bins = Data.freq/1.0e6
                else :
                    raise NotImplementedError('Frequency resampling.')
                    freq_bins = calc_bins(params['centre'][2],
                                 params['freq_spacing'], shape[2], 'lower')

                # Allowcate memory for the map and pointing counts.
                map = ma.zeros(shape, dtype=float)
                counts = sp.zeros(shape, dtype=int)

            # Figure out the pointing pixel index and the frequency indicies.
            Data.calc_pointing()
            ra_inds = calc_inds(Data.ra, params['centre'][0], shape[0],
                                params['spacing']/sp.cos(params['centre'][1]))
            dec_inds = calc_inds(Data.dec, params['centre'][1], shape[0],
                                params['spacing'])
            data = Data.data
            
            ## This code resamples the data along the frequency axis.
            #freq_inds = calc_inds(Data.freq/1.0e6, params['centre'][2], shape[2],
            #                      params['freq_spacing'])
            #data_resampled = ma.empty((dims[0], dims[1], dims[2], shape[2]))
            #for new_f_ind in range(shape[2]) :
            #    inds, = sp.where(sp.equal(freq_inds, new_f_ind))
            #    if len(inds) == 0 :
            #        data_resampled[:,:,:,new_f_ind] = ma.masked
            #    else :
            #        # Looks like ma.medain has a bug and the following line
            #        # fails.  Will have to loop instead.
            #        #data_resampled[:,:,:,new_f_ind] = ma.median(Data.data[
            #        #                              :,:,:,freq_inds[inds]], 3)
            #        for ii in range(dims[0]) :
            #            for jj in range(dims[1]) :
            #                for kk in range(dims[2]) :
            #                    data_resampled[ii,jj,kk,new_f_ind] = (
            #                        ma.median(Data.data[
            #                        ii,jj,kk,freq_inds[inds]]))

            for time_ind in range(dims[0]) :
                # Take cal on and cal off, XX and YY to be independant
                # measurements (only really valid for calibrated data).
                for jj in range(dims[2]) :
                    inds = sp.where(data_resampled.get_mask()[time_ind,0,jj,:])
                    map[ra_inds[time_ind], dec_inds[time_ind], inds] += (
                        data_resampled[time_ind,0,jj,inds])
                    counts[ra_inds[time_ind], dec_inds[time_ind], inds] += 1
                    inds = sp.where(data_resampled.get_mask()[time_ind,3,jj,:])
                    map[ra_inds[time_ind], dec_inds[time_ind], inds] += (
                        data_resampled[time_ind,3,jj,inds])
                    counts[ra_inds[time_ind], dec_inds[time_ind], inds] += 1
                        

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

    inds = sp.array((pointing - centre) / spacing + shape/2.0, dtype=int)
    
    return inds

def calc_bins(centre, shape, spacing=1., edge='lower') :
    """Calculates bin edges.

    Calculates a 1D array of regularly spaced bins.  User specifies whethar
    lower, middle, or upper side bins are to be calculated.
    """
    
    if edge == 'lower' :
        offset = 0.
    elif edge == 'middle' :
        offset = spacing/2.0
    elif edge == 'upper' :
        offset = float(spacing)
    else :
        raise ValueError("Argument edge must be either 'lower', 'middle' or "
                         "'upper'")
    
    bins = (spacing * sp.arange(-shape/2.0, shape/2.0, dtype=float)
            + centre + offset)
    
    return bins


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        mk_map_simple(str(sys.argv[1]))
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        mk_map_simple()
                


