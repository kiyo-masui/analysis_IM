"""Make a dirty map (and noise matrix) from data.

Loops over polarizations and only consideres only the 0th cal state (if you
want something else, change it in the time stream).
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

prefix = 'dm_'
# Parameters prefixed with 'dm_' when read from file.
params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./testoutput_",
               # What data to process within each file.
               'scans' : (),
               'IFs' : (0,),
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               # In pixels.
               'map_shape' : (5, 5),
               'pixel_spacing' : 0.5, # degrees
               # What polarizations to process.
               'polarizations' : (),
               # What noise model to use.  Options are:
               #  'grid' : just grid the data and devid by number of counts.
               #  'diag_file' : measure the time varience of each file.
               #  'disjoint_scans' : treat each scan as highly correlated
               #                     internally but uncorrelated between scans.
               # 'diag*' and 'grid' produce a noise file in the same format as
               # the map. 'disjoint_scans' produces a full noise covarience
               # (n_pix^2 numbers), each polarization in its onwn file.
               'noise_model' : 'grid',
               # Measured noise parameters.
               'noise_parameters_input_root' : 'None'
               }

class DirtyMapMaker(object) :
    """Converts time stream data into a dirty map.
    
    This module reads in multiple time stream files and turns them into a 
    (dirty) map. There are several options for how the noise is treated.
    'grid' and 'diag_file' both assume diagonal noise, while 'disjoint_scans'
    assumes that data entries can only be compared within a scan.
    
    This last one uses a full N^2 covarience matrix which requires a rather
    large amount of memory. It is not reccomended that you run this on one
    anything but the smallest map sizes, unless you have several Gigs of RAM.
    This is meant to be run on machines with > 24 GB of RAM.

    The net result is noise wieghted map, aka the dirty map.  The equivalend
    algebraic operation is. P^T N^(-1) d.  The map inverse noise matrix is 
    also produced: C_n^(-1) = P^T N^(-1) P.
    """
   
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback

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
        map_shape = params['map_shape']
        spacing = params['pixel_spacing']
        algorithm = params['noise_model']
        noise_root = params['noise_parameters_input_root']
        ra_spacing = -spacing/sp.cos(params['field_centre'][1]*sp.pi/180.)
        if not algorithm in ('grid', 'diag_file', 'disjoint_scans') :
            raise ValueError('Invalid noise model: ' + algorithm)
        if len(params['IFs']) != 1 :
            raise ce.FileParameterTypeError('Can only process a single IF.')

        # Set up to iterate over the pol states.
        npol = 2 # This will be reset when we read the first data block.
        pol_ind = 0

        all_file_names = []

        while pol_ind < npol :
            # Flag for the first block processed (will allowcate memory on the 
            # first iteration).
            first_block = True
            # Loop over the files to process.
            try :
                for file_middle in params['file_middles'] :
                    input_fname = (params['input_root'] + file_middle +
                                   params['input_end'])
                    # Read in the data, and loop over data blocks.
                    Reader = fitsGBT.Reader(input_fname, feedback=self.feedback)
                    Blocks = Reader.read(params['scans'], params['IFs'])
                    
                    # Calculate the time varience at each frequency.  This will
                    # be used as weights in most algorithms.
                    if not algorithm == 'grid' :
                        if not noise_root == 'None':
                            # We have measured variance.
                            noise_pars = sp.load(noise_root + file_middle 
                                                 + ".npy")
                            var = noise_pars[params['IFs'][0], pol_ind, 0, :]
                        else :
                            # We need to measure the variance.
                            var = tools.calc_time_var_file(Blocks, pol_ind, 0)
                            # Convert from masked array to array.
                            var = var.filled(9999.)
                    else :
                        var = 1.
                    weight = 1/var
                    
                    for Data in Blocks :
                        dims = Data.dims
                        # On first pass set up the map parameters.
                        if first_block :
                            shape = map_shape + (dims[-1],)
                            Data.calc_freq()
                            centre_freq = Data.freq[dims[-1]//2]
                            delta_freq = Data.field['CDELT1']
                            if pol_ind==0 :
                                # Figure out the length of the polarization
                                # loop.
                                npol = dims[1]
                                # Accumulate the data history.
                                history = hist.History(Data.history)
                            # Get the current polarization integer.
                            this_pol = Data.field['CRVAL4'][pol_ind]
                            # Check that we even want to make a dirty map for
                            # this polarization.
                            if ((not utils.polint2str(this_pol) in
                                params['polarizations']) and
                                params['polarizations']) :
                                # Break to the end of the polarization loop.
                                raise ce.NextIteration()
                            # Allowcate memory for the map.
                            map_data = sp.zeros(shape, dtype=float)
                            map_data = algebra.make_vect(map_data,
                                            axis_names=('ra', 'dec', 'freq'))
                            # Allowcate memory for the inverse map noise.
                            if algorithm in ('grid', 'diag_file') :
                                noise_inv = sp.zeros(shape, dtype=float)
                                noise_inv = algebra.make_mat(noise_inv,
                                            axis_names=('ra', 'dec', 'freq'),
                                            row_axes=(0,1,2), col_axes=(0,1,2))
                            elif algorithm in ('disjoint_scans', 'ds_grad') :
                                # At each frequency use full N^2 noise matrix, 
                                # but assume each frequency has uncorrelated
                                # noise. This is a big matrix so make sure it
                                # is reasonable.
                                size = shape[0]^2*shape[1]^2*shape[2]
                                if size > 4e9 : # 16 GB
                                    raise RunTimeError('Map size too big. '
                                                       'Asked for a lot '
                                                       'of memory.')
                                noise_inv = sp.zeros(shape[0:2] + shape,
                                                     dtype=sp.float32)
                                noise_inv = algebra.make_mat(noise_inv,
                                                axis_names=('ra', 'dec', 'ra',
                                                            'dec', 'freq'),
                                                row_axes=(0,1,4), 
                                                col_axes=(2,3,4))
                                # Allowcate memory for temporary data. Hold the
                                # number of times each pixel in this scan is
                                # hit. Factor of 2 longer in time in case some
                                # scans are longer than first block (guppi).
                                pixel_hits = sp.empty((2*dims[0], dims[-1]))
                            first_block = False
                        else :
                            if pol_ind==0 :
                                history.merge(Data)
                        # Figure out the pointing pixel index and the frequency 
                        # indicies.
                        Data.calc_pointing()
                        ra_inds = tools.calc_inds(Data.ra, 
                                    params['field_centre'][0], shape[0], 
                                    ra_spacing)
                        dec_inds = tools.calc_inds(Data.dec, 
                                             params['field_centre'][1], 
                                             shape[1], params['pixel_spacing'])
                        data = Data.data[:,pol_ind,0,:]
                        if algorithm in ('grid', 'diag_file') :
                            add_data_2_map(data, ra_inds, dec_inds, map_data, 
                                           noise_inv, weight)
                        elif algorithm in ('disjoint_scans', ) :
                            add_data_2_map(data - ma.mean(data, 0), ra_inds, 
                                           dec_inds, map_data, None, weight)
                            pixel_hits[:] = 0
                            pixel_list = pixel_counts(data, ra_inds, dec_inds, 
                                            pixel_hits, map_shape=shape[0:2])
                            add_scan_noise(pixel_list, pixel_hits, var,
                                           noise_inv)
                        # End Blocks for loop.
                    # End file name for loop.
                # Now write the dirty maps out for this polarization. 
                # Use memmaps for this since we want to reorganize data 
                # and write at the same time.
                # New maps will have the frequency axis as slowly varying, for
                # future efficiency.
                map_file_name = (params['output_root'] + 'dirty_map_' +
                                 utils.polint2str(this_pol) + '.npy')
                mfile = algebra.open_memmap(map_file_name, mode='w+',
                                             shape=(shape[2],) + shape[:2])
                map_mem = algebra.make_vect(mfile, axis_names=('freq', 'ra',
                                                                'dec'))
                # And the noise matrix.
                noise_file_name = (params['output_root'] + 'noise_inv_' +
                                 utils.polint2str(this_pol) + '.npy')
                if algorithm in ('disjoint_scans', 'ds_grad') :
                    mfile = algebra.open_memmap(noise_file_name, mode='w+',
                                                shape=(shape[2],) +
                                                shape[:2]*2)
                    noise_mem = algebra.make_mat(mfile, axis_names=('freq', 
                                'ra', 'dec', 'ra', 'dec'), row_axes=(0, 1, 2), 
                                col_axes=(0, 3, 4))
                else :
                    mfile = algebra.open_memmap(noise_file_name, mode='w+',
                                                 shape=(shape[2],) +
                                                shape[:2])
                    noise_mem = algebra.make_mat(mfile, axis_names=('freq', 
                        'ra', 'dec'), row_axes=(0, 1, 2), col_axes=(0, 1, 2))
                # Give the data arrays axis information.
                map_mem.set_axis_info('freq', centre_freq, delta_freq)
                map_mem.set_axis_info('ra', params['field_centre'][0], 
                                      ra_spacing)
                map_mem.set_axis_info('dec', params['field_centre'][1], 
                                      params['pixel_spacing'])
                noise_mem.set_axis_info('freq', centre_freq, delta_freq)
                noise_mem.set_axis_info('ra', params['field_centre'][0], 
                                        ra_spacing)
                noise_mem.set_axis_info('dec', params['field_centre'][1], 
                                      params['pixel_spacing'])
                # Copy the data to the memory maps after rearranging.  
                # The roll_axis should return a view, so this should
                # be memory efficient.
                map_mem[...] = sp.rollaxis(map_data, -1)
                noise_mem[...] = sp.rollaxis(noise_inv, -1)

                # Free up all that memory and flush memory maps to file.
                del mfile, map_mem, noise_mem, map_data, noise_inv 
                
                # Save the file names for the history.
                all_file_names.append(kiyopy.utils.abbreviate_file_path(
                    map_file_name))
                all_file_names.append(kiyopy.utils.abbreviate_file_path(
                    noise_file_name))
            except ce.NextIteration :
                pass
            pol_ind += 1
            # End polarization for loop.
        history.add("Made dirty map.", all_file_names)
        h_file_name = (params['output_root'] + 'history.hist')
        history.write(h_file_name)

def add_data_2_map(data, ra_inds, dec_inds, map, noise_i=None, weight=1) :
    """Add a data masked array to a map.
    
    This function also adds the weight to the noise matrix for diagonal noise.
    """

    ntime = len(ra_inds)
    shape = sp.shape(map)
    if len(dec_inds) != ntime or len(data[:,0]) != ntime :
        raise ValueError('Time axis of data, ra_inds and dec_inds must be'
                         ' same length.')
    if not noise_i is None and map.shape != noise_i.shape :
        raise ValueError('Inverse noise array must be the same size as the map'
                         ' or None.')

    for time_ind in range(ntime) :
        if (ra_inds[time_ind] >= 0 and ra_inds[time_ind] < shape[0] and
            dec_inds[time_ind] >= 0 and dec_inds[time_ind] < shape[1]) :
            # Get unmasked
            unmasked_inds = sp.logical_not(ma.getmaskarray(data[time_ind,:]))
            ind_map = (ra_inds[time_ind], dec_inds[time_ind], unmasked_inds)
            map[ind_map] += (weight*data)[time_ind, unmasked_inds]
            if not noise_i is None :
                if not hasattr(weight, '__iter__') :
                    noise_i[ind_map] += weight
                else :
                    noise_i[ind_map] += weight[unmasked_inds]

def pixel_counts(data, ra_inds, dec_inds, pixel_hits, map_shape=(-1,-1)) :
    """Counts the hits on each unique pixel.

    Returns pix_list, a list of tuples, each tuple is a (ra,dec) index on a 
    map pixel hit on this scan.  The list only contains unique entries.  The
    array pixel_hits (preallowcated for performance), is
    filled with the number of hits on each of these pixels as a function of
    frequency index. Only the entries pixel_hits[:len(pix_list), :] 
    are meaningful.
    """

    if ra_inds.shape != dec_inds.shape or ra_inds.ndim != 1 :
        raise ValueError('Ra and Dec arrays not properly shaped.')
    if (pixel_hits.shape[-1] != data.shape[-1] or pixel_hits.shape[0] < 
        len(ra_inds)) :
        raise ValueError('counts not allowcated to right shape.')

    pix_list = []
    for ii in range(len(ra_inds)) :
        pix = (ra_inds[ii], dec_inds[ii])
        if ((map_shape[0] > -1 and pix[0] >= map_shape[0]) or 
            (map_shape[1] > -1 and pix[1] >= map_shape[1]) or
            pix[0] < 0 or pix[1] < 0) :
            continue
        elif not pix in pix_list :
            pix_list.append(pix)
        unmasked_freqs = sp.logical_not(ma.getmaskarray(data)[ii,:])
        pixel_hits[pix_list.index(pix), unmasked_freqs] += 1

    return pix_list

def add_scan_noise(pixels, pixel_hits, variance, noise_inv) :
    """Adds a scan to the map inverse noise matrix.

    Passed a list of map indices that were hit this scan and an array of the
    number of times each pixel was hit at each frequency.  This is all added to
    noise_inv, inversly wieghted by variance.  Governing formula found in
    Kiyo's notes November 27, 2010.
    """
    
    # Calculate number of good pointings are each frequency.
    n_pix = len(pixels)
    pointings = sp.sum(pixel_hits[:n_pix,:], 0)
    f_inds = pointings > 0
    point = pointings[f_inds]
    if hasattr(variance, '__iter__') :
        var = variance[f_inds]
    else :
        var = variance
    for ii, p1 in enumerate(pixels) :
        noise_inv[p1 + p1 + (f_inds,)] += sp.array(pixel_hits[ii,f_inds], 
                                                   dtype=float)/var
        for jj, p2 in enumerate(pixels) :
            noise_inv[p1 + p2 + (f_inds,)] -= (sp.array(pixel_hits[ii,f_inds], 
                                   dtype=float)*pixel_hits[jj,f_inds]/point/var)


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        DirtyMap(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        DirtyMap().execute()

