"""A map maker.

Loops over polarizations and only consideres the 0th cal (if you want something
else, change it in the time stream).
"""

import copy

import scipy as sp
import numpy.ma as ma
import numpy.linalg as linalg

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
               'output_root' : "./testoutput_",
               'output_end' : ".fits",
               # What data to process within each file.
               'scans' : (),
               'IFs' : (0,),
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               'map_shape' : (5, 5),
               'pixel_spacing' : 0.5, # degrees
               # What time streams to include in map. Should sum to 1.
               #'cal_weights' : (0.5, 0.5),
               #'pol_weights' : (0.5, 0.0, 0.0, 0.5),
               # What noise model to use.  Options are:
               #  'grid' : just grid the data and devid by number of counts.
               #  'diag_file' : measure the time varience of each file.
               #  'disjoint_scans' : treat each scan as highly correlated
               #                     internally but uncorrelated between scans.
               # 'diag*' and 'grid' produce a noise file in the same format as
               # the map. 'disjoint_scans' produces a full noise covarience
               # (n_pix^2 numbers), each polarization in its onwn file.
               'noise_model' : 'grid'
               }

class MapMaker(object) :
    """Converts time stream data into a Map.
    
    This module reads in multiple time stream files and turns them into a map.
    There are several options for how the noise is treated.  'grid' and
    'diag_file' both assume diagonal noise, while 'disjoint_scans' assumes that
    data entries can only be compared within a scan.
    
    This last one uses a full N^2 covarience matrix which requires a rather
    large amount of memory. It is not reccomended that you run this on one
    anything but the smallest map sizes, unless you have several Gigs of RAM.
    This is meant to be run on machines with > 24 GB of RAM.
    """
   
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix='mm_')
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
        ra_spacing = -spacing/sp.cos(params['field_centre'][1]*sp.pi/180.)
        if not algorithm in ('grid', 'diag_file', 'disjoint_scans') :
            raise ValueError('Invalid noise model: ' + algorithm)
        if len(params['IFs']) != 1 :
            raise ce.FileParameterTypeError('Can only process a single IF.')
        
        # Set up to interate over pols.
        first_pol = True
        map_list = []
        if algorithm in ('grid', 'diag_file', 'disjoint_scans') :
            noise_list = []
        npol = 2 # This will be reset when we read the first data block.
        pol_ind = 0

        while pol_ind < npol :
            # Flag for the first block processed (will allowcate memory on the 
            # first iteration).
            first_block = True
            # Loop over the files to process.
            for file_middle in params['file_middles'] :
                input_fname = (params['input_root'] + file_middle +
                               params['input_end'])
                # Read in the data, and loop over data blocks.
                Reader = fitsGBT.Reader(input_fname)
                Blocks = Reader.read(params['scans'], params['IFs'])
                
                # Calculate the time varience at each frequency.  This will be
                # used as weights in most algorithms.
                if not algorithm == 'grid' :
                    var = tools.calc_time_var_file(Blocks, 0, 0).filled(999.)
                else :
                    var = 1.
                weight = 1/var

                for Data in Blocks :
                    dims = Data.dims
                    Data.calc_freq()
                    # On first pass set up the map parameters.
                    if first_block :
                        shape = map_shape + (dims[-1],)
                        if first_pol :
                            npol = dims[1]
                        Map = tools.set_up_map(Data, params['field_centre'], 
                                            shape[0:2], (ra_spacing, spacing))
                        # Also store the polarization.
                        Map.set_field('POL',
                                      Data.field['CRVAL4'][pol_ind], (), 'I')
                        # Will store the map data outside of Map for now.
                        map_data = sp.zeros(shape, dtype=float)
                        # How many pointings each pixel gets.
                        counts = sp.zeros(shape, dtype=float)
                        # Allowcate memory for the inverse map noise.
                        if algorithm in ('grid', 'diag_file') :
                            noise_inv = sp.zeros(shape, dtype=float)
                        elif algorithm in ('disjoint_scans', 'ds_grad') :
                            # At each frequency use full N^2 noise matrix, but
                            # assume each frequency has uncorrelated noise.
                            # This is a big matrix so make sure it is
                            # reasonable.
                            size = shape[0]^2*shape[1]^2*shape[2]
                            if size > 4e9 : # 16 GB
                                raise RuntimError('Map size too big.  Asked '
                                                  'for a lot of memory.')
                            noise_inv = sp.zeros(shape[0:2] + shape,
                                                 dtype=sp.float32)
                            # Allowcate memory for temporary data. Hold the
                            # number of times each pixel in this scan is hit.
                            pixel_hits = sp.empty((dims[0], dims[-1]))
                        first_block = False
                    else :
                        Map.history = data_map.merge_histories(Map, Data)

                    # Figure out the pointing pixel index and the frequency 
                    # indicies.
                    Data.calc_pointing()
                    ra_inds = tools.calc_inds(Data.ra, 
                               params['field_centre'][0], shape[0], ra_spacing)
                    dec_inds = tools.calc_inds(Data.dec, 
                                         params['field_centre'][1], 
                                         shape[1], params['pixel_spacing'])
                    data = Data.data[:,pol_ind,0,:]
                    if algorithm in ('grid', 'diag_file') :
                        add_data_2_map(data, ra_inds, dec_inds, map_data, 
                                       noise_inv, counts, weight)
                    elif algorithm in ('disjoint_scans', ) :
                        add_data_2_map(data - ma.mean(data, 0), ra_inds, 
                                       dec_inds, map_data, None, counts, weight)
                        pixel_hits[:] = 0
                        pixel_list = pixel_counts(data, ra_inds, dec_inds, 
                                        pixel_hits, map_shape=shape[0:2])
                        add_scan_noise(pixel_list, pixel_hits, var, noise_inv)
            
            if algorithm in ('grid', 'diag_file') :
                uncovered_inds = counts==0
                noise_inv[uncovered_inds] = 1.0
                map_data = (map_data/noise_inv)
                noise = 1/noise_inv
            elif algorithm in ('disjoint_scans',) :
                nf = shape[-1]
                if (not sp.all(sp.isfinite(noise_inv)) or 
                    not sp.all(sp.isfinite(map_data))) :
                    raise RuntimeError("Didn't do a good enough job of "
                                       "avioding devide by 0.")
                if self.feedback > 1 :
                    print 'Solving for optimal map.'
                    print 'Discarding near zero eigenvalues at each frequency:'
                for f_ind in xrange(nf) :
                    inds_i, inds_j = sp.where(counts[:,:,f_ind]!=0)
                    map_vect = map_data[inds_i, inds_j, f_ind]
                    noise_inv_mat = noise_inv[inds_i, inds_j, :, :, f_ind]
                    noise_inv_mat = noise_inv_mat[:, inds_i, inds_j]
                    if noise_inv_mat.size == 0 :
                        continue
                    # This is where the magic happens.  Solve for the optimal
                    # map from the noise inverse and the dirty map.
                    L, S = linalg.eigh(noise_inv_mat)
                    cutoff = 1.0e-6
                    singular = (L < cutoff*max(L))
                    discarded = float(len(sp.nonzero(singular)))/len(singular)
                    if self.feedback > 1 :
                        print "%d: %6.2f"%(f_ind, discarded*100)+"% discarded"
                    L[singular] = 1.
                    L_inv = 1/L
                    L_inv[singular] = 0.
                    map_rot = sp.dot(S.T, map_vect)
                    map_rot[singular] = 0.
                    map_vect = sp.dot(S, L_inv*map_rot)
                    map_vect -= sp.mean(map_vect)
                    map_data[inds_i, inds_j, f_ind] = map_vect
                uncovered_inds = counts==0
                diag_noise_inv = noise_inv.diagonal(0,1,3).diagonal(0,0,1)
                diag_noise_inv = diag_noise_inv.swapaxes(0,2)
                diag_noise_inv[uncovered_inds] = 1.0
                noise = 1/diag_noise_inv

            Map.set_data(map_data)
            Map.data[uncovered_inds] = ma.masked
            Map.add_history('Gridded data with map_maker_simple.', 
                            ('Algorithm: ' + algorithm,))
            Noise = copy.deepcopy(Map)
            Noise.set_data(noise)
            Noise.data[uncovered_inds] = ma.masked
            
            Noise.verify()
            noise_list.append(Noise)
                
            Map.verify()
            map_list.append(Map)
            pol_ind += 1
            first_pol = False

        fits_map.write(map_list, 
                       params['output_root'] + 'map' + params['output_end'])
        if algorithm in ('grid', 'diag_file', 'disjoint_scans') :
            fits_map.write(noise_list, params['output_root'] + 'noise' 
                           + params['output_end'])

def add_data_2_map(data, ra_inds, dec_inds, map, noise_i, counts, weight=1) :
    """Add a data masked array to a map."""

    ntime = len(ra_inds)
    shape = sp.shape(map)
    if len(dec_inds) != ntime or len(data[:,0]) != ntime :
        raise ValueError('Time axis of data, ra_inds and dec_inds must be'
                         ' same length.') 
    if sp.shape(counts) != sp.shape(map) or len(map[0,0,:]) != len(data[0,:]) :
        raise ValueError('Map-counts shape  mismatch or data frequency axis '
                         'length mismatch.')

    for time_ind in range(ntime) :
        if (ra_inds[time_ind] >= 0 and ra_inds[time_ind] < shape[0] and
            dec_inds[time_ind] >= 0 and dec_inds[time_ind] < shape[1]) :
            # Get unmasked
            unmasked_inds = sp.logical_not(ma.getmaskarray(data[time_ind,:]))
            ind_map = (ra_inds[time_ind], dec_inds[time_ind], unmasked_inds)
            counts[ind_map] += 1
            map[ind_map] += (weight*data)[time_ind, unmasked_inds]
            if not noise_i is None :
                if not hasattr(weight, '__iter__') :
                    noise_i[ind_map] += weight
                else :
                    noise_i[ind_map] += weight[unmasked_inds]

def pixel_counts(data, ra_inds, dec_inds, counts, map_shape=(-1,-1)) :
    """Counts the hits on each unique pixel.

    Returns pix_list, a list of tuples, each tuple is a (ra,dec) index on a 
    map pixel hit on this scan.  The list only contains unique entries.  The
    array counts (preallowcated for performance), is
    filled with the number of hits on each of these pixels as a function of
    frequency index. Only the entries counts[:len(pix_list), :] are meaningful.
    """

    if ra_inds.shape != dec_inds.shape or ra_inds.ndim != 1 :
        raise ValueError('Ra and Dec arrays not properly shaped.')
    if counts.shape[-1] != data.shape[-1] or counts.shape[0] < len(ra_inds):
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
        counts[pix_list.index(pix), unmasked_freqs] += 1

    return pix_list

def add_scan_noise(pixels, counts, variance, noise_inv) :
    """Adds a scan to the map inverse noise matrix.

    Passed a list of map indices that were hit this scan and an array of the
    number of times each pixel was hit at each frequency.  This is all added to
    noise_inv, inversly wieghted by variance.  Governing formula found in
    Kiyo's notes November 27, 2010.
    """
    
    # Calculate number of good pointings are each frequency.
    n_pix = len(pixels)
    pointings = sp.sum(counts[:n_pix,:], 0)
    f_inds = pointings > 0
    point = pointings[f_inds]
    if hasattr(variance, '__iter__') :
        var = variance[f_inds]
    else :
        var = variance
    for ii, p1 in enumerate(pixels) :
        noise_inv[p1 + p1 + (f_inds,)] += sp.array(counts[ii,f_inds], 
                                                   dtype=float)/var
        for jj, p2 in enumerate(pixels) :
            noise_inv[p1 + p2 + (f_inds,)] -= (sp.array(counts[ii,f_inds], 
                                   dtype=float)*counts[jj,f_inds]/point/var)


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        MapMaker(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        MapMaker().execute()
                


