"""Module measure the correlation of blocks of time stream data with the map.

The idea here is to force a constant *relative* calibration from session to
session.
"""

import shelve
import multiprocessing as mp
import time as time_module

import numpy as np
import scipy.linalg as linalg
import numpy.ma as ma

from kiyopy import parse_ini, utils
import kiyopy.pickle_method
import kiyopy.utils
import kiyopy.custom_exceptions as ce
import core.fitsGBT
from utils import misc

# XXX
import matplotlib.pyplot as plt


params_init = {
               # IO.
               "input_root" : "./testdata/",
               "file_middles" : ("testfile_guppi_combined",),
               "input_end" : ".fits",
               "output_root" : "./",
               "output_filename" : "map_correlations.shelve",
               "scans" : (),
               "IFs" : (),
               'map_input_root' : './',
               'map_type' : 'clean_map_',
               'map_bands' : (),
               'map_polarizations' : (),
               'interpolation' : 'nearest',
               'smooth_modes_subtract' : 2,
               }

prefix = 'mc_'

class Measure(object) :
    """Measures the correlation of the data with the map.
    """

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                      prefix=prefix, feedback=feedback)
        self.feedback = feedback

        # Read in the map files.
        map_fnames_start = (self.params['map_input_root']
                            + self.params['map_type'])
        self.maps = []
        for band in self.params['map_bands']:
            this_band_maps = []
            for pol in self.params['map_polarizations']:
                map_file_name = (map_fnames_start + pol + '_' + str(band)
                                 + '.npy')
                map = algebra.load(map_file_name)
                map = algebra.make_vect(map)
                this_band_maps.append(map)
            self.maps.append(this_band_maps)


    def execute(self, nprocesses=1) :
                
        params = self.params
        kiyopy.utils.mkparents(params['output_root'] + 
                               params['output_filename'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        file_middles = params['file_middles']
        # Store the correlation and normalization session by session in a
        # dictionary.
        corr_dict = {}
        norm_dict = {}
        # Also store the frequency axis.
        freq_dict = {}
        # Finally, a gain dictionary.
        gain_dict = {}
        # Loop though all the files and accumulate the correlation and the
        # normalization.  Files in the same session are summed together.
        for middle in file_middles:
            key, corr, norm, freq = self.process_file(middle)
            if corr_dict.has_key(key):
                if corr_dict[key].shape != corr.shape:
                    msg = ("All data needs to have the same band and"
                           "polarization structure.")
                    raise ce.DataError(msg)
                corr_dict[key] += corr
                norm_dict[key] += norm
                if not np.allclose(freq_dict[key], freq):
                    raise ce.DataError("Frequency structure not consistant.")
            else:
                corr_dict[key] = corr
                norm_dict[key] = norm
                freq_dict[key] = freq
        # Now that we have the correlation summed for all files in each
        # session, normalize it and store it as an output.
        output_fname = params['output_root'] + params["output_filename"]        
        out_db = shelve.open(output_fname)
        # Loop through all the data divisions and processes them one at a
        # time.
        for key in corr_dict.iterkeys():
            corr = corr_dict[key]
            norm = norm_dict[key]
            # Normalize.
            corr[norm==0] = 0
            norm[norm==0] = 1
            gains = corr / norm
            gain_dict[key] = gains
            #plt.semilogy(h, '.')
            #plt.figure()
            #for ii in range(1,5):
            #    plt.plot(v[:,-ii])
            #plt.show()
            out_db[key + '.gains'] = gains
            out_db[key + '.freq'] = freq_dict[key]
        out_db.close()
        #### Apply the calibration to the data. ####
        for middle in file_middles:
            key = get_key(middle)
            gain = gain_dict[key]
            self.calibrate_file(middle, gain)

    def process_file(self, middle):
        params = self.params
        file_name = (params['input_root'] + middle
                     + params['input_end'])
        band_inds = params["IFs"]
        Reader = core.fitsGBT.Reader(file_name, feedback=self.feedback)
        n_bands = len(Reader.IF_set)
        if not band_inds:
            band_inds = range(n_bands)
        # Number of bands we acctually process.
        n_bands_proc = len(band_inds)
        if not band_inds:
            band_inds = range(n_bands)
        # Number of bands we acctually process.
        n_bands_proc = len(band_inds)
        # Get the key that will group this file with other files.
        key = get_key(middle)
        # Read one block to figure out how many polarizations and channels
        # there are.
        Data = Reader.read(0,0)
        n_pol = Data.dims[1]
        n_cal = Data.dims[2]
        n_chan = Data.dims[3]
        # Allowcate memory for the outputs.
        corr = np.zeros((n_bands_proc, n_pol, n_cal, n_chan),
                         dtype=float)
        norm = np.zeros(corr.shape, dtype=int)
        freq = np.empty((n_bands_proc, n_chan))
        for ii in range(n_bands_proc):
            Blocks = Reader.read((), ii)
            Blocks[0].calc_freq()
            freq[ii,:] = Blocks[0].freq
            # We are going to look for an exact match in for the map 
            # frequencies. This could be made more general since the sub_map
            # function can handle partial overlap, but this will be fine for
            # now.
            for band_maps in self.maps:
                maps_freq = band_maps[0].get_axis('freq')
                if np.allclose(maps_freq, freq[ii,:]):
                    maps = band_maps
                    break
            else:
                raise NotImplementedError('No maps with frequency axis exactly'
                                          ' matching data.')
            # Check the polarization axis. If the same number of maps where
            # passed, check that the polarizations are in order.  If only one
            # map was passed, correlate all data polarizations against it.
            data_pols = Blocks[0].field['CRVAL4'].copy()
            if len(band_maps == 1):
                maps_to_correlate = band_maps * len(data_pols)
            else:
                for ii in range(len(data_pols)):
                    if (misc.polint2str(data_pols[ii])
                        != self.params['map_polarizations'][ii]):
                        msg = ('Map polarizations not in same order'
                               ' as data polarizations.')
                        raise NotImplementedError(map)
                maps_to_correlate = band_maps
            # Now process each block.
            for Data in Blocks:
                this_corr, this_norm = get_correlation(Data,
                                maps_to_correlate,
                                interpolation=params['interpolation'],
                                modes_subtract=params['smooth_modes_subtract'])
                corr[ii,...] += corr
                norm[ii,...] += norm
        return key, corr, norm, freq

    def calibrate_file(middle, gains):
        # This function is largely cut and pasted from process file. I should
        # really combine the code into an iterator but that's a lot of work.
        # Alternativly, I could make a meta function and pass a function to it.
        params = self.params
        file_name = (params['input_root'] + middle
                     + params['input_end'])
        band_inds = params["IFs"]
        Reader = core.fitsGBT.Reader(file_name, feedback=self.feedback)
        n_bands = len(Reader.IF_set)
        if not band_inds:
            band_inds = range(n_bands)
        # Number of bands we acctually process.
        n_bands_proc = len(band_inds)
        if not band_inds:
            band_inds = range(n_bands)
        # Number of bands we acctually process.
        n_bands_proc = len(band_inds)
        # Get the key that will group this file with other files.
        key = get_key(middle)
        # Read one block to figure out how many polarizations and channels
        # there are.
        Data = Reader.read(0,0)
        n_pol = Data.dims[1]
        n_cal = Data.dims[2]
        n_chan = Data.dims[3]
        # Allowcate memory for the outputs.
        corr = np.zeros((n_bands_proc, n_pol, n_cal, n_chan),
                         dtype=float)
        norm = np.zeros(corr.shape, dtype=int)
        freq = np.empty((n_bands_proc, n_chan))
        for ii in range(n_bands_proc):
            Blocks = Reader.read((), ii)
            Blocks[0].calc_freq()
            freq[ii,:] = Blocks[0].freq
            # We are going to look for an exact match in for the map 
            # frequencies. This could be made more general since the sub_map
            # function can handle partial overlap, but this will be fine for
            # now.
            for band_maps in self.maps:
                maps_freq = band_maps[0].get_axis('freq')
                if np.allclose(maps_freq, freq[ii,:]):
                    maps = band_maps
                    break
            else:
                raise NotImplementedError('No maps with frequency axis exactly'
                                          ' matching data.')
            # Check the polarization axis. If the same number of maps where
            # passed, check that the polarizations are in order.  If only one
            # map was passed, correlate all data polarizations against it.
            data_pols = Blocks[0].field['CRVAL4'].copy()
            if len(band_maps == 1):
                maps_to_correlate = band_maps * len(data_pols)
            else:
                for ii in range(len(data_pols)):
                    if (misc.polint2str(data_pols[ii])
                        != self.params['map_polarizations'][ii]):
                        msg = ('Map polarizations not in same order'
                               ' as data polarizations.')
                        raise NotImplementedError(map)
                maps_to_correlate = band_maps
            # Now process each block.
            for Data in Blocks:
                this_corr, this_norm = get_correlation(Data,
                                maps_to_correlate,
                                interpolation=params['interpolation'],
                                modes_subtract=params['smooth_modes_subtract'])
                corr[ii,...] += corr
                norm[ii,...] += norm
        return key, corr, norm, freq



def get_key(middle):
    # For now just use the session number as the key.
    separate = middle.split('/')[-1]
    sess_num = separate.split('_')[0]
    key = sess_num
    return key

def get_correlation(Data, maps, interpolation='nearest', modes_subtract=2):
    "Correlates the maps with the data."
    
    n_pols = Data.dims[1]
    if len(maps) != n_pols:
        raise ValueError("Supplied wrong number of maps.")
    # Get the time array (for slope mode subtraction).
    Data.calc_time()
    time = Data.time
    # Initialize outputs.
    correlation = np.zeros(Data.dims[1:], dtype=float)
    normalization = np.zeros(Data.dims[1:], dtype=float)
    for ii in range(n_pols):
        map = maps[ii]
        if map.shape[0] != Data.dims[3]:
            raise RuntimeError("Map and data frequency axes not the same"
                               " length.")
        Data.calc_pointing()
        Data.calc_freq()
        # Figure out which pointings (times) are inside the map bounds.
        map_ra = map.get_axis('ra')
        map_dec = map.get_axis('dec')
        on_map_inds = np.logical_and(
            np.logical_and(Data.ra >= min(map_ra), Data.ra <= max(map_ra)),
            np.logical_and(Data.dec >= min(map_dec), Data.dec <= max(map_dec)))
        # Convert map to a time domain array.
        submap = np.empty((np.sum(on_map_inds), map.shape[0]), dtype=float)
        kk = 0
        for jj in range(len(on_map_inds)) :
            if on_map_inds[jj] :
                submap[kk,:] = map.slice_interpolate([1, 2], 
                        [Data.ra[jj], Data.dec[jj]], kind=interpolation)
                kk += 1
        # Now get the corresponding data.
        subdata = Data.data[on_map_inds,ii,:,:]
        un_mask = np.logical_not(ma.getmaskarray(subdata))
        subdata = subdata.filled(0)
        # Broadcast map data up to the same shape as the time stream (add cal
        # axis).
        submap = np.zeros_like(subdata) + submap[:,None,:]
        # Get rid of the low frequency, smooth components by subtracting out
        # basis polynomials.
        # Generate basis polynomials that are orthnormal given the mask.
        on_map_time = time[on_map_inds]
        polys = misc.ortho_poly(on_map_time[:,None,None], modes_subtract,
                                un_mask, 0)
        # Subtract out of the data.
        mags = np.sum(subdata * un_mask * polys, 1)
        to_subtract = np.sum(mags[:,None,...] * polys, 0)
        subdata -= to_subtract
        # Subtract out of the map.
        mags = np.sum(submap * un_mask * polys, 1)
        to_subtract = np.sum(mags[:,None,...] * polys, 0)
        submap -= to_subtract
        # Calculate the correlation and the normalization.
        correlation[ii,:,:] = np.sum(submap * un_mask * subdata, 0)
        normalization[ii,:,:] = np.sum(submap * un_mask * submap, 0)
    return correlation, normalization
