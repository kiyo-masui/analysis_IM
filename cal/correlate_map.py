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

# XXX
import matplotlib.pyplot as plt


params_init = {
               # IO.
               "input_root" : "./testdata/",
               "file_middles" : ("testfile_guppi_combined",),
               "input_end" : ".fits",
               "output_root" : "./",
               "output_filename" : "foreground_modes.shelve",
               "scans" : (),
               "IFs" : (),
               'map_input_root' : './',
               'map_type' : 'clean_map_',
               'map_bands' : (),
               'map_polarizations' : (),
               'interpolation' : 'nearest'
               }

prefix = 'tf_'

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
            norm[norm==0] = 1
            corr /= norm
            #plt.semilogy(h, '.')
            #plt.figure()
            #for ii in range(1,5):
            #    plt.plot(v[:,-ii])
            #plt.show()
            out_db[key + '.gains'] = eigen_vects
            out_db[key + '.freq'] = freq_dict[key]
        out_db.close()


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
                if sp.allclose(maps_freq, freq[ii,:]):
                    maps = band_maps
                    break
            else:
                raise NotImplementedError('No maps with frequency axis exactly'
                                          ' matching data.')
            # Check the polarization axis.
            data_pols = Blocks[0].field['CRVAL4'].copy()
            for ii in range(len(data_pols)):
                if (misc.polint2str(data_pols[ii])
                    != self.params['map_polarizations'][ii]):
                    msg = ('Map polarizations not in same order'
                           ' as data polarizations.')
                    raise NotImplementedError(map)
            # Now process each block.
            for Data in Blocks:
                this_corr, this_norm = get_correlation(Data, band_maps,
                                        interpolation=params['interpolation'])
                corr[ii,...] += corr
                norm[ii,...] += norm
        return key, corr, norm, freq


def get_key(middle):
    # For now just use the session number as the key.
    separate = middle.split('/')[-1]
    sess_num = separate.split('_')[0]
    key = sess_num
    return key

def get_correlation(Data, maps, interpolation='nearest'):
    "Correlates the maps with the data."
    
    n_pols = Data.dims[1]
    if len(maps) != n_pols:
        raise ValueError("Supplied wrong number of maps.")
    # Get the time array (for slope mode subtraction).
    Data.calc_time()
    time = Data.time
    # Initialize outputs.
    correlation = sp.zeros(Data.dims[1:], dtype=float)
    normalization = sp.zeros(Data.dims[1:], dtype=float)
    for ii in range(n_pols):
        map = maps[ii]
        if map.size[0] != Data.dims[3]:
            raise RuntimeError("Map and data frequency axes not the same
                               " length.")
        Data.calc_pointing()
        Data.calc_freq()
        # Figure out which pointings (times) are inside the map bounds.
        map_ra = Map.get_axis('ra')
        map_dec = Map.get_axis('dec')
        on_map_inds = sp.logical_and(
            sp.logical_and(Data.ra >= min(map_ra), Data.ra <= max(map_ra)),
            sp.logical_and(Data.dec >= min(map_dec), Data.dec <= max(map_dec)))
        # Convert map to a time domain array.
        submap = sp.empty((sp.sum(on_map_inds), map.shape[0]), dtype=float)
        jj = 0
        for ii in range(len(on_map_inds)) :
            if on_map_inds[ii] :
                submap[jj,:] = map.slice_interpolate([1, 2], 
                        [Data.ra[ii], Data.dec[ii]], kind=interpolation)
                jj += 1
        # Now get the corresponding data.
        subdata = Data.data[on_map_inds,ii,:,:]
        un_mask = sp.logical_not(ma.getmaskarray(subdata))
        subdata = subdata.filled(0)
        # Broadcast map data up to the same shape as the time stream (add cal
        # axis).
        submap = sp.zeros_like(subdata) + submap[:,None,:]
        # We do not want to correlate the mean mode and the slope mode, so
        # subtract these out.
        # First the mean mode.
        submap -= sp.sum(submap * unmask, 0) / sp.sum(un_mask, 0)
        subdata -= sp.sum(subdata * unmask, 0) / sp.sum(un_mask, 0)
        # Now for the slope.  No need to ensure the slope mode that is 
        # orthoganol to the mean mode (which is not ensured because of the
        # mask) because we are zeroing both. 
        
        # Calculate the correlation and the normalization.
        correlation[ii,:,:] = sp.sum(submap * un_mask * subdata, 0)
        normilization[ii,:,:] = sp.sum(submap * un_mask * submap, 0)
    return correlation, normalization
