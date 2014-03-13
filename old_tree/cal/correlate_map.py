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

from core import algebra
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
               "output_end" : ".fits",
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
               'diff_gain_cal_only' : False,
               # Low pass filtering, options are 'gaussian' and 'edge'.
               'filter_type': 'edge'
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
        n_new = nprocesses-1  # How many new processes to spawn at once.
        n_files = len(file_middles)
        if n_new > 0:
            # Multiprocessed version.
            process_list = range(n_new)
            pipe_list = range(n_new)
            for ii in xrange(n_files + n_new) :
                if ii >= n_new :
                    # First end a process before starting a new one.
                    key, corr, norm, freq =  pipe_list[ii%n_new].recv()
                    if corr_dict.has_key(key):
                        if corr_dict[key].shape != corr.shape:
                            msg = ("All data needs to have the same band and"
                                   "polarization structure.")
                            raise ce.DataError(msg)
                        corr_dict[key] += corr
                        norm_dict[key] += norm
                        if not np.allclose(freq_dict[key], freq):
                            raise ce.DataError("Frequency structure not "
                                               "consistant.")
                    else:
                        corr_dict[key] = corr
                        norm_dict[key] = norm
                        freq_dict[key] = freq
                if ii < n_files :
                    # Start a new process.
                    Here, Far = mp.Pipe()
                    pipe_list[ii%n_new] = Here
                    process_list[ii%n_new] = mp.Process(
                        target=self.process_file, args=(file_middles[ii], Far))
                    process_list[ii%n_new].start()
        else:
            # Single process.
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
            corr[norm==0] = 1
            norm[norm==0] = 1
            gains = corr / norm
            gain_dict[key] = gains
            #plt.figure()
            #if params['diff_gain_cal_only']:
            #    plt.plot(freq_dict[key][0,:],
            #             (gains[0,0,0,:] + gains[0,0,1,:])/2., '.')
            #    plt.plot(freq_dict[key][0,:],
            #             (gains[0,3,0,:] + gains[0,3,1,:])/2., '.')
            #else:
            #    plt.plot(freq_dict[key][0,:], gains[0,0,0,:], '.')
            #plt.title(key)
            #plt.xlabel('Frequency (Hz)')
            #plt.ylabel('Correlation amplitude')
            out_db[key + '.gains'] = gains
            out_db[key + '.freq'] = freq_dict[key]
        #if not params['diff_gain_cal_only']:
        #    plt.show()
        out_db.close()
        #### Apply the calibration to the data. ####
        for middle in file_middles:
            key = get_key(middle)
            gain = gain_dict[key]
            freq = freq_dict[key]
            self.calibrate_file(middle, gain, freq)

    def process_file(self, middle, Pipe=None):
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
        norm = np.zeros(corr.shape, dtype=float)
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
            if len(band_maps) == 1:
                maps_to_correlate = band_maps * len(data_pols)
            else:
                for jj in range(len(data_pols)):
                    if (misc.polint2str(data_pols[jj])
                        != self.params['map_polarizations'][jj]):
                        msg = ('Map polarizations not in same order'
                               ' as data polarizations.')
                        raise NotImplementedError(map)
                maps_to_correlate = band_maps
            # Now process each block.
            for Data in Blocks:
                this_corr, this_norm = get_correlation(Data,
                                maps_to_correlate,
                                interpolation=params['interpolation'],
                                modes_subtract=params['smooth_modes_subtract'],
                                filter_type=params['filter_type'])
                # Check that the answer we got is sane, if not, throw away the
                # this set.
                tmp_corr = this_corr.copy()
                tmp_norm = this_norm.copy()
                tmp_corr[tmp_norm == 0] = 1.
                tmp_norm[tmp_norm == 0] = 1.
                tmp_gains = tmp_corr / tmp_norm
                if np.all(tmp_gains < 2) and np.all(tmp_gains > 0.5):
                    corr[ii,...] += this_corr
                    norm[ii,...] += this_norm
                else:
                    pass
        if Pipe is None:
            return key, corr, norm, freq
        else:
            Pipe.send((key, corr, norm, freq))

    def calibrate_file(self, middle, gain, freq):
        # This function is largely cut and pasted from process file. I should
        # really combine the code into an iterator but that's a lot of work.
        # Alternativly, I could make a meta function and pass a function to it.
        params = self.params
        file_name = (params['input_root'] + middle
                     + params['input_end'])
        # Output parameters.
        Writer = core.fitsGBT.Writer(feedback=self.feedback)
        out_filename = (params['output_root'] + middle
                        + params['output_end'])
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
            if len(band_maps) == 1:
                maps_to_correlate = band_maps * len(data_pols)
            else:
                for jj in range(len(data_pols)):
                    if (misc.polint2str(data_pols[jj])
                        != self.params['map_polarizations'][jj]):
                        msg = ('Map polarizations not in same order'
                               ' as data polarizations.')
                        raise NotImplementedError(map)
                maps_to_correlate = band_maps
            # Now process each block.
            for Data in Blocks:
                if params['diff_gain_cal_only']:
                    if tuple(Data.field['CRVAL4']) != (-5, -7, -8, -6):
                        msg = ("Expected polarizations to be ordered "
                               "(XX, XY, YX, YY).")
                        raise NotImplementedError(msg)
                    Data.data[:,0,:,:] /= gain[ii,0,:,:]
                    Data.data[:,3,:,:] /= gain[ii,3,:,:]
                    cross_gain = np.sqrt(gain[ii,0,:,:] * gain[ii,3,:,:])
                    Data.data[:,1,:,:] /= cross_gain
                    Data.data[:,2,:,:] /= cross_gain
                else:
                    Data.data /= gain[ii,...]
                Writer.add_data(Data)
        # Write the data back out.
        utils.mkparents(out_filename)
        Writer.write(out_filename)
        

def get_key(middle):
    # For now just use the session number as the key.
    separate = middle.split('/')[-1]
    separate = separate.split('_')
    sess_num = separate[0]
    scan_num = int(separate[-1].split('-')[0])
    scan_set_num = scan_num//8
    div = scan_set_num % 2
    round = scan_set_num // 10
    #key = sess_num + "_" + str(round) + "_" + str(div)
    key = sess_num
    return key

def get_correlation(Data, maps, interpolation='nearest', modes_subtract=2,
                    filter_type='edge'):
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
        # If none of the scan is on the map, skip.
        if not np.any(on_map_inds):
            continue
        # Convert map to a time domain array.
        submap = np.empty((np.sum(on_map_inds), map.shape[0]), dtype=float)
        kk = 0
        for jj in range(len(on_map_inds)) :
            if on_map_inds[jj] :
                submap[kk,:] = map.slice_interpolate([1, 2], 
                        [Data.ra[jj], Data.dec[jj]], kind=interpolation)
                kk += 1
        n_time_on = kk
        # Set up filter parameters.
        if filter_type == 'edge':
            total_modes = modes_subtract
            subtract_weights = np.ones(total_modes)
        if filter_type == 'gaussian' or filter_type == 'gaussian/edge':
            total_modes = min((4*modes_subtract, n_time_on))
            # Gaussian taper to filtering.
            subtract_weights = np.exp(-(np.arange(total_modes, dtype=float)
                                      / modes_subtract)**2 / 2.)
            if filter_type == 'gaussian/edge':
                subtract_weights[0:2] = 1.
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
        # Test if the mask is the same for all slices.  If it is, that greatly
        # reduces the work as we only have to generate one set.
        all_masks_same = True
        for jj in range(n_time_on):
            if np.all(un_mask[jj,...] == un_mask[jj,0,0]):
                continue
            else:
                all_masks_same = False
                break
        if all_masks_same:
            polys = misc.ortho_poly(on_map_time, total_modes,
                                    un_mask[:,0,0], 0)
            polys.shape = (total_modes, n_time_on, 1, 1)
        else:
            polys = misc.ortho_poly(on_map_time[:,None,None], total_modes,
                                    un_mask, 0)
        # Subtract out of the data.
        mags_data = np.sum(subdata * un_mask * polys, 1)
        # If using a taper, add that in.
        mags_data *= subtract_weights[:,None,None]
        to_subtract_data = np.sum(mags_data[:,None,...] * polys, 0)
        subdata -= to_subtract_data
        # Subtract out of the map.
        mags_map = np.sum(submap * un_mask * polys, 1)
        mags_map *= subtract_weights[:,None,None]
        to_subtract_map = np.sum(mags_map[:,None,...] * polys, 0)
        submap -= to_subtract_map
        # Calculate the correlation and the normalization.
        corr = np.sum(submap * un_mask * subdata, 0)
        norm = np.sum(submap * un_mask * submap, 0)
        # Calculate inverse reduced Chi-squared and weight this measurement by
        # it. Be carefull about the 0 information case (all masked).
        filled_norm = norm.copy()
        # No information.
        bad_inds = np.logical_or(norm == 0, np.sum(un_mask, 0) <
                                 2 * modes_subtract)
        filled_norm[bad_inds] = 1
        amp = corr / filled_norm
        fit = submap * amp
        inv_chi_sq = np.sum((subdata - fit)**2 * un_mask, 0)
        inv_chi_sq[bad_inds] = 1.
        inv_chi_sq = np.sum(un_mask, 0) / inv_chi_sq
        inv_chi_sq[bad_inds] = 0
        corr *= inv_chi_sq
        norm *= inv_chi_sq
        # Store results in output arrays
        correlation[ii,:,:] = corr
        normalization[ii,:,:] = norm
        if False and (ii == 0) and int(Data.field['SCAN']) == 86:
            print correlation[ii,:,:]
            print normalization[ii,:,:]
            print inv_chi_sq
            print correlation[ii,:,:] / normalization[ii,:,:]
            plt.plot(on_map_time, (subdata * un_mask)[:,0,0], '.')
            plt.plot(on_map_time, (fit * un_mask)[:,0,0], '.')
            #plt.plot((subdata * un_mask)[:,0,0] - (submap * un_mask)[:,0,0],
            #           '.')
            #plt.plot((subdata * un_mask)[:,0,3] - (submap * un_mask)[:,0,3],
            #           '.')
            plt.show()
    return correlation, normalization
