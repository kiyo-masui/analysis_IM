"""Module that subtracts map data from time_stream data, leaving only noise.
"""

import cPickle

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import kiyopy.utils as ku
import base_single
import map.tools
from core import algebra

class Subtract(base_single.BaseSingle) :
    """Pipeline module subtracts a map from time stream data.

    This module reads in a map and times stream data.  It then subtracts the
    signal part off of each time bin using the pointing information and the
    map.  This should leave only noise and map residuals (signal that isn't in
    the map).
    """
    
    prefix = 'sm_'
    params_init = {
                   'map_file' : 'testfile_map.fits',
                   'solve_for_gain' : False,
                   # Empty string to not write an output.
                   'gain_output_end' : '',
                   'interpolation' : 'nearest'
                   # XXX: What about if I start subtracing off a linear piece.
                   # TODO: Option to rescaleoutput data to map gain.
                   }

    # Add extra stuff to the constructor.
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
        base_single.BaseSingle.__init__(self, parameter_file_or_dict,
                                        feedback)
        # Read in the calibration file.
        map_file_name = self.params['map_file']
        self.Map = algebra.load(map_file_name)
        self.Map = algebra.make_vect(self.Map)

    def action(self, Data) :
        if (not self.params['solve_for_gain'] or
            self.params['gain_output_end'] is '') :
            sub_map(Data, self.Map, self.params['solve_for_gain'],
                    interpolation=self.params['interpolation'])
        else :
            block_gain = {}
            Data.calc_freq()
            block_gain['freq'] = sp.copy(Data.freq)
            block_gain['time'] = Data.field['DATE-OBS'][0]
            block_gain['scan'] = Data.field['SCAN']
            block_gain['gain'] = sub_map(Data, self.Map, True)
            self.gain_list.append(block_gain)

        Data.add_history('Subtracted map from data.', 
            ('Map file: ' + ku.abbreviate_file_path(self.params['map_file']),))
        return Data

    # Overwrite the base single process file method so we can also pickle the
    # gains we've solved for.
    def process_file(self, file_ind) :
        # Initialize a list to hold the gains calculated in sub_map.
        self.gain_list = []
        # Do the normal thing from the base class.
        base_single.BaseSingle.process_file(self, file_ind)
        # Pickle the gain list.
        if (self.params['solve_for_gain'] and
            not self.params['gain_output_end'] is '') :
            file_middle = self.params['file_middles'][file_ind]
            gain_fname = (self.params['output_root'] + file_middle
                          + self.params['gain_output_end'])
            f = open(gain_fname, 'w')
            cPickle.dump(self.gain_list, f, 0)

def sub_map(Data, Maps, correlate=False, pols=(), make_plots=False,
            interpolation='nearest') :
    """Subtracts a Map out of Data."""
    
    # Import locally since many machines don't have matplotlib.
    if make_plots :
        import matplotlib.pyplot as plt
    
    # Convert pols to an interable.
    if pols is None :
        pols = range(Data.dims[1])
    elif not hasattr(pols, '__iter__') :
        pols = (pols, )
    elif len(pols) == 0 :
        pols = range(Data.dims[1])
    # If solving for gains, need a place to store them.
    if correlate :
        out_gains = sp.empty((len(pols),) + Data.dims[2:4])
    for pol_ind in pols :
        # Check if there one map was passed or multiple.
        if isinstance(Maps, list) or isinstance(Maps, tuple) :
            if len(Maps) != len(pols) :
                raise ValueError("Must provide one map, or one map per "
                                 "polarization.")
            Map = Maps[pol_ind]
        else :
            Map = Maps
        if not Map.axes == ('freq', 'ra', 'dec') :
            raise ValueError("Expected map axes to be ('freq', 'ra', 'dec').")
        Data.calc_pointing()
        Data.calc_freq()
        # Map Parameters.
        centre = (Map.info['freq_centre'], Map.info['ra_centre'],
                  Map.info['dec_centre'])
        shape = Map.shape
        spacing = (Map.info['freq_delta'], Map.info['ra_delta'], 
                   Map.info['dec_delta'])
        # Nearest code is depricated.  We could just use the general code.
        if interpolation == 'nearest' :
            # These indices are the length of the time axis. Integer indicies.
            ra_ind = map.tools.calc_inds(Data.ra, centre[1], shape[1],
                                         spacing[1])
            dec_ind = map.tools.calc_inds(Data.dec, centre[2], shape[2],
                                          spacing[2])
            # Exclude indices that are off map or out of band. Boolian indices.
            on_map_inds = sp.logical_and(
                                 sp.logical_and(ra_ind>=0, ra_ind<shape[1]),
                                 sp.logical_and(dec_ind>=0, dec_ind<shape[2]))
            # Make an array of map data the size of the time stream data.
            submap = Map[:, ra_ind[on_map_inds], dec_ind[on_map_inds]]
        else :
            map_ra = Map.get_axis('ra')
            map_dec = Map.get_axis('dec')
            on_map_inds = sp.logical_and(
                sp.logical_and(Data.ra > min(map_ra), Data.ra < max(map_ra)),
                sp.logical_and(Data.dec > min(map_dec), Data.dec<max(map_dec)))
            submap = sp.empty((Map.shape[0], sp.sum(on_map_inds)), dtype=float)
            jj = 0
            for ii in range(len(on_map_inds)) :
                if on_map_inds[ii] :
                    submap[:, jj] = Map.slice_interpolate([1, 2], 
                            [Data.ra[ii], Data.dec[ii]], kind=interpolation)
                    jj += 1
        # Length of the data frequency axis.
        freq_ind = map.tools.calc_inds(Data.freq, centre[0], shape[0], 
                                       spacing[0])
        in_band_inds = sp.logical_and(freq_ind >= 0, freq_ind < shape[0])
        submap = submap[freq_ind[in_band_inds], ...]
        # Broadcast to the same shape and combine.
        covered_inds = sp.logical_and(on_map_inds[:, sp.newaxis], 
                                      in_band_inds[sp.newaxis, :])
        # submap is the size of the data that is on the map.  Expand to full 
        # size of data.
        subdata = sp.zeros(sp.shape(covered_inds))
        subdata[covered_inds] = sp.rollaxis(submap, 1, 0).flatten()
        subdata[sp.logical_not(covered_inds)] = 0.0
        # Now start using the actual data.  Loop over cal and pol indicies.
        for cal_ind in range(Data.dims[2]) :
            data = Data.data[:,pol_ind, cal_ind, :]
            data[sp.logical_not(covered_inds)] = ma.masked
            # Find the common good indicies.
            un_mask = sp.logical_not(data.mask)
            # Find the number of good indicies at each frequency.
            counts = sp.sum(un_mask, 0)
            counts[counts == 0] = -1
            # Subtract out the mean from the map.
            tmp_subdata = (subdata - sp.sum(un_mask*subdata, 0)/counts)
            # Correlate to solve for an unknown gain.
            if correlate :
                tmp_data = data.filled(0.0)
                tmp_data = (tmp_data - sp.sum(un_mask*data, 0)
                            / counts)
                gain = (sp.sum(un_mask*tmp_subdata*tmp_data, 0) / 
                        sp.sum(un_mask*tmp_subdata*tmp_subdata, 0))
                gain[counts == -1] = 0.0
                out_gains[pol_ind,cal_ind,:] = gain
            else :
                gain = 1.0
            # Now do the subtraction and mask the off map data.  We use the
            # mean subtracted map, to preserve data mean.
            if make_plots :
                plt.figure()
                plt.plot(ma.mean((gain*tmp_subdata), -1), '.b')
                plt.plot(ma.mean((tmp_subdata), -1), '.r')
                plt.plot(ma.mean((data - ma.mean(data, 0)), -1), '.g')
                #plt.plot(ma.mean((data), -1), '.g')
                #plt.plot((gain*tmp_subdata)[:, 45], '.b')
                #plt.plot((data - ma.mean(data, 0))[:, 45], '.g')
            data[...] -= gain*tmp_subdata
    if correlate :
        return out_gains

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Subtract(str(sys.argv[1])).execute()
