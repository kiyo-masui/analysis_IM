"""Module that subtracts map data from time_stream data, leaving only noise.
"""

import scipy as sp
import numpy.ma as ma

import kiyopy.custom_exceptions as ce
import base_single
import map.tools
from core import fits_map

class Subtract(base_single.BaseSingle) :
    """Pipeline module subtracts a map from time stream data.

    This module reads in a map and times stream data.  It then subtracts the
    signal part off of each time bin using the pointing information and the
    map.  This should leave only noise and map residuals (signal that isn't in
    the map).
    """
    
    prefix = 'sm_'
    params_init = {
                   'map_file' : 'testfile_map.fits'
                   # XXX: how to deal with time mean? (often this doesn't
                   # matter.)
                   }

    # Add extra stuff to the constructor.
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
        base_single.BaseSingle.__init__(self, parameter_file_or_dict,
                                        feedback)
        # Read in the calibration file.
        map_file_name = self.params['map_file']
        self.Map = fits_map.read(map_file_name, 1, feedback=self.feedback)

    def action(self, Data) :
        sub_map(Data, self.Map)
        Data.add_history('Subtracted map from data.', 
                         ('Map file: ' + self.params['map_file'],))

        return Data

def sub_map(Data, Map) :
    """Subtracts a Map out of Data."""

    # Some dimension checks.  Eventually may want to have a tuple of maps, one
    # for each polaization.  For now, only use I.
    if ((len(Data.field['CRVAL4']) > 1) or 
        (hasattr(Map, '__iter__') and len(Map) > 1)):
        raise NotImplementedError('Multiple polarizations not supported.')
    pol_ind = 0
    if hasattr(Map, '__iter__') :
        Map = Map[0]
    if Map.field['POL'].item() != Data.field['CRVAL4'][0] :
        raise ce.DataError("Polarization types don't match.")
        
    Data.calc_pointing()
    Data.calc_freq()
    centre, shape, spacing = map.tools.get_map_params(Map)
    nt_data = Data.dims[0]
    nf_data = Data.dims[3]
    # These indices are the length of the time axis. Integer indicies.
    ra_ind = map.tools.calc_inds(Data.ra, centre[0], shape[0], spacing[0])
    dec_ind = map.tools.calc_inds(Data.dec, centre[1], shape[1], spacing[1])
    # Length of the frequency axis.
    freq_ind = map.tools.calc_inds(Data.freq, centre[2], shape[2], spacing[2])
    # Exclude indices that are off map or out of band. Boolian indices.
    on_map_inds = sp.logical_and(sp.logical_and(ra_ind>=0, ra_ind<shape[0]),
                                 sp.logical_and(dec_ind>=0, dec_ind<shape[1]))
    in_band_inds = sp.logical_and(freq_ind >= 0, freq_ind < shape[2])
    # Broadcast to the same shape and combine.
    covered_inds = sp.logical_and(on_map_inds[:, sp.newaxis], 
                                  in_band_inds[sp.newaxis, :])
    # Subtract the map from the data.  Mask the covered indicies.
    submap = Map.data[ra_ind[on_map_inds], dec_ind[on_map_inds], :]
    submap = submap[:, freq_ind[in_band_inds]]
    for cal_ind in range(Data.dims[2]) :
        Data.data[:,pol_ind, cal_ind, :][covered_inds] = (
            Data.data[:,pol_ind, cal_ind, :][covered_inds] - submap.flatten())
        Data.data[:,pol_ind, cal_ind, :][sp.logical_not(covered_inds)] = (
            ma.masked)
    

    
    
