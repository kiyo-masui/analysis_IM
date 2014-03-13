"""Module subtracts previousely measured foreground modes from the time stream.
"""

import shelve

import numpy.ma as ma
import numpy as np

import kiyopy.custom_exceptions as ce
import base_single
from foreground import ts_measure

import matplotlib.pyplot as plt


class Subtract(base_single.BaseSingle):
    """Pipeline module that subtracts foreground modes from the time stream
    """

    prefix = 'fs_'

    params_init = {'foreground_file': 'default_fname',
                   'n_modes_subtract': 1
                  }
    
    def __init__(self, parameter_file_or_dict=None, feedback=2):
        
        # Call the base_single init.
       	base_single.BaseSingle.__init__(self, parameter_file_or_dict,
                                        feedback)
        # Now open the database file that has the measured foregrounds.
        foreground_file = self.params['foreground_file']
        self.foreground_db = shelve.open(foreground_file, 'r')
        
    def action(self, Data):
        
        params = self.params
        data = Data.data
        n_pol = data.shape[1]
        n_cal = data.shape[2]
        n_chan = data.shape[3]
        # Figure out the data base key for this file.
        file_ind = self.file_ind
        file_middle = self.params['file_middles'][file_ind]
        key = ts_measure.get_key(file_middle)
        # Get the eigenvectors.
        # Shape is (n_bands, n_pols, n_cals, n_chan, n_chan).
        evects = self.foreground_db[key + '.vects']
        if evects.shape[1:] != (n_pol, n_cal, n_chan, n_chan):
            raise ce.DataError("Data and foreground array shapes don't match.")
        # Shape is (n_bands, n_chan).
        freq = self.foreground_db[key + '.freq']
        # Figure out which band the current data block is from.
        Data.calc_freq()
        this_freq = Data.freq
        band_ind = -1
        for ii in range(freq.shape[0]):
            if np.allclose(this_freq, freq[ii,:]):
                band_ind = ii
        if band_ind == -1:
            raise ce.DataError("Frequency axis does not match any IF.")
        evects = evects[ii,...]
        # Subtract the time mean.
        data[...] -= ma.mean(data, 0)
        # Now loop through the polarizations and subtract out the appropiate
        # foreground vectors.
        for ii in range(n_pol):
            for jj in range(n_cal):
                for kk in range(params['n_modes_subtract']):
                    # Get the relevant data.
                    vector = evects[ii,jj,:,-1 - kk]
                    this_data = data[:,ii,jj,:]
                    this_mask = np.logical_not(ma.getmaskarray(this_data))
                    filled_data = this_data.filled(0)
                    # Calculate the amplitude of the foregrounds at each time,
                    # properly treating the mask.
                    amps = np.sum(filled_data * vector, -1)
                    norms = np.sum(this_mask * vector**2, -1)
                    bad_times = norms < 1.e-2
                    norms[bad_times] = 1.
                    amps[bad_times] = 0
                    amps /= norms
                    # Now subtract the forgrounds out.
                    this_data[...] -= amps[:,None] * vector
                    this_data[bad_times,:] = ma.masked

        Data.add_history('Subtracted foreground modes.', 
                         ("Subtracted %d modes." % params['n_modes_subtract']))
        return Data



# If this file is run from the command line, execute the main function.
if __name__=="__main__":
    import sys
    Subtract(str(sys.argv[1])).execute()
