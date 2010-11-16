
"""Program that calculates the correlation function across frequency slices.
"""

import scipy as sp
import matplotlib.pyplot as plt

from kiyopy import parse_ini
import kiyopy.custom_exceptions as ce
from core import fits_map

params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_map",),
               'input_end' : ".fits",
               'output_root' : "./testoutput",
               # What frequencies to correlate:
               'freq' : (15, 40, 65, 90, 115),
               # Angular lags at which to calculate the correlation.  Upper
               # edge bins in degrees.
               'lags' : tuple(sp.arange(0.1, 1., 0.1))
               }
prefix = 'fs_'

class FreqSlices(object) :
    
    def __init__(self, parameter_file_or_dict=None) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix)

    def execute(self) :
        params = self.params
        # Write parameter file.
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        # Get the map data from file.
        if len(params['file_middles']) != 1 :
            raise ce.FileParameterTypeError('For now can only process one'
                                            'file.')
        file_name = (params['input_root'] + params['file_middles'][0] + 
                     params['input_end'])
        Map = fits_map.read(file_name)
        Map.calc_axes()
        # Set up correlation function.
        lags = sp.array(params['lags'])
        freq1 = sp.array(params['freq'], dtype=int)
        freq2 = sp.array(params['freq'], dtype=int)
        nlags = len(lags)
        nf = len(freq1)
        corr = sp.zeros((nf, nf, nlags), dtype=float)
        counts = sp.zeros(corr.shape, dtype=int)
        # Noting that if DEC != 0, then a degree of RA is less than a degree.
        ra_fact = sp.cos(sp.pi*Map.lat[len(Map.lat)//2]/180.0)
        
        # Double loop over all three coordinates. Using only a 2 space indent
        # because this loop is 6 deep.
        for ira1 in range(Map.dims[0]) :
          # Get to cut one loop in half to avoid repeats.
          for ira2 in range(ira1, Map.dims[0]) :
            for jdec1 in range(Map.dims[1]) :
              # This loop must start at zero to get all pairs.
              for jdec2 in range(Map.dims[1]) :
                if (sp.any(Map.data.mask[ira1,jdec1,freq1]) or 
                    sp.any(Map.data.mask[ira2,jdec2,freq2])) :
                        continue
                # Figure out the lag then which bin it falls into.
                dra = (Map.long[ira2] - Map.long[ira1])*ra_fact
                ddec = Map.lat[jdec2] - Map.lat[jdec1]
                lag = sp.sqrt(dra**2 + ddec**2)
                lag_ind = 0
                for this_lag_ind, this_lag_bin in enumerate(lags) :
                    if lag > this_lag_bin :
                       lag_ind = this_lag_ind + 1
                    else :
                        break
                # Take the outer product accross frequencies.
                if lag_ind < nlags :
                    data1 = Map.data[ira1,jdec1,freq1]
                    data2 = Map.data[ira2,jdec2,freq2]
                    this_corr = (data1[:,sp.newaxis] * data2[sp.newaxis,:])
                    corr[:,:,lag_ind] += this_corr.filled(0.)
                    counts[:,:,lag_ind] += sp.logical_not(this_corr.mask)
        corr /= counts - 1
        norms = corr[:,:,0].diagonal()
        norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
        corr /= norms[:,:,sp.newaxis]
        
        self.lags = lags
        self.corr = corr
        self.counts = counts
        self.norms = norms
        self.freq = Map.freq[freq1]

    def plot(self, finds1, finds2) :
        """After performing analysis, make plots."""

        markers = ['o','s','^','h','p','v']
        colours = ['b','g','r','c','m','y']
        nmarkers = len(markers)
        ncolours = len(colours)
        
        for find1 in finds1 :
            # Plot the correlations.
            plt.figure()
            for ii, find2 in enumerate(finds2) :
                plt.plot(self.lags, self.corr[find1,find2,:],
                         linestyle='None',
                         marker=markers[ii%nmarkers],
                         markerfacecolor=colours[ii%ncolours])
            plt.xlabel('lag (degrees)')
            plt.ylabel('normalized correlation with ' + 
                       str(int(round(self.freq[find1]/1.e6))) + ' MHz bin')
            #Plot the normalizations to give an absolute scale.
            plt.figure()
            for ii, find2 in enumerate(finds2) :
                plt.plot(self.freq[find2]/1.e6, self.norms[find1,find2], 
                         linestyle='None',
                         marker=markers[ii%nmarkers],
                         markerfacecolor=colours[ii%ncolours])
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('Normalization of correlation with ' +
                       str(int(round(self.freq[find1]/1.e6))) + 
                       ' MHz bin (K^2)')


        

# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        FreqSlices(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        FreqSlices().execute()
        


