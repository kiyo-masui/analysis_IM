
"""Program that calculates the correlation function across frequency slices.
"""

import copy

import scipy as sp
import numpy.ma as ma
from numpy import linalg
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
               'lags' : (0.1, 0.2)
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
        if len(params['file_middles']) == 1 :
            file_name = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end'])
            (Map1, Noise1) = fits_map.read(file_name)
            Map2 = copy.deepcopy(Map1)
            Noise2 = copy.deepcopy(Noise1)
        elif len(params['file_middles']) == 2 :
            file_name = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end'])
            (Map1, Noise1) = fits_map.read(file_name)
            file_name = (params['input_root'] + params['file_middles'][1] + 
                         params['input_end'])
            (Map2, Noise2) = fits_map.read(file_name)
        else :
            raise ce.FileParameterTypeError('For now can only process one'
                                            'file.')
        # If any eigenmodes have been marked for subtraction, subtract them
        # out.
        freq = sp.array(params['freq'], dtype=int)
        if hasattr(self, 'freq_eig_modes') :
            for ira in range(Map1.dims[0]) :
                for jdec in range(Map1.dims[1]) :
                    if sp.any(Map1.data.mask[ira,jdec,freq]) :
                        continue
                    else :
                        for v in self.freq_eig_modes :
                            v.shape = freq.shape
                            amp = sp.sum(v*Map1.data[ira,jdec,freq])
                            Map1.data[ira,jdec,freq] -= amp*v
            for ira in range(Map2.dims[0]) :
                for jdec in range(Map2.dims[1]) :
                    if sp.any(Map2.data.mask[ira,jdec,freq]) :
                        continue
                    else :
                        for v in self.freq_eig_modes :
                            v.shape = freq.shape
                            amp = sp.sum(v*Map2.data[ira,jdec,freq])
                            Map2.data[ira,jdec,freq] -= amp*v

        # Inverse noise weight the datapa.
        if (not sp.alltrue(Map1.data.mask == Noise1.data.mask) 
            or not sp.alltrue(Map2.data.mask == Noise2.data.mask)) :
            raise ce.DataError('Expected Map and Noise to have the same mask.')
        # Deweight noisy pixels by 1/(C_s+C_n).  Assume var(map)-mean(N) ~ C_s.
        var_f = (ma.var(Map1.data.reshape((Map1.dims[0]*Map1.dims[1],
                                         Map1.dims[2])), 0) 
                 - ma.mean(Noise1.data.reshape((Map1.dims[0]*Map1.dims[1],
                                              Map1.dims[2])), 0))
        Noise1.data += var_f
        var_f = (ma.var(Map2.data.reshape((Map2.dims[0]*Map2.dims[1],
                                         Map2.dims[2])), 0) 
                 - ma.mean(Noise2.data.reshape((Map2.dims[0]*Map2.dims[1],
                                              Map2.dims[2])), 0))
        Noise2.data += var_f
        Map1.data /= Noise1.data
        Map1.calc_axes()
        Map2.data /= Noise2.data
        Map2.calc_axes()
        # Set up correlation function.
        lags = sp.array(params['lags'])
        freq1 = freq
        freq2 = freq
        map1 = Map1.data[:,:,freq1]
        map2 = Map2.data[:,:,freq2]
        noise1 = Noise1.data[:,:,freq1]
        noise2 = Noise2.data[:,:,freq2]
        nlags = len(lags)
        nf = len(freq1)
        corr = sp.zeros((nf, nf, nlags), dtype=float)
        counts = sp.zeros(corr.shape, dtype=float)
        # Noting that if DEC != 0, then a degree of RA is less than a degree.
        ra_fact = sp.cos(sp.pi*Map1.lat[len(Map1.lat)//2]/180.0)
        
        # Double loop over all three coordinates. Using only a 2 space indent
        # because this loop is 6 deep.
        for ira1 in range(Map1.dims[0]) :
          for ira2 in range(Map2.dims[0]) :
            for jdec1 in range(Map1.dims[1]) :
              # This loop must start at zero to get all pairs.
              for jdec2 in range(Map2.dims[1]) :
                # any instead of all, making correlation <= 1.
                if (sp.any(map1.mask[ira1,jdec1,:]) or 
                    sp.any(map2.mask[ira2,jdec2,:])) :
                        continue
                # Figure out the lag then which bin it falls into.
                dra = (Map2.long[ira2] - Map1.long[ira1])*ra_fact
                ddec = Map2.lat[jdec2] - Map1.lat[jdec1]
                lag = sp.sqrt(dra**2 + ddec**2)
                lag_ind = 0
                for this_lag_ind, this_lag_bin in enumerate(lags) :
                    if lag > this_lag_bin :
                       lag_ind = this_lag_ind + 1
                    else :
                        break
                # Take the outer product accross frequencies.
                if lag_ind < nlags :
                    data1 = map1[ira1,jdec1,:]
                    data2 = map2[ira2,jdec2,:]
                    n1 = 1./noise1[ira1,jdec1,:] 
                    n2 = 1./noise2[ira2,jdec2,:] 
                    this_corr = (data1[:,sp.newaxis] * data2[sp.newaxis,:])
                    this_weight = n1[:,sp.newaxis] * n2[sp.newaxis,:]
                    corr[:,:,lag_ind] += this_corr.filled(0.)
                    counts[:,:,lag_ind] += this_weight
        corr /= counts
        norms = corr[:,:,0].diagonal()
        norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
        corr /= norms[:,:,sp.newaxis]
        corr = (corr + sp.swapaxes(corr, 0, 1))/2.0
        
        self.lags = lags
        self.corr = corr
        self.counts = counts
        self.norms = norms
        self.freq1 = Map1.freq[freq1]
        self.freq2 = Map2.freq[freq1]
        self.freq_inds = freq

    def get_freq_eigenmodes(self, n) :
        """After generating the covarience, find the n dominante frequency
        modes.  Mark them to be removed if execute is called again."""
        
        [h, V] = linalg.eig(self.corr[:,:,0]*self.norms)
        hs = list(h)
        hs.sort()
        vectors = []
        for ii in range(n) :
            ind, = sp.where(h == hs[-ii-1])
            vectors.append(V[:,ind])
        
        self.freq_eig_values = h
        self.freq_eig_modes = vectors

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

    def plot_eig(self) :
        """Plots the eigen values and prints out some statistics."""
        
        n = len(self.freq_inds)
        plt.semilogy(abs(sp.sort(self.freq_eig_values/n)), marker='o', 
                     linestyle='None')
        print 'Mean noise: ', sp.sum(self.freq_eig_values)/n
        print 'Largest eigenvalues/n : ', 
        print sp.sort(self.freq_eig_values/n)[-10:]

        

# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        FreqSlices(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        FreqSlices().execute()
        


