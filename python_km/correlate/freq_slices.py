
"""Program that calculates the correlation function across frequency slices.
"""

import copy

import scipy as sp
import numpy.ma as ma
from numpy import linalg
import matplotlib.pyplot as plt

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from core import fits_map, utils
import map.tools

params_init = {
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile",),
               'input_end_map' : "_map.fits",
               'input_end_noise' : "_noise.fits",
               'output_root' : "./testoutput",
               # What frequencies to correlate:
               'freq' : (),
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
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        # Get the map data from file.
        if len(params['file_middles']) == 1 :
            map_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_map'])
            noise_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_noise'])
            Map1 = fits_map.read(map_file)
            Noise1 = fits_map.read(noise_file)
            Map2 = copy.deepcopy(Map1)
            Noise2 = copy.deepcopy(Noise1)
        elif len(params['file_middles']) == 2 :
            map_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_map'])
            noise_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_noise'])
            Map1 = fits_map.read(map_file)
            Noise1 = fits_map.read(noise_file)
            map_file = (params['input_root'] + params['file_middles'][1] + 
                         params['input_end_map'])
            noise_file = (params['input_root'] + params['file_middles'][1] + 
                         params['input_end_noise'])
            Map2 = fits_map.read(map_file)
            Noise2 = fits_map.read(noise_file)
        else :
            raise ce.FileParameterTypeError('For now can only process one'
                                            'file.')
        # If any eigenmodes have been marked for subtraction, subtract them
        # out.
        freq = sp.array(params['freq'], dtype=int)
        subflag = False
        if (hasattr(self, 'freq_Lsvd_modes') and 
            hasattr(self, 'freq_Rsvd_modes')) :
            print 'Subtracting ' + str(len(self.freq_Lsvd_modes)),
            print 'frequency svd-modes.'
            Lmodes = self.freq_Lsvd_modes
            Rmodes = self.freq_Rsvd_modes
            subflag = True
        elif hasattr(self, 'freq_eig_modes') :
            print 'Subtracting ' + str(len(self.freq_eig_modes)),
            print 'frequency eigen-modes.'
            Lmodes = self.freq_eig_modes
            Rmodes = self.freq_eig_modes
            subflag = True
        if subflag :
            for ira in range(Map1.dims[0]) :
                for jdec in range(Map1.dims[1]) :
                    if sp.any(Map1.data.mask[ira,jdec,freq]) :
                        continue
                    else :
                        for v in Lmodes :
                            v.shape = freq.shape
                            amp = sp.sum(v*Map1.data[ira,jdec,freq])
                            Map1.data[ira,jdec,freq] -= amp*v
            map_out = (params['output_root'] + params['file_middles'][0] + 
                         '_cleaned' + params['input_end_map'])
            fits_map.write(Map1, map_out)
            for ira in range(Map2.dims[0]) :
                for jdec in range(Map2.dims[1]) :
                    if sp.any(Map2.data.mask[ira,jdec,freq]) :
                        continue
                    else :
                        for v in Rmodes :
                            v.shape = freq.shape
                            amp = sp.sum(v*Map2.data[ira,jdec,freq])
                            Map2.data[ira,jdec,freq] -= amp*v
            map_out = (params['output_root'] + params['file_middles'][1] + 
                         '_cleaned' + params['input_end_map'])
            fits_map.write(Map2, map_out)

        # Inverse noise weight the data.
        if (not sp.alltrue(Map1.data.mask == Noise1.data.mask) 
            or not sp.alltrue(Map2.data.mask == Noise2.data.mask)) :
            raise ce.DataError('Expected Map and Noise to have the same mask.')
        # Deweight noisy pixels by 1/(C_s+C_n).  Assume var(map)-mean(N) ~ C_s.
        #var_f = (ma.var(Map1.data.reshape((Map1.dims[0]*Map1.dims[1],
        #                                 Map1.dims[2])), 0) 
        #         - ma.mean(Noise1.data.reshape((Map1.dims[0]*Map1.dims[1],
        #                                      Map1.dims[2])), 0))
        #Noise1.data += var_f
        #var_f = (ma.var(Map2.data.reshape((Map2.dims[0]*Map2.dims[1],
        #                                 Map2.dims[2])), 0) 
        #         - ma.mean(Noise2.data.reshape((Map2.dims[0]*Map2.dims[1],
        #                                      Map2.dims[2])), 0))
        #Noise2.data += var_f
        Map1.data /= Noise1.data
        Map1.calc_axes()
        Map2.data /= Noise2.data
        Map2.calc_axes()
        # Set up correlation function.
        lags = sp.array(params['lags'])
        if len(lags)==0 :
            centre, shape, spacing = map.tools.get_map_params(Map1)
            # Assume the pixel width in real degrees is the dec spacing.
            wid = spacing[1]
            # Choose lags at fixed grid distance.
            lags = sp.array([0.0, 1.0, sp.sqrt(2.0), 2.0, sp.sqrt(5.0),
                             sp.sqrt(8.0), 3.0, sp.sqrt(10.0)])
            lags = (lags + 0.01)*wid
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
                    corr[:,:,lag_ind] += this_corr
                    counts[:,:,lag_ind] += this_weight
        corr /= counts
        norms = corr[:,:,0].diagonal()
        norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
        corr /= norms[:,:,sp.newaxis]
        
        self.lags = lags
        self.corr = corr
        self.counts = counts
        self.norms = norms
        self.freq1 = Map1.freq[freq1]
        self.freq2 = Map2.freq[freq1]
        self.allfreq = Map1.freq
        self.freq_inds = freq

    def get_freq_eigenmodes(self, n) :
        """After generating the covarience, find the n dominante frequency
        modes.  Mark them to be removed if execute is called again."""

        corr = (self.corr + sp.swapaxes(self.corr, 0, 1))/2.0
        [h, V] = linalg.eig(corr[:,:,0]*self.norms)
        hs = list(abs(h))
        hs.sort()
        vectors = []
        for ii in range(n) :
            ind, = sp.where(abs(h) == hs[-ii-1])
            if len(ind) > 1 :
                raise NotImplementedError('2 eigenvalues bitwise equal.')
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

    def plot_freq(self, norms=False) :
        
        nf = len(self.freq_inds)
        #freq_diffs = sp.sort(self.allfreq - min(self.allfreq))
        n_bins = 12
        factor = 1.5
        start = 2.1e6
        freq_diffs = sp.empty(n_bins)
        freq_diffs[0] = start
        for ii in range(1,n_bins) :
           freq_diffs[ii] = factor*freq_diffs[ii-1] 
        n_diffs = len(freq_diffs)
        corrf = sp.zeros(n_diffs)
        countsf = sp.zeros(n_diffs, dtype=int)
        for ii in range(nf) :
            for jj in range(nf) :
                if norms :
                    thiscorr = self.corr[ii,jj,0]*self.norms[ii,jj]
                else :
                    thiscorr = self.corr[ii,jj,0]
                df = abs(self.freq1[ii] - self.freq2[jj])
                for kk in range(1, n_diffs-1) :
                    if (abs(freq_diffs[kk]-df) <= abs(freq_diffs[kk-1]-df)
                        and abs(freq_diffs[kk]-df) < abs(freq_diffs[kk+1]-df)) :
                        d_ind = kk
                if abs(freq_diffs[0]-df) < abs(freq_diffs[1]-df) :
                    d_ind = 0
                if abs(freq_diffs[-1]-df) <= abs(freq_diffs[-2]-df) :
                    d_ind = n_diffs - 1
                corrf[d_ind] += thiscorr
                countsf[d_ind] +=1
        corrf /= countsf
        plt.semilogx(freq_diffs/1e6, corrf*1e6, linestyle='None',
                 marker='o', markerfacecolor='b')
        plt.xlabel('Frequency Lag (MHz)')
        plt.ylabel('Correlation (mK$^2$)')

    def get_freq_svd_modes(self, n) :
        """Same as get freq eigenmodes, but treats left and right maps
        separatly with an SVD.
        """
        U, s, V = linalg.svd(self.corr[:,:,0]*self.norms)
        V = V.T
        hs = list(s)
        hs.sort()
        Lvectors = []
        Rvectors = []
        for ii in range(n) :
            ind, = sp.where(abs(s) == hs[-ii-1])
            if len(ind) > 1 :
                raise NotImplementedError('2 eigenvalues bitwise equal.')
            Lvectors.append(U[:,ind])
            Rvectors.append(V[:,ind])
        
        self.freq_svd_values = s
        self.freq_Lsvd_modes = Lvectors
        self.freq_Rsvd_modes = Rvectors

    def plot_svd(self) :
        """Plots the eigen values and prints out some statistics."""
        
        n = len(self.freq_inds)
        plt.semilogy(abs(sp.sort(self.freq_svd_values/n)), marker='o', 
                     linestyle='None')
        print 'Mean noise: ', sp.sum(self.freq_svd_values)/n
        print 'Largest eigenvalues/n : ', 
        print sp.sort(self.freq_svd_values/n)[-10:]

def degrade_corr_res(self, beam=0.4) :
    """Brings all freqencies to the same resolution."""

    if self.params['lags'] != () :
        raise ValueError("Expeceted automatically genreated lags.")
    wid = self.lags[0]/1.01
    real_lags = (self.lags/wid -0.01)*wid
    lag_weights = sp.array([1.0, 4.0, 4.0, 4.0, 8.0, 4.0, 4.0, 8.0])
    # Array of the square gid distances of each of the lags.
    sq_ind_lag = sp.array([0, 1, 2, 4, 5, 8, 9, 10])
    sig = (beam/2.2355)**2
    sig1 = sig - (utils.get_beam(self.freq1)/2.355)**2
    sig2 = sig - (utils.get_beam(self.freq2)/2.355)**2
    if min(sig1) < 0 or min(sig2) < 0 :
        raise ValueError("Final beam must be bigger than lowest resolution.")
    n1 = len(self.freq1)
    n2 = len(self.freq2)
    new_corr = sp.zeros((n1,n2,1))
    corr = self.corr*self.norms[:,:,sp.newaxis]

    # Make an array of all lag indicies squared.
    one_d_inds = sp.arange(-3, 3+1, dtype=int)
    n = len(oned)
    inds_sq_1 = sp.resize(sp.reshape(one_d_inds, (n,1,1,1))**2 
                          + sp.reshape(one_d_inds, (1,n,1,1))**2,
                          (n,n,n,n)).flatten()
    inds_sq_2 = sp.resize(sp.reshape(one_d_inds, (1,1,n,1))**2 
                          + sp.reshape(one_d_inds, (1,1,1,n))**2,
                          (n,n,n,n)).flatten()
    inds_dif_sq = ((sp.reshape(one_d_inds,(n,1,1,1)) 
                    - sp.reshape(one_d_inds,(1,n,1,1)))**2
                   +  (sp.reshape(one_d_inds,(1,1,n,1)) 
                       - sp.reshape(one_d_inds,(1,1,1,n)))**2).flatten()
    max_lag_sq = max(sq_ind_lag)
    in_bounds = sp.logical_and(sp.logical_and(inds_sq_1 <= max_lag_sq,
                    inds_sq_1 <= max_lag_sq), inds_dif_sq <= max_lag_sq)
    inds_sq_1 = inds_sq_1[in_bounds]
    inds_sq_2 = inds_sq_2[in_bounds]
    inds_dif_sq = inds_dif_sq[in_bounds]
    for ii in range(n1) :
        for jj in range(n2) :
            
            this_norm = 0
            s1 = sig1[ii]
            s2 = sig2[jj]
            b1 = sp.exp(-real_lags**2/2/s1)/(2*sp.pi*s1)
            b2 = sp.exp(-real_lags**2/2/s2)/(2*sp.pi*s2)
            for kk in range(len(real_lags)) :
                new_corr[ii,jj,0] += b1*b2*lag_weights[kk]*corr[ii,jj,kk]
                this_norm += b1*b2*lag_weights[kk]
            new_corr[ii,jj,0] /=  this_norm * (sp.sqrt(s1/s2) + sp.sqrt(s1/s2))
    corr = new_corr
    norms = corr[:,:,0].diagonal()
    norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
    corr /= norms[:,:,sp.newaxis]
    self.norms = norms
    self.corr = corr


def subtract_svd_modes_corr(self) :
    """Removes svd modes from corr (as opposed to the map)."""

    corr = self.corr*self.norms[:,:,sp.newaxis]
    for ii in range(len(self.freq_Lsvd_modes)):
        corr = corr - (self.freq_svd_values[ii] 
                       * self.freq_Lsvd_modes[ii][:,sp.newaxis]
                       * self.freq_Rsvd_modes[ii][sp.newaxis,:])
    norms = corr[:,:,0].diagonal()
    norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
    corr /= norms[:,:,sp.newaxis]
    self.norms = norms
    self.corr = corr


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        FreqSlices(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        FreqSlices().execute()
        


