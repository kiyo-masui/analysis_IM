
"""Program that calculates the correlation function across frequency slices.
"""

import copy

import scipy as sp
import numpy.ma as ma
from numpy import linalg
from numpy import random
import matplotlib.pyplot as plt

from kiyopy import parse_ini
import kiyopy.utils
import kiyopy.custom_exceptions as ce
from core import utils, algebra
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
        if len(params['file_middles']) == 1 :
            fn = params['file_middles'][0]
            params['file_middles'] = (fn, fn)
        # Get the map data from file.
        #if len(params['file_middles']) == 1 :
        #    map_file = (params['input_root'] + params['file_middles'][0] + 
        #                 params['input_end_map'])
        #    noise_file = (params['input_root'] + params['file_middles'][0] + 
        #                 params['input_end_noise'])
        #    Map1 = algebra.load(map_file)
        #    print "Loading noise."
        #    Noise1 = algebra.make_mat(algebra.open_memmap(noise_file, mode='r'))
        #    Noise1 = Noise1.mat_diag()
        #    print "Done."
        #    Map2 = copy.deepcopy(Map1)
        #    Noise2 = copy.deepcopy(Noise1)
        if len(params['file_middles']) == 2 :
            map_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_map'])
            noise_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_noise'])
            Map1 = algebra.make_vect(algebra.load(map_file))
            print "Loading noise."
            Noise1 = algebra.make_mat(algebra.open_memmap(noise_file, mode='r'))
            Noise1 = Noise1.mat_diag()
            print "Done."
            map_file = (params['input_root'] + params['file_middles'][1] + 
                         params['input_end_map'])
            noise_file = (params['input_root'] + params['file_middles'][1] + 
                         params['input_end_noise'])
            Map2 = algebra.make_vect(algebra.load(map_file))
            print "Loading noise."
            Noise2 = algebra.make_mat(algebra.open_memmap(noise_file, mode='r'))
            Noise2 = Noise2.mat_diag()
            print "Done."
        else :
            raise ce.FileParameterTypeError('For now can only process one'
                                            ' or two files.')
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
            for ira in range(Map1.shape[1]) :
                for jdec in range(Map1.shape[2]) :
                    # if sp.any(Map1.data.mask[ira,jdec,freq]) :
                    #    continue
                    # else :
                    for v in Lmodes :
                        # v.shape = freq.shape
                        v = v.reshape(freq.shape)
                        # amp = sp.sum(v*Map1.data[ira,jdec,freq])
                        amp = sp.dot(v, Map1[freq,ira,jdec])
                        Map1[freq,ira,jdec] -= amp*v
                        map_out = (params['output_root'] + params['file_middles'][0] + 
                                   '_cleaned' + params['input_end_map'])
            # fits_map.write(Map1, map_out)
            algebra.save(map_out, Map1)
            
            for ira in range(Map2.shape[1]) :
                for jdec in range(Map2.shape[2]) :
                    # if sp.any(Map1.data.mask[ira,jdec,freq]) :
                    #    continue
                    # else :
                    for v in Lmodes :
                        # v.shape = freq.shape
                        v = v.reshape(freq.shape)
                        # amp = sp.sum(v*Map1.data[ira,jdec,freq])
                        amp = sp.dot(v, Map2[freq,ira,jdec])
                        Map2[freq,ira,jdec] -= amp*v
                        map_out = (params['output_root'] + params['file_middles'][1] + 
                                   '_cleaned' + params['input_end_map'])
            # fits_map.write(Map2, map_out)
            algebra.save(map_out, Map2)

        # Inverse noise weight the data.
        #if (not sp.alltrue(Map1.data.mask == Noise1.data.mask) 
        #    or not sp.alltrue(Map2.data.mask == Noise2.data.mask)) :
        #    raise ce.DataError('Expected Map and Noise to have the same mask.')
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
        #map1.calc_axes()
        #Map2.calc_axes()
        # Set up correlation function.
        lags = sp.array(params['lags'])
        if len(lags)==0 :
            # Assume the pixel width in real degrees is the dec spacing.
            wid = abs(map1.info['dec_delta'])
            # Choose lags at fixed grid distance.
            lags = sp.array([0.0, 1.0, sp.sqrt(2.0), 2.0, sp.sqrt(5.0),
                             sp.sqrt(8.0), 3.0, sp.sqrt(10.0)])
            lags = (lags + 0.01)*wid
        freq1 = freq
        freq2 = freq

        f1 = Map1.get_axis('freq')
        f2 = Map2.get_axis('freq')
        r1 = Map1.get_axis('ra')
        r2 = Map2.get_axis('ra')
        d1 = Map1.get_axis('dec')
        d2 = Map2.get_axis('dec')

        map1 = Map1[freq1,:,:]
        map2 = Map2[freq2,:,:]
        noise1 = Noise1[freq1,:,:]
        noise2 = Noise2[freq2,:,:]

        # Noise weight
        map1 *= noise1
        map2 *= noise2

        nlags = len(lags)
        nf = len(freq1)
        corr = sp.zeros((nf, nf, nlags), dtype=float)
        counts = sp.zeros(corr.shape, dtype=float)
        # Noting that if DEC != 0, then a degree of RA is less than a degree.
        ra_fact = sp.cos(sp.pi*Map1.info['dec_centre'] / 180.0) * Map1.info['ra_delta']
        #dec_fact = Map1.info['dec_delta']

        # Double loop over all three coordinates. Using only a 2 space indent
        # because this loop is 6 deep.
        print "Starting Correlation."
        for ira1 in range(map1.shape[1]) :
            for ira2 in range(map2.shape[1]) :
                for jdec1 in range(map1.shape[2]) :
                    # This loop must start at zero to get all pairs.
                    for jdec2 in range(map2.shape[2]) :
                        # any instead of all, making correlation <= 1.
                        #if (sp.any(map1.mask[ira1,jdec1,:]) or 
                        #    sp.any(map2.mask[ira2,jdec2,:])) :
                        #continue
                        # Figure out the lag then which bin it falls into.
                        dra = (r2[ira2] - r1[ira1])* ra_fact
                        ddec = (d2[jdec2] - d1[jdec1])
                        lag = sp.sqrt(dra**2 + ddec**2)
                        lag_ind = 0
                        for this_lag_ind, this_lag_bin in enumerate(lags) :
                            if lag > this_lag_bin :
                                lag_ind = this_lag_ind + 1
                            else :
                                break
                        # Take the outer product accross frequencies.
                        if lag_ind < nlags :
                            data1 = map1[:,ira1,jdec1]
                            data2 = map2[:,ira2,jdec2]
                            n1 = noise1[:,ira1,jdec1] 
                            n2 = noise2[:,ira2,jdec2] 
                            this_corr = (data1[:,sp.newaxis] * data2[sp.newaxis,:])
                            this_weight = n1[:,sp.newaxis] * n2[sp.newaxis,:]
                            corr[:,:,lag_ind] += this_corr
                            counts[:,:,lag_ind] += this_weight
        corr /= counts
        norms = corr[:,:,0].diagonal()
        if sp.any(norms<0) and not subflag :
            raise RunTimeError("Something went horribly wrong")
        norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
        corr /= norms[:,:,sp.newaxis]
        
        self.lags = lags
        self.corr = corr
        self.counts = counts
        self.norms = norms
        self.freq1 = f1[freq1]
        self.freq2 = f2[freq2]
        self.allfreq = f1
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
            vectors.append(V[:,ind[0]])
        
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
                       str(int(round(self.freq1[find1]/1.e6))) + ' MHz bin')
            #Plot the normalizations to give an absolute scale.
            plt.figure()
            for ii, find2 in enumerate(finds2) :
                plt.plot(self.freq2[find2]/1.e6, self.norms[find1,find2], 
                         linestyle='None',
                         marker=markers[ii%nmarkers],
                         markerfacecolor=colours[ii%ncolours])
            plt.xlabel('Frequency (MHz)')
            plt.ylabel('Normalization of correlation with ' +
                       str(int(round(self.freq1[find1]/1.e6))) + 
                       ' MHz bin (K^2)')

    def plot_eig(self) :
        """Plots the eigen values and prints out some statistics."""
        
        n = len(self.freq_inds)
        plt.semilogy(abs(sp.sort(self.freq_eig_values/n)), marker='o', 
                     linestyle='None')
        print 'Mean noise: ', sp.sum(self.freq_eig_values)/n
        print 'Largest eigenvalues/n : ', 
        print sp.sort(self.freq_eig_values/n)[-10:]

    def plot_freq(self, norms=False, lag_ind=0) :

        # Set up binning.    
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
        # Allowcate memory.
        corrf = sp.zeros(n_diffs)
        countsf = sp.zeros(n_diffs, dtype=int)
        for ii in range(nf) :
            for jj in range(nf) :
                if norms :
                    thiscorr = self.corr[ii,jj,lag_ind]*self.norms[ii,jj]
                else :
                    thiscorr = self.corr[ii,jj,lag_ind]
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
                 marker='o')
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
            Lvectors.append(U[:,ind[0]])
            Rvectors.append(V[:,ind[0]])
        
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

def plot_contour(self, norms=False, lag_inds=(0)) :
    
    lag_inds = list(lag_inds)
    # Set up binning.    
    nf = len(self.freq_inds)
    #freq_diffs = sp.sort(self.allfreq - min(self.allfreq))
    n_bins = 12
    factor = 2.0
    start = 2.1e6
    freq_diffs = sp.empty(n_bins)
    freq_diffs[0] = 0.5
    freq_diffs[1] = start
    for ii in range(2,n_bins) :
       freq_diffs[ii] = factor*freq_diffs[ii-1]
    n_diffs = len(freq_diffs)
    # Allowcate memory.
    corrf = sp.zeros((n_diffs, len(lag_inds)))
    countsf = sp.zeros(n_diffs, dtype=int)
    for ii in range(nf) :
        for jj in range(nf) :
            if norms :
                thiscorr = (self.corr[ii,jj,lag_inds] * 
                            self.norms[ii,jj,sp.newaxis])
            else :
                thiscorr = self.corr[ii,jj,lag_inds]
            df = abs(self.freq1[ii] - self.freq2[jj])
            for kk in range(1, n_diffs-1) :
                if (abs(freq_diffs[kk]-df) <= abs(freq_diffs[kk-1]-df)
                    and abs(freq_diffs[kk]-df) < abs(freq_diffs[kk+1]-df)) :
                    d_ind = kk
            if abs(freq_diffs[0]-df) < abs(freq_diffs[1]-df) :
                d_ind = 0
            if abs(freq_diffs[-1]-df) <= abs(freq_diffs[-2]-df) :
                d_ind = n_diffs - 1
            corrf[d_ind, :] += thiscorr
            countsf[d_ind] +=1
    corrf /= countsf[:, sp.newaxis]
    pdat = sp.sign(corrf)*sp.sqrt(abs(corrf))*1e3
    a = plt.figure()
    #a.set_figwidth(a.get_figwidth()/3.0)
    f = plt.contourf(self.lags[lag_inds], (freq_diffs)/1e6, pdat)
    f.ax.set_xscale('log')
    f.ax.set_yscale('log')
    plt.axis('scaled')
    plt.xlim((0.05, 0.9))
    plt.ylim((0.8, 100))
    plt.xlabel("angular lag, $\sigma$ (degrees, 34$\cdotp$Mpc/h)")
    plt.ylabel("frequency lag, $\pi$ (MHz, 4.5$\cdotp$Mpc/h)")
    c = plt.colorbar(f)
    c.ax.set_ylabel("correlation (mK)")

def plot_collapsed(self, norms=False, lag_inds=(0), save_old=False,
                   plot_old=False) :
    
    lag_inds = list(lag_inds)
    # Set up binning.    
    nf = len(self.freq_inds)
    freq_diffs = sp.arange(0.1e6, 100e6, 2e6)
    n_diffs = len(freq_diffs)
    # Allowcate memory.
    corrf = sp.zeros((n_diffs, len(lag_inds)))
    countsf = sp.zeros(n_diffs, dtype=int)
    for ii in range(nf) :
        for jj in range(nf) :
            if norms :
                thiscorr = (self.corr[ii,jj,lag_inds] * 
                            self.norms[ii,jj,sp.newaxis])
            else :
                thiscorr = self.corr[ii,jj,lag_inds]
            df = abs(self.freq1[ii] - self.freq2[jj])
            for kk in range(1, n_diffs-1) :
                if (abs(freq_diffs[kk]-df) <= abs(freq_diffs[kk-1]-df)
                    and abs(freq_diffs[kk]-df) < abs(freq_diffs[kk+1]-df)) :
                    d_ind = kk
            if abs(freq_diffs[0]-df) < abs(freq_diffs[1]-df) :
                d_ind = 0
            if abs(freq_diffs[-1]-df) <= abs(freq_diffs[-2]-df) :
                d_ind = n_diffs - 1
            corrf[d_ind, :] += thiscorr
            countsf[d_ind] +=1
    corrf /= countsf[:, sp.newaxis]
    # Now collapse to 1 axis:
    a_fact = 34.0
    f_fact = 4.5
    nbins = 8
    lags = sp.empty(nbins)
    lags[0] = 6.0
    lags[1] = 10.0
    for ii in range(2, nbins) :
        lags[ii] = 1.5*lags[ii-1]
    R = self.lags[lag_inds]
    R = (a_fact*R[:, sp.newaxis])**2
    R = R + (f_fact*freq_diffs[sp.newaxis, :]/1.0e6)**2
    R = sp.sqrt(R)
    corr = sp.zeros(nbins)
    counts = sp.zeros(nbins, dtype=int)
    lag_weights = sp.arange(1, len(lag_inds) + 1)**2
    for ii in range(len(lag_inds)) :
        for jj in range(n_diffs) :
            dR = R[ii, jj]
            for kk in range(1, nbins - 1) :
                if (abs(lags[kk]-dR) <= abs(lags[kk-1]-dR)
                    and abs(lags[kk]-dR) < abs(lags[kk+1]-dR)) :
                    ind = kk
            if abs(lags[0]-dR) < abs(lags[1]-dR) :
                ind = 0
            if abs(lags[-1]-dR) <= abs(lags[-2]-dR) :
                ind = nbins - 1
            counts[ind] += lag_weights[ii]
            corr[ind] += lag_weights[ii]*corrf[jj, ii]
    n_boot = 1000
    corrb = ma.zeros((nbins, n_boot))
    countsb = ma.zeros((nbins, n_boot), dtype=int)
    for qq in range(n_boot) :
        for mm in range(n_diffs*len(lag_inds)) :
            ii = random.random_integers(0, len(lag_inds)-1)
            jj = random.random_integers(0, n_diffs-1)
            dR = R[ii, jj]
            for kk in range(1, nbins - 1) :
                if (abs(lags[kk]-dR) <= abs(lags[kk-1]-dR)
                    and abs(lags[kk]-dR) < abs(lags[kk+1]-dR)) :
                    ind = kk
            if abs(lags[0]-dR) < abs(lags[1]-dR) :
                ind = 0
            if abs(lags[-1]-dR) <= abs(lags[-2]-dR) :
                ind = nbins - 1
            countsb[ind, qq] += lag_weights[ii]
            corrb[ind, qq] += lag_weights[ii]*corrf[jj, ii]
    corr = corr/counts
    corrb = corrb/countsb
    pdat = sp.sign(corr)*sp.sqrt(abs(corr))*1e3
    pdatb = sp.sign(corrb)*sp.sqrt(abs(corrb))*1e3
    pdatb_mean = sp.mean(pdatb, -1)
    pdatb_sig = sp.std(pdatb, -1)
    a = plt.figure()
    ax = plt.gca()
    ax.set_yscale("log")
    ax.set_xscale("log")
    errors = sp.empty((2, nbins))
    errors[:,:] = pdatb_sig
    if save_old :
        self.old_pdat = pdat
        self.old_errors = errors
    elin = 2
    msize = 6
    f = plt.errorbar(lags[pdat-errors[0,:]>0], pdat[pdat-errors[0,:]>0],
                     errors[:,pdat-errors[0,:]>0], linestyle='None', marker='o', 
                     color='b', elinewidth=elin, markersize=msize)
    f = plt.errorbar(lags[pdat+errors[1,:]<0], -pdat[pdat+errors[1,:]<0],
                     errors[:,pdat+errors[1,:]<0], linestyle='None', marker='o', 
                     color='r', elinewidth=elin, markersize=msize)
    inds = sp.logical_and(pdat-errors[0,:]<=0, pdat > 0)
    if sp.any(inds) :
        vals = pdat[inds] + 2*errors[1, inds]
        es = sp.zeros((2, len(vals)))
        es[0,:] = 0.25*abs(vals)
        f = plt.errorbar(lags[inds], vals,
                         es, linestyle='None', marker='None', 
                         color='b', lolims=True, elinewidth=elin, 
                         markersize=msize)
    inds = sp.logical_and(pdat+errors[1,:]>=0, pdat < 0)
    if sp.any(inds) :
        vals = pdat[inds] - 2*errors[0, inds]
        es = sp.zeros((2, len(vals)))
        es[0,:] = 0.25*abs(vals)
        f = plt.errorbar(lags[inds], -vals,
                         es, linestyle='None', marker='None', 
                         color='r', lolims=True, elinewidth=elin,
                         markersize=msize)
    t_lags = sp.arange(0.1,100,0.1)
    r0 = 5.5
    rb = 7.0
    t = (sp.sqrt(((rb + t_lags)/r0)**(-1.8)))
    t = t*0.15/t[0]
    f = plt.plot(t_lags, t, marker='None', color='k', linestyle='-')
    if plot_old :
        elin = 0.4
        msize = 6
        mfc = 'w'
        pdat = self.old_pdat
        errors = self.old_errors
        f = plt.errorbar(lags[pdat-errors[0,:]>0], pdat[pdat-errors[0,:]>0],
                         errors[:,pdat-errors[0,:]>0], linestyle='None', 
                         marker='o', color='b', elinewidth=elin, mfc=mfc, 
                         markersize=msize)
        f = plt.errorbar(lags[pdat+errors[1,:]<0], -pdat[pdat+errors[1,:]<0],
                         errors[:,pdat+errors[1,:]<0], linestyle='None',
                         marker='o', color='r', elinewidth=elin, mfc=mfc,
                         markersize=msize)
        inds = sp.logical_and(pdat-errors[0,:]<=0, pdat > 0)
        if sp.any(inds) :
            vals = pdat[inds] + 2*errors[1, inds]
            es = sp.zeros((2, len(vals)))
            es[0,:] = 0.25*abs(vals)
            f = plt.errorbar(lags[inds], vals,
                             es, linestyle='None', marker='None', 
                             color='b', lolims=True, elinewidth=elin,
                             markersize=msize, mfc=mfc)
        inds = sp.logical_and(pdat+errors[1,:]>=0, pdat < 0)
        if sp.any(inds) :
            vals = pdat[inds] - 2*errors[0, inds]
            es = sp.zeros((2, len(vals)))
            es[0,:] = 0.25*abs(vals)
            f = plt.errorbar(lags[inds], -vals,
                             es, linestyle='None', marker='None', 
                             color='r', lolims=True, elinewidth=elin,
                             markersize=msize, mfc=mfc)
    plt.axis([4, 100, 0.01, 200.0])
    plt.xlabel('lag (Mpc/h)')
    plt.ylabel('correlation (mK)')

def degrade_corr_res(self, beam=0.4) :
    """Brings all freqencies to the same resolution in correlation."""

    # XXX: Try this instead of inteegrating, us the analytic factors for if the
    # intrinsic cerrelation funciton is a delta.
    if self.params['lags'] != () or self.corr.shape[-1] != len(self.lags):
        raise ValueError("Expeceted automatically generated lags.")
    wid = self.lags[1]/1.01
    real_lags = (self.lags/wid - 0.01)*wid
    lag_weights = sp.array([1.0, 4.0, 4.0, 4.0, 8.0, 4.0, 4.0, 8.0])*wid**2
    # Array of the square gid distances of each of the lags.
    #sq_ind_lag = sp.array([0, 1, 2, 4, 5, 8, 9, 10])
    sig = (beam/2.355)**2
    sig1 = sig - (utils.get_beam(self.freq1)/2.355)**2
    sig2 = sig - (utils.get_beam(self.freq2)/2.355)**2
    if min(sig1) < 0 or min(sig2) < 0 :
        raise ValueError("Final beam must be bigger than lowest resolution.")
    n1 = len(self.freq1)
    n2 = len(self.freq2)
    new_corr = sp.zeros((n1,n2,1))
    corr = self.corr*self.norms[:,:,sp.newaxis]

    for ii in range(n1) :
        for jj in range(n2) :
            s1 = sig1[ii]
            s2 = sig2[jj]
            window = sp.exp(-real_lags**2/(2*(s1+s2)))/(2*sp.pi*(s1+s2))
            new_corr[ii,jj,0] = sp.sum(window*lag_weights*corr[ii,jj,:])
            norm = sp.sum(window*lag_weights)
            new_corr[ii,jj,0] /= norm
    corr = new_corr
    norms = corr[:,:,0].diagonal()
    norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
    corr /= norms[:,:,sp.newaxis]
    self.norms = norms
    self.corr = corr

def fake_degrade_corr_res(self, beam=0.4) :
    """Brings all frequencies to the same resolution by scaling the 0 lag
    correation as if it was a delta function."""

    sig = (beam/2.355)**2
    sig1 = (utils.get_beam(self.freq1)/2.355)**2
    sig2 = (utils.get_beam(self.freq2)/2.355)**2
    if min(sig1) < 0 or min(sig2) < 0 :
        raise ValueError("Final beam must be bigger than lowest resolution.")
    n1 = len(self.freq1)
    n2 = len(self.freq2)
    new_corr = sp.zeros((n1,n2,1))
    corr = self.corr*self.norms[:,:,sp.newaxis]

    for ii in range(n1) :
        for jj in range(n2) :
            s1 = sig1[ii]
            s2 = sig2[jj]
            new_corr[ii,jj,0] = corr[ii,jj,0]*(s1+s2)/(2*sig)

    corr = new_corr
    norms = corr[:,:,0].diagonal()
    norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
    corr /= norms[:,:,sp.newaxis]
    self.norms = norms
    self.corr = corr

def subtract_svd_modes_corr(self) :
    """Removes svd modes from corr (as opposed to the map)."""

    corr = self.corr*self.norms[:,:,sp.newaxis]
    for ii in range(len(self.freq_Lsvd_modes)) :
        corr = corr - (self.freq_svd_values[ii] 
                       * self.freq_Lsvd_modes[ii][:,sp.newaxis,sp.newaxis]
                       * self.freq_Rsvd_modes[ii][sp.newaxis,:,sp.newaxis])
    norms = corr[:,:,0].diagonal()
    norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
    corr /= norms[:,:,sp.newaxis]
    self.norms = norms
    self.corr = corr

def generate_fake_corr(self) :
    """generate a fake correlation matrix."""

    if self.params['lags'] != () :
        raise ValueError("Expeceted automatically genreated lags.")
    wid = self.lags[1]/1.01
    real_lags = (self.lags/wid - 0.01)*wid
    sig1 = (utils.get_beam(self.freq1)/2.355)**2
    sig2 = (utils.get_beam(self.freq2)/2.355)**2
    if min(sig1) < 0 or min(sig2) < 0 :
        raise ValueError("Final beam must be bigger than lowest resolution.")
    n1 = len(self.freq1)
    n2 = len(self.freq2)
    n_lags = len(real_lags)
    new_corr = sp.zeros((n1,n2,n_lags))

    for ii in range(n1) :
        for jj in range(n2) :
            s1 = sig1[ii]
            s2 = sig2[jj]
            window = sp.exp(-real_lags**2/(2*(s1+s2)))/(2*sp.pi*(s1+s2))
            new_corr[ii,jj,:] = window
    norms = new_corr[:,:,0].diagonal()
    norms = sp.sqrt(norms[:,sp.newaxis]*norms[sp.newaxis,:])
    new_corr /= norms[:,:,sp.newaxis]
    self.norms *= norms
    self.corr = new_corr


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        FreqSlices(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        FreqSlices().execute()
        


