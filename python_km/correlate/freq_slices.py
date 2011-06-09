
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
from map import beam

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
               'lags' : (0.1, 0.2),
               'convolve' : False,
               'sub_weighted_mean' : False,
               'modes' : 10,
               'first_pass_only' : False,
               'make_plots' : False,
               'factorizable_noise' : False
               }
prefix = 'fs_'

class MapPair(object) :
    """Pair of maps that are processed together and cross correlated.

    Parameters
    ----------
    Map1
    Map2
    Noise_inv1
    Noise_inv2
    freq
    """

    def __init__(self, Map1, Map2, Noise_inv1, Noise_inv2, freq) :
        self.Map1 = Map1
        self.Map2 = Map2
        self.Noise_inv1 = Noise_inv1
        self.Noise_inv2 = Noise_inv2
        self.freq = freq

    def degrade_resolution(self) :
        """Convolves the maps down to the lowest resolution.

        Also convolves the noise, making sure to deweight pixels near the edge
        as well.  Converts noise to factorizable form by averaging.
        """
        Noise1 = self.Noise_inv1
        Noise2 = self.Noise_inv2

        # Get the beam data.
        gfreq=self.Map1.get_axis("freq")
        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792, 
                 0.281176247549, 0.270856788455, 0.26745856078, 
                 0.258910010848, 0.249188429031])
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
        freq_data *= 1.0e6
        beam_diff=sp.sqrt(max(1.1*beam_data)**2-(beam_data)**2)
        b = beam.GaussianBeam(beam_diff,freq_data)
        # Convolve to a common resolution.
        self.Map2=b.apply(self.Map2)
        self.Map1=b.apply(self.Map1)
        # Reduce noise to factorizable.
        
        
        # This block of code needs to be split off into a function and applied
        # twice (so we are sure to do the same thing to each).
        Noise1[Noise1<1.e-30]=1.e-30
        Noise1=1./Noise1
        Noise1=b.apply(Noise1,cval=1.e30)

        Noise1=ma.array(Noise1)
        Noise1f=sp.mean(Noise1,0)
        Noise1[Noise1>1.e20]=ma.masked
        Noise1=Noise1/Noise1f[None,:,:]
        Noise1=sp.mean(sp.mean(Noise1,1),1)[:,None,None]*Noise1f[None,:,:]
        Noise1=(1/Noise1).filled(0)

        Noise2[Noise2<1.e-30]=1.e-30
        Noise2=1/Noise2
        Noise2=b.apply(Noise2,cval=1.e30)
        Noise2=ma.array(Noise2)
        Noise2f=sp.mean(Noise2,0)
        Noise2[Noise2>1.e20]=ma.masked
        Noise2=Noise2/Noise2f[None,:,:]
        Noise2=sp.mean(sp.mean(Noise2,1),1)[:,None,None]*Noise2f[None,:,:]
        Noise2=(1/Noise2).filled(0)

        self.Noise_inv1 = Noise1
        self.Noise_inv2 = Noise2

    def make_noise_factorizable(self) :
        # TODO.
        pass

    def subtract_weighted_mean(self) :
        """Subtracts the weighted mean from each frequency slice."""
        means1 = sp.sum(sp.sum(self.Noise_inv1*self.Map1, -1), -1)
        means1 /= sp.sum(sp.sum(self.Noise_inv1, -1), -1)
        means1.shape += (1, 1)
        self.Map1 -= means1
        means2 = sp.sum(sp.sum(self.Noise_inv2*self.Map2, -1), -1)
        means2 /= sp.sum(sp.sum(self.Noise_inv2, -1), -1)
        means2.shape += (1, 1)
        self.Map2 -= means2

    def subtract_frequency_modes(self, modes1, modes2=None, fname1=None, 
                                 fname2=None) :
        """Subtract frequency mode from the map.

        Parameters
        ---------
        modes1 : list of 1D arrays.
            Arrays must be the same length as self.freq.  Modes to subtract out
            of the map one.
        modes2 : list of 1D arrays.
            Modes to subtract out of map 2.  If `None` set to `modes1`.
        fname1 : string
            File name to save the mode amplitude maps too for map1.  If None,
            don't save.
        fname2 : string
            File name to save the mode amplitude maps too for map2.  If None,
            don't save.
        """

        if modes2 == None :
            modes2 = modes1

        Map1 = self.Map1
        Map2 = self.Map2
        freq = self.freq
        
        # First map.
        if fname1 :
            outmap = sp.empty((len(modes1),)+Map1.shape[1:])
            outmap = algebra.make_vect(outmap,axis_names=('freq', 'ra', 'dec'))
            outmap.copy_axis_info(Map1)
        for ira in range(Map1.shape[1]) :
            for jdec in range(Map1.shape[2]) :
                # if sp.any(Map1.data.mask[ira,jdec,freq]) :
                #    continue
                # else :
                for i,v in enumerate(modes1) :
                    # v.shape = freq.shape
                    v = v.reshape(freq.shape)
                    # amp = sp.sum(v*Map1.data[ira,jdec,freq])
                    amp = sp.dot(v, Map1[freq,ira,jdec])
                    Map1[freq,ira,jdec] -= amp*v
                    if fname1 :
                        outmap[i,ira,jdec] = amp
        if fname1 :
            algebra.save(fname1, outmap)
        
        # Second map.
        if fname2 :
            outmap = sp.empty((len(modes1),)+Map1.shape[1:])
            outmap = algebra.make_vect(outmap,axis_names=('freq', 'ra', 'dec'))
            outmap.copy_axis_info(Map2)
        for ira in range(Map2.shape[1]) :
            for jdec in range(Map2.shape[2]) :
                # if sp.any(Map1.data.mask[ira,jdec,freq]) :
                #    continue
                # else :
                for v in modes1 :
                    # v.shape = freq.shape
                    v = v.reshape(freq.shape)
                    amp = sp.dot(v, Map2[freq,ira,jdec])
                    Map2[freq,ira,jdec] -= amp*v
                    if fname2 :
                        outmap[i,ira,jdec] = amp
        if fname2 :
            algebra.save(fname2, outmap)

        # TODO: Noise needs to be Get rid of ignored frequencies.



    def correlate(self, lags=()) :
        """Calculate the cross correlation function of the maps.

        The cross correlation function is a funciton of f1, f2 and angular lag.
        
        The angular lag bins are passed, all pairs of frequencies are
        calculated.

        Parameters
        ----------
        lags : array like
            Angular lags bins (upper side bin edges).

        Returns
        -------
        corr : array
            Normalized such that the 0 lag, f1=f2 compenents are unity.
        """

        Map1 = self.Map1
        Map2 = self.Map2
        Noise1 = self.Noise_inv1
        Noise2 = self.Noise_inv2

        # Set up correlation function.
        if not lags :
            # Assume the pixel width in real degrees is the dec spacing.
            wid = abs(map1.info['dec_delta'])
            # Choose lags at fixed grid distance.
            lags = sp.array([0.0, 1.0, sp.sqrt(2.0), 2.0, sp.sqrt(5.0),
                             sp.sqrt(8.0), 3.0, sp.sqrt(10.0)])
            lags = (lags + 0.01)*wid
        freq1 = self.freq
        freq2 = self.freq

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
        ra_fact = (sp.cos(sp.pi*Map1.info['dec_centre'] / 180.0)
                   * Map1.info['ra_delta'])
        #dec_fact = Map1.info['dec_delta']

        # Calculate the pairwise lags.
        dra = (r1[:,None] - r2[None,:]) * ra_fact
        ddec = d1[:,None] - d2[None,:]
        lag = dra[:,None,:,None]**2 + ddec[None,:,None,:]**2
        lag = sp.sqrt(lag)
        # Bin this up.
        lag_inds = sp.digitize(lag.flatten(), lags)
        print "Starting Correlation."
        for if1 in range(len(freq1)) :
            for jf2 in range(len(freq2)) :
                # Calculate the pairwise products.
                data1 = map1[if1,:,:]
                data2 = map2[jf2,:,:]
                weights1 = noise1[if1,:,:]
                weights2 = noise2[jf2,:,:]
                dprod = data1[...,None,None] * data2[None,None,...]
                wprod = weights1[...,None,None] * weights2[None,None,...]
                for klag in range(nlags) :
                    mask = lag_inds == klag
                    corr[if1,jf2,klag] += sp.sum(dprod.flatten()[mask])
                    counts[if1,jf2,klag] += sp.sum(wprod.flatten()[mask])
        corr /= counts
        
        return corr


class NewSlices(object) :
    """Pipeline module that scripts together the cross correlation of maps.
    """
 
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
        # Get the map data from file as well as the noise inverse.  The Noise
        # data is huge matrix, but for now we are just using the diagonal.
        if len(params['file_middles']) == 1 :
            fn = params['file_middles'][0]
            params['file_middles'] = (fn, fn)
        if len(params['file_middles']) == 2 :
            map_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_map'])
            noise_file = (params['input_root'] + params['file_middles'][0] + 
                         params['input_end_noise'])
            Map1 = algebra.make_vect(algebra.load(map_file))
            print "Loading noise."
            
            # XXX For now make the noise a copy of the maps: loads faster.  Just
            # for testing.
            Noise_inv1 = algebra.make_vect(algebra.load(map_file))

            # Un comment these eventually.
            #Noise_inv1 = algebra.make_mat(algebra.open_memmap(noise_file,
            #                                                  mode='r'))
            #Noise_inv1 = Noise_inv1.mat_diag()
            #print "Done."
            
            map_file = (params['input_root'] + params['file_middles'][1] + 
                         params['input_end_map'])
            noise_file = (params['input_root'] + params['file_middles'][1] + 
                         params['input_end_noise'])
            
            Map2 = algebra.make_vect(algebra.load(map_file))
            print "Loading noise."
            # XXX For now make the noise a copy of the maps: loads faster.  Just
            # for testing.
            Noise_inv2 = algebra.make_vect(algebra.load(map_file))
            
            #Noise_inv2 = algebra.make_mat(algebra.open_memmap(noise_file,
            #                                                  mode='r'))
            #Noise_inv2 = Noise_inv2.mat_diag()
            #print "Done."
        else :
            raise ce.FileParameterTypeError('For now can only process one'
                                            ' or two files.')
        freq = sp.array(params['freq'], dtype=int)
        # Create pair of maps.
        Pair = MapPair(Map1, Map2, Noise_inv1, Noise_inv2, freq)
        # Hold a reference in self.
        self.Pair = Pair
        if params["convolve"] :
            Pair.degrade_resolution()
        if params['factorizable_noise'] :
            Pair.make_noise_factorizable()
        if params['sub_weighted_mean'] :
            Pair.subtract_weighted_mean()
        # Correlate the maps.
        self.fore_corr = Pair.correlate(params['lags'])
        self.fore_Pair = copy.deepcopy(Pair)

        # Subtract Foregrounds.
        # TODO: Provide a list of integers for params["modes"] so we can try a
        # few different numbers of modes to subtract.  This is computationally
        # expensive.
        vals, modes1, modes2 = get_freq_svd_modes(self.fore_corr, 
                                                  params['modes'])
        self.vals = vals
        self.modes1 = modes1
        self.modes2 = modes2
        # TODO: Add option to save the mode maps.
        Pair.subtract_frequency_modes(modes1, modes2)
        if params['first_pass_only'] :
            self.Pair = Pair
            return
        # Correlate the cleaned maps.
        # Here we could calculate the power spectrum instead eventually.
        corr = Pair.correlate(params['lags'])
        self.corr = corr

        if params["make_plots"] :
            self.make_plots()
            
    def make_plots(self) :
        plt.figure()
        plot_svd(self.vals)


def get_freq_svd_modes(corr, n) :
    """Same as get freq eigenmodes, but treats left and right maps
    separatly with an SVD.
    """
    U, s, V = linalg.svd(corr[:,:,0])
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
    
    return s, Lvectors, Rvectors

def plot_svd(vals) :
    """Plots the svd values and prints out some statistics."""
        
    n = len(vals)
    plt.semilogy(abs(sp.sort(-vals/n)), marker='o', 
                 linestyle='None')
    print 'Mean noise: ', sp.sum(vals)/n
    print 'Largest eigenvalues/n : ', 
    print sp.sort(vals/n)[-10:]

def rebin_core_freq_lag(corr, freq) :
    #TODO
    pass








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
            fn = params['file_middles'][0]
            params['file_middles'] = (fn, fn)
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
        print Map1.size
        print sp.sum(Map1 < 0), sp.sum(Map2 < 0)
        print sp.sum(Noise1 < 0), sp.sum(Noise2 < 0)
        print sp.sum(Noise1 == 0), sp.sum(Noise2 == 0)
        freq = sp.array(params['freq'], dtype=int)
        # Convolving down to a common resolution and reducing to factorized
        # noise.
        if params["convolve"] :
            # Get the beam data.
            gfreq=Map1.get_axis("freq")
            fwhm = utils.get_beam(gfreq)
            fwhm[:12]=fwhm[13]
            beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792, 
                     0.281176247549, 0.270856788455, 0.26745856078, 
                     0.258910010848, 0.249188429031])
            freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                                 dtype=float)
            freq_data *= 1.0e6
            beam_diff=sp.sqrt(max(1.1*beam_data)**2-(beam_data)**2)
            b = beam.GaussianBeam(beam_diff,freq_data)
            # Convolve to a common resolution.
            print "convolving"
            Map2=b.apply(Map2)
            Map1=b.apply(Map1)
            # Reduce noise to factorizable.
            Noise1[Noise1<1.e-30]=1.e-30
            Noise1=1/Noise1
            Noise1=b.apply(Noise1,cval=1.e30)
            Noise1_m=ma.array(Noise1)
            Noise1f=sp.mean(Noise1_m,0)
            Noise1_m[Noise1_m>1.e20]=ma.masked
            Noise1_m=Noise1_m/Noise1f[None,:,:]
            Noise1_m=(sp.mean(sp.mean(Noise1_m,1),1)[:,None,None] 
                      * Noise1f[None,:,:])
            Noise1_m=(1/Noise1_m).filled(0)
            Noise1[...] = Noise1_m

            
            Noise2[Noise2<1.e-30]=1.e-30
            Noise2=1/Noise2
            Noise2=b.apply(Noise2,cval=1.e30)
            Noise2_m=ma.array(Noise2)
            Noise2f=sp.mean(Noise2_m,0)
            Noise2_m[Noise2_m>1.e20]=ma.masked
            Noise2_m=Noise2_m/Noise2f[None,:,:]
            Noise2_m=(sp.mean(sp.mean(Noise2_m,1),1)[:,None,None]
                      * Noise2f[None,:,:])
            Noise2_m=(1/Noise2_m).filled(0)
            Noise2[...] = Noise2_m
        # Subtracting the frequency weighted mean from each slice.
        print sp.sum(Map1 < 0), sp.sum(Map2 < 0)
        print sp.sum(Noise1 < 0), sp.sum(Noise2 < 0)
        if params['sub_weighted_mean'] :
            means1 = sp.sum(sp.sum(Noise1*Map1, -1), -1)
            print means1
            means1 /= sp.sum(sp.sum(Noise1, -1), -1)
            print means1
            means1.shape += (1, 1)
            Map1 -= means1
            means2 = sp.sum(sp.sum(Noise2*Map2, -1), -1)
            print means2
            means2 /= sp.sum(sp.sum(Noise2, -1), -1)
            print means2
            means2.shape += (1, 1)
            Map2 -= means2
        print sp.sum(Map1 < 0), sp.sum(Map2 < 0)
        print sp.sum(Noise1 < 0), sp.sum(Noise2 < 0)
        # If any eigenmodes have been marked for subtraction, subtract them
        # out.
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
            outmap = sp.empty((len(Lmodes),)+Map1.shape[1:])
            outmap = algebra.make_vect(outmap,axis_names=('freq', 'ra', 'dec'))
            outmap.copy_axis_info(Map1)
            for ira in range(Map1.shape[1]) :
                for jdec in range(Map1.shape[2]) :
                    # if sp.any(Map1.data.mask[ira,jdec,freq]) :
                    #    continue
                    # else :
                    for i,v in enumerate(Lmodes) :
                        # v.shape = freq.shape
                        v = v.reshape(freq.shape)
                        # amp = sp.sum(v*Map1.data[ira,jdec,freq])
                        amp = sp.dot(v, Map1[freq,ira,jdec])
                        Map1[freq,ira,jdec] -= amp*v
                        outmap[i,ira,jdec] = sp.sum(amp)
            map_out = (params['output_root'] + params['file_middles'][0] + 
                                   '_cleaned' + params['input_end_map'])
            noise_out = (params['output_root'] + params['file_middles'][0] + 
                       '_cleaned' + params['input_end_noise'])
            map_out_modes = (params['output_root'] + params['file_middles'][0] + 
                                   '_Lmodes' + params['input_end_map'])
            # fits_map.write(Map1, map_out)
            algebra.save(map_out, Map1)
            algebra.save(noise_out, Noise1)
            algebra.save(map_out_modes, outmap)
            
            outmap = sp.empty((len(Rmodes),)+Map1.shape[1:])
            outmap = algebra.make_vect(outmap,axis_names=('freq', 'ra', 'dec'))
            outmap.copy_axis_info(Map2)
            for ira in range(Map2.shape[1]) :
                for jdec in range(Map2.shape[2]) :
                    # if sp.any(Map1.data.mask[ira,jdec,freq]) :
                    #    continue
                    # else :
                    for v in Lmodes :
                        # v.shape = freq.shape
                        v = v.reshape(freq.shape)
                        amp = sp.dot(v, Map2[freq,ira,jdec])
                        Map2[freq,ira,jdec] -= amp*v
                        outmap[i,ira,jdec] = sp.sum(amp)
            map_out = (params['output_root'] + params['file_middles'][1] + 
                       '_cleaned' + params['input_end_map'])
            noise_out = (params['output_root'] + params['file_middles'][1] + 
                       '_cleaned' + params['input_end_noise'])
            map_out_modes = (params['output_root'] + params['file_middles'][0] + 
                                   '_Rmodes' + params['input_end_map'])
            # fits_map.write(Map2, map_out)
            algebra.save(map_out, Map2)
            algebra.save(noise_out, Noise2)
            algebra.save(map_out_modes, outmap)

        print sp.sum(Map1 < 0), sp.sum(Map2 < 0)
        print sp.sum(Noise1 < 0), sp.sum(Noise2 < 0)
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
        ra_fact = (sp.cos(sp.pi*Map1.info['dec_centre'] / 180.0)
                   * Map1.info['ra_delta'])
        #dec_fact = Map1.info['dec_delta']

        # Calculate the pairwise lags.
        dra = (r1[:,None] - r2[None,:]) * ra_fact
        ddec = d1[:,None] - d2[None,:]
        lag = dra[:,None,:,None]**2 + ddec[None,:,None,:]**2
        lag = sp.sqrt(lag)
        # Bin this up.
        lag_inds = sp.digitize(lag.flatten(), lags)
        print "Starting Correlation."
        for if1 in range(len(freq1)) :
            for jf2 in range(len(freq2)) :
                # Calculate the pairwise products.
                data1 = map1[if1,:,:]
                data2 = map2[jf2,:,:]
                weights1 = noise1[if1,:,:]
                weights2 = noise2[jf2,:,:]
                dprod = data1[...,None,None] * data2[None,None,...]
                wprod = weights1[...,None,None] * weights2[None,None,...]
                for klag in range(nlags) :
                    mask = lag_inds == klag
                    corr[if1,jf2,klag] += sp.sum(dprod.flatten()[mask])
                    counts[if1,jf2,klag] += sp.sum(wprod.flatten()[mask])
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
        plt.semilogy(abs(sp.sort(-self.freq_svd_values/n)), marker='o', 
                     linestyle='None')
        print 'Mean noise: ', sp.sum(self.freq_svd_values)/n
        print 'Largest eigenvalues/n : ', 
        print sp.sort(self.freq_svd_values/n)[-10:]

def plot_contour(self, norms=False, lag_inds=(0)) :
    """Used as a method of FreqSlices.  A function instead for debugging."""
    
    lag_inds = list(lag_inds)
    # Set up binning.    
    nf = len(self.freq_inds)
    #freq_diffs = sp.sort(self.allfreq - min(self.allfreq))
    n_bins = 20
    factor = 1.5
    #start = 2.1e6
    freq_diffs = sp.empty(n_bins)
    freq_diffs[0] = 0.0001
    freq_diffs[1] = 2.5*200.0/256
    freq_diffs[2] = 4.5*200.0/256
    freq_diffs[3] = 6.5*200.0/256
    for ii in range(4,n_bins) :
       freq_diffs[ii] = factor*freq_diffs[ii-1]
    freq_diffs *= 1e6
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
    """Used as a method of FreqSlices.  A function instead for debugging."""
    
    lag_inds = list(lag_inds)
    # Set up binning.    
    nf = len(self.freq_inds)
    freq_diffs = sp.arange(0.1e6, 100e6, 200.0/256*1e6)
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
    # Cosmology dependant conersion to MPc.
    a_fact = 34.0
    f_fact = 4.5
    nbins = 10
    lags = sp.empty(nbins)
    lags[0] = 2.0
    lags[1] = 4.0
    for ii in range(1, nbins) :
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
    inds = pdat-errors[0,:]>0
    if sp.any(inds) :
        f = plt.errorbar(lags[inds], pdat[inds],
                         errors[:,inds], linestyle='None', marker='o', 
                         color='b', elinewidth=elin, markersize=msize)
    inds = pdat+errors[1,:]<0
    if sp.any(inds) :
        f = plt.errorbar(lags[inds], -pdat[inds], errors[:,inds], 
                         linestyle='None', marker='o', color='r', 
                         elinewidth=elin, markersize=msize)
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
        inds = pdat-errors[0,:]>0
        if sp.any(inds) :
            f = plt.errorbar(lags[inds], pdat[inds],
                             errors[:,inds], linestyle='None', 
                             marker='o', color='b', elinewidth=elin, mfc=mfc, 
                             markersize=msize)
        inds = pdat+errors[1,:]<0
        if sp.any(inds) :
            f = plt.errorbar(lags[inds], -pdat[inds],
                             errors[:,inds], linestyle='None',
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
    plt.axis([1.5, 100, 0.01, 500.0])
    plt.xlabel('lag (Mpc/h)')
    plt.ylabel('correlation (mK)')


# For running this module from the command line
if __name__ == '__main__' :
    import sys
    if len(sys.argv) == 2 :
        FreqSlices(str(sys.argv[1])).execute()
    elif len(sys.argv) > 2:
        print 'Maximum one argument, a parameter file name.'
    else :
        FreqSlices().execute()
        


