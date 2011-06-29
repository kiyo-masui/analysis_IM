
"""Program that calculates the correlation function across frequency slices.
"""

import copy
import multiprocessing
import time
import os.path as path
import cPickle
import random as rn

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
               'factorizable_noise' : False,
               'no_weights' : False
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
        
        # Give infinite noise to unconsidered frequencies (This doesn't affect
        # anything but the output maps).
        n = Noise_inv1.shape[0]
        for ii in range(n) :
            if not ii in freq :
                Noise_inv1[ii,...] = 0
                Noise_inv2[ii,...] = 0

        # Set attributes.
        self.Map1 = Map1
        self.Map2 = Map2
        self.Noise_inv1 = Noise_inv1
        self.Noise_inv2 = Noise_inv2
        self.freq = freq
        # To be set by methods:
        # NOTE fore_Pairs will not have these set, only Pairs.
        self.fore_corr = 0
        self.corr = 0
        self.fore_counts = 0
        self.counts = 0
        self.vals = 0
        self.modes1 = 0
        self.modes2 = 0
        # For multiprocessing function since this can't be passed in.
        self.lags = 0
        # For saving, to keep track of each mapname.
        self.Map1_name = ''
        self.Map2_name = ''
        # Which section [A,B,C,D...] the maps is from.
        self.Map1_code = ''
        self.Map2_code = ''
        # To have information about where things are.
        self.params = {}

    def set_names(self, name1, name2):
        """Set Map1_name to name1 and Map2_name to name2. 
        Also the map codes. Note it is hardcoded for 4 maps right now."""
        self.Map1_name = name1
        self.Map2_name = name2
        # Set Map1's section.
        if "A" in name1:
            self.Map1_code = "A"
        elif "B" in name1:
            self.Map1_code = "B"
        elif "C" in name1:
            self.Map1_code = "C"
        elif "D" in name1:
            self.Map1_code = "D"
        else:
            print "Maps can only be names by sections A,B,C or D."
            raise
        # Set Map2's section.
        if "A" in name2:
            self.Map2_code = "A"
        elif "B" in name2:
            self.Map2_code = "B"
        elif "C" in name2:
            self.Map2_code = "C"
        elif "D" in name2:
            self.Map2_code = "D"
        else:
            print "Maps can only be names by sections A,B,C or D."
            raise

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
        Noise1[Noise1<1.e-30] = 1.e-30
        Noise1 = 1./Noise1
        Noise1 = b.apply(Noise1, cval=1.e30)
        Noise1 = 1./Noise1
        Noise1[Noise1<1.e-20] = 0

        Noise2[Noise2<1.e-30] = 1.e-30
        Noise2 = 1/Noise2
        Noise2 = b.apply(Noise2, cval=1.e30)
        Noise2 = 1./Noise2
        Noise2[Noise2<1.e-20] = 0
        
        self.Noise_inv1 = algebra.as_alg_like(Noise1, self.Noise_inv1)
        self.Noise_inv2 = algebra.as_alg_like(Noise2, self.Noise_inv2)

    def make_noise_factorizable(self) :
        """Convert noise weights such that the factor into a function a
        frequecy times a function of pixel by taking means over the origional
        weights.
        """
        
        def make_factorizable(Noise) :
            # Take the reciprical.
            Noise[Noise<1.e-30] = 1.e-30
            Noise = 1./Noise
            Noise = ma.array(Noise)
            # Get the freqency averaged noise per pixel.  Propagate mask in any
            # frequency to all frequencies.
            for ii in range(Noise.shape[0]) :
                if sp.all(Noise[ii,...]>1.e20):
                    Noise[ii,...] = ma.masked
            Noise_fmean = ma.mean(Noise, 0)
            Noise_fmean[Noise_fmean>1.e20] = ma.masked
            # Get the pixel averaged noise in each frequency.
            Noise[Noise>1.e20] = ma.masked
            Noise /= Noise_fmean
            Noise_pmean = ma.mean(ma.mean(Noise, 1), 1)
            # Combine.
            Noise = Noise_pmean[:,None,None] * Noise_fmean[None,:,:]
            Noise = (1./Noise).filled(0)

            return Noise

        Noise_inv1 = make_factorizable(self.Noise_inv1)
        Noise_inv2 = make_factorizable(self.Noise_inv2)
        self.Noise_inv1 = algebra.as_alg_like(Noise_inv1, self.Noise_inv1)
        self.Noise_inv2 = algebra.as_alg_like(Noise_inv2, self.Noise_inv2)

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

        # Zero out all the infinit noise pixels (0 weight).
        self.Map1[self.Noise_inv1<1.e-20] = 0
        self.Map2[self.Noise_inv2<1.e-20] = 0

    def subtract_frequency_modes(self, modes1, modes2=None, fname1=None, 
                                 fname2=None) :
        """Subtract frequency mode from the map. fnames will have the
        output root added to the beginning. !!!DOES NOT SAVE MODES NOW, 
        SAVES THE MAPS AFTER SUBRACTION.

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

        if not (fname1 is None):
            other_mapname1 = self.Map2_code
        if not (fname2 is None):
            other_mapname2 = self.Map1_code
        
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
            save_file = self.params['output_root'] + fname1 + \
                            "_cleaned_clean_map_I_with_" + \
                            other_mapname1 + ".npy"
            #algebra.save(fname1, outmap)
            algebra.save(save_file, Map1) # outmap)
        
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
            save_file = self.params['output_root'] + fname2 + \
                            "_cleaned_clean_map_I_with_" + \
                            other_mapname2 + ".npy"
            #algebra.save(fname1, outmap)
            algebra.save(save_file, Map2) # outmap)


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
                print str(if1) + " : " + str(jf2)
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

        # Should probably return counts as well, since this can be used as
        # noise weight.
        return corr, counts


class NewSlices(object) :
    """Pipeline module that scripts together the cross correlation of maps.
    Now: params, Pair, fore_corr, corr, fore_Pair, vals, modes1, modes2.
    Want: Pairs, corr_final, core_std.
    """
 
    def __init__(self, parameter_file_or_dict=None) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix)

    def execute(self) :
        params = self.params
        freq = sp.array(params['freq'], dtype=int)
        lags = sp.array(params['lags'])
        # Write parameter file.
        kiyopy.utils.mkparents(params['output_root'])
        parse_ini.write_params(params, params['output_root'] + 'params.ini',
                               prefix=prefix)
        # Get the map data from file as well as the noise inverse.  The Noise
        # data is huge matrix, but for now we are just using the diagonal.
        if len(params['file_middles']) == 1 :
            fn = params['file_middles'][0]
            params['file_middles'] = (fn, fn)
        if len(params['file_middles']) >= 2 :
#            raise ce.FileParameterTypeError('For now can only process one'
#                                            ' or two files.')
            # Deal with multiple files.
            num_maps = len(params['file_middles'])
            Maps = []
            Noise_invs = []
            # Load all maps and noises once.
            for ii in range(0, num_maps):
                map_file = (params['input_root'] + 
                    params['file_middles'][ii] + params['input_end_map'])
                print "Loading map %d of %d." %(ii+1, num_maps)
                Map = algebra.make_vect(algebra.load(map_file))
                Maps.append(Map)
                if not params["no_weights"] :
                    noise_file = (params['input_root'] + 
                        params['file_middles'][ii] + params['input_end_noise'])
                    print "Loading noise %d of %d." %(ii+1, num_maps)
    #                Noise_invs.append(abs(algebra.make_vect(algebra.load(map_file))))
                    Noise_inv = algebra.make_mat(algebra.open_memmap(noise_file,
                                                                      mode='r'))
                    Noise_inv = Noise_inv.mat_diag()
                else :
                    Noise_inv = algebra.ones_like(Map)
                    #Noise_inv[...] += (100 + 
                    #    ii*sp.arange(Noise_inv.shape[0])[:,None,None])
                    #Noise_inv[...] += (100 +
                    #    ii*sp.arange(Noise_inv.shape[1])[None,:,None])

                Noise_invs.append(Noise_inv)
            Pairs = []
            # Make pairs with deepcopies to not make mutability mistakes. 
            for ii in range(0, num_maps):
                for jj in range(0, num_maps):
                    if (jj > ii):
                        Map1 = copy.deepcopy(Maps[ii])
                        Map2 = copy.deepcopy(Maps[jj])
                        Noise_inv1 = copy.deepcopy(Noise_invs[ii])
                        Noise_inv2 = copy.deepcopy(Noise_invs[jj])
                        Pair = MapPair(Map1, Map2,
                            Noise_inv1, Noise_inv2, freq)
                        Pair.lags = lags
                        Pair.params = params
                        Pair.set_names(params['file_middles'][ii],
                                       params['file_middles'][jj])
                        Pairs.append(Pair)
            num_map_pairs = len(Pairs)
            print "%d Map Pairs created from %d maps." %(len(Pairs), num_maps)
        # Hold a reference in self.
        self.Pairs = Pairs


        if params["convolve"] :
            for Pair in Pairs:
                Pair.degrade_resolution()
        if params['factorizable_noise'] :
            for Pair in Pairs:
                Pair.make_noise_factorizable()
        if params['sub_weighted_mean'] :
            for Pair in Pairs:
                Pair.subtract_weighted_mean()

        # Correlate the maps.
        # Should not return anything, fore_corr should be saved in MapPair.
# Before multiprocessing.
#        self.fore_Pairs = []
#        for Pair in Pairs:
#            fore_corr, fore_counts = Pair.correlate(lags)
#            Pair.fore_corr = fore_corr
#            Pair.fore_counts = fore_counts
#            self.fore_Pairs.append(copy.deepcopy(Pair))
#        return
# After :
        fore_Pairs = []
        processes_list = []
#        for Pair in Pairs:
        for ii in range(0, num_map_pairs):
#            fore_muliproc(Pair)
            p = multiprocessing.Process(target=multiproc,
                                        args=([Pairs[ii], ii, False]))
            processes_list.append(p)
            p.start()
            #print(parent_conn.recv())

#        while p.is_alive():
        while True in [p.is_alive() for p in processes_list]:
            print "processing"
            time.sleep(5)
        # just to be safe
        time.sleep(1)
        print "Loading Map Pairs back into program [a min or so wait]"
        file_name = "Map_Pair_for_freq_slices_fore_corr_"
        for count in range(0, num_map_pairs):
#        while path.exists(file_name+str(count)+".pkl"):
            print "Loading correlation for Pair %d" %(count)
            f = open(file_name+str(count)+".pkl", "r")
            Pairs[count].fore_corr = cPickle.load(f)
            fore_Pairs.append(Pairs[count])
            f.close()
        self.fore_Pairs = copy.deepcopy(fore_Pairs)
        # To have fore_corr in self.Pairs, too, to not need
        # self.fore_Pairs later.
        self.Pairs = copy.deepcopy(fore_Pairs)
        print "gung ho!"


                
        Pairs = self.Pairs
        # Subtract Foregrounds.
        # TODO: Provide a list of integers for params["modes"] so we can try a
        # few different numbers of modes to subtract.  This is computationally
        # expensive.
        for Pair in Pairs:
            # Since these values are different for diff maps,
            # can there be a 'final' one? [Like an average?]
            vals, modes1, modes2 = get_freq_svd_modes(Pair.fore_corr, 
                                                  params['modes'])
            Pair.vals = vals
            Pair.modes1 = modes1
            Pair.modes2 = modes2
 #       self.vals = vals
 #       self.modes1 = modes1
 #       self.modes2 = modes2

        # TODO: Add option to save the mode maps.
#        for Pair in Pairs:
        for ii in range(0, len(Pairs)):
            
            Pairs[ii].subtract_frequency_modes(Pairs[ii].modes1, 
                Pairs[ii].modes2, Pairs[ii].Map1_name, Pairs[ii].Map2_name)
        if params['first_pass_only'] :
            self.Pairs = Pairs
            return
        # Correlate the cleaned maps.
        # Here we could calculate the power spectrum instead eventually.

# Before multiprocessing:
#        for Pair in Pairs:
#            corr, counts = Pair.correlate(lags)
#            Pair.corr = corr
#            Pair.counts = counts
# After:
        # Same as previous multiproc, but this time, you get the final
        # correlation for each Pair.
        temp_Pair_list = []
        processes_list = []
        for ii in range(0, num_map_pairs):
            p = multiprocessing.Process(target=multiproc,
                                        args=([Pairs[ii], ii, True]))
            processes_list.append(p)
            p.start()

        while True in [p.is_alive() for p in processes_list]:
            print "processing"
            time.sleep(5)
        # just to be safe
        time.sleep(1)
        print "Loading Map Pairs back into program [a min or so wait]"
        file_name = "Map_Pair_for_freq_slices_corr_"
        for count in range(0, num_map_pairs):
            print "Loading correlation for Pair %d" %(count)
            f = open(file_name+str(count)+".pkl", "r")
            Pairs[count].corr = cPickle.load(f)
            temp_Pair_list.append(Pairs[count])
            f.close()
        self.Pairs = copy.deepcopy(temp_Pair_list)
        print "gung ho!"
        

        # Get the average correlation and its standard deviation.
        corr_list = []
        for Pair in self.Pairs:
            corr_list.append(Pair.corr)
        self.corr_final, self.corr_std = get_corr_and_std_3D(corr_list)

        if params["make_plots"] :
            print "Plots not supported for multiple pairs."
 #           self.make_plots()
        return
            
    def make_plots(self) :
        plt.figure()
        plot_svd(self.vals)
        plt.show()

def multiproc(Pair, pair_number, final):
    """Do the correlation for a MapPair Pair. Pair is then saved
    to a numbered file corresponding to pair_number. This is done since
    this function gets multiprocessed and must have a common information
    ground somewhere. If final is false, the fore correlation is being done.
    If final is True, the final correlation is being done. final is used
    to save Pair to the appropriate file."""
    print "I am starting."
#    print Pair.lags
#    print type(Pair.lags)
#    fore_core, fore_counts = Pair.correlate(Pair.lags)
#    Pair.shout(44)
    control_correlation(Pair, Pair.lags, final)
    if final:
        file_name = "Map_Pair_for_freq_slices_corr_" + \
                        str(pair_number) + ".pkl"
        to_save = Pair.corr
    else:
        file_name = "Map_Pair_for_freq_slices_fore_corr_" + \
                        str(pair_number) + ".pkl"
        to_save = Pair.fore_corr
    f = open(file_name, "w")
    print "Writing to: ",
    print file_name
    cPickle.dump(to_save,f)
    f.close()
    print "I have have gotten here. "
    return

def control_correlation(Pair, lags, final):
    """For multiprocessing since it cannot deal with return values.
    Save the return values of correlate(lags) in self.corr and .counts
    if final is True. Else, save return values in self.fore_corr and
    .fore_counts."""
     
    if final:
        Pair.corr, Pair.counts = Pair.correlate(lags)
    else:
        Pair.fore_corr, Pair.fore_counts = Pair.correlate(lags)

def get_corr_and_std_3D(corr_list):
    '''Return the average correlation and the std for each point in
    the 3D correlation. corr_list is a list of 3D correlation matrices
    of length > 0.'''
    # Get average.
    corr_sum = corr_list[0]
    for ii in range(1, len(corr_list)):
        corr_sum += corr_list[ii]
    corr_avg = corr_sum/float(len(corr_list))
    # Get std. Note it will be all 0 if only one pair was used.
    xx,yy,zz = corr_list[0].shape
    corr_std = sp.zeros((xx,yy,zz))
    for ii in range(0,xx):
        for jj in range(0,yy):
            for kk in range(0,zz):
                value_list = []
                for corr in corr_list:
                    value_list.append(corr[ii,jj,kk])
                corr_std[ii,jj,kk] = sp.std(value_list)
    return corr_avg, corr_std

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

def subtract_modes_corr(corr, n) :
    """Similar to get_freq_svd_modes, but the modes are subtracted and from the
    correlation and the cleaned correlation is returned."""

    corr = copy.deepcopy(corr)

    s, Lvectors, Rvectors = get_freq_svd_modes(corr, n)
    for ii in range(n) :
        corr -= s[ii] * Lvectors[ii][:,None,None] * Rvectors[ii][None,:,None]

    return corr

def plot_svd(vals) :
    """Plots the svd values and prints out some statistics."""
        
    n = len(vals)
    plt.semilogy(abs(sp.sort(-vals/n)), marker='o', 
                 linestyle='None')
    print 'Mean noise: ', sp.sum(vals)/n
    print 'Largest eigenvalues/n : ', 
    print sp.sort(vals/n)[-10:]

def normalize_corr(corr):
    """Return the normalized 3D correlation along the diagonal. 
    [3D meaning f,f_prime,lag]."""
    # Get dimensions.
    fs,f_primes,lags = corr.shape
    corr_norm = sp.zeros((fs,f_primes,lags))
    # Each angular lag is separate.
    for lag in range(0, lags):
        # At each pixel, get the geometric mean of the 2 freqs along the
        # diagonal [ie. at C[f,f,lag] and C[f_prime,f_prime,lag]). Divide
        # the pixel by that value to make it normalized.
        for f in range(0, fs):
            print lag,
            print f
            for f_prime in range(0, f_primes):
                value = corr[f,f,lag] * corr[f_prime,f_prime,lag]
                factor = sp.sqrt(value)
                corr_norm[f,f_prime,lag] = corr[f,f_prime,lag] / factor
    return corr_norm
        
    

def rebin_corr_freq_lag(corr, freq1, freq2=None, weights=None, nfbins=20) :
    """Collapses frequency pair correlation function to frequency lag.
    
    Basically this constructs the 2D correlation function.

    This function takes in a 3D corvariance matrix which is a function of
    frequency and frequency prime and returns the same matrix but as a function
    of frequency - frequency prime (2D).

    freq1 is the actualy freq values not the indeces [700mh not 0,1,2...]
    for weights see counts in correlate.
    """
    
    if freq2 is None :
        freq2 = freq1
    # Default is equal weights.
    if weights is None :
        weights = sp.zeros_like(corr) + 1.0
    corr *= weights
    
    nf1 = corr.shape[0]
    nf2 = corr.shape[1]
    nlags = corr.shape[2]
    # Frequency bin size.
    df = min(abs(sp.diff(freq1)))
    # Frequency bin upper edges.
    fbins = sp.arange(1, nfbins+1)*df
    # Allowcate memory for outputs.
    out_corr = sp.zeros((nfbins, nlags))
    out_weights = sp.zeros((nfbins, nlags))
    
    # Loop over all frequency pairs and bin by lag.
    for ii in range(nf1) :
        print ii
        for jj in range(nf2) :
            f_lag = abs(freq1[ii] - freq2[jj])
            bin_ind = sp.digitize([f_lag], fbins)[0]
            if bin_ind < nfbins:
                out_corr[bin_ind,:] += corr[ii,jj,:]
                out_weights[bin_ind,:] += weights[ii,jj,:]
    # Normalize dealing with 0 weight points explicitly.
    bad_inds = out_weights < 1.0e-20
    out_weights[bad_inds] = 1.0
    out_corr/=out_weights
    out_weights[bad_inds] = 0.0

    return out_corr, out_weights








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
        


