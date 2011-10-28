r"""Program that calculates the correlation function across frequency slices.
"""
# TODO: figure out data model for storing lags, etc. with class instance
# TODO: needs to continue to work with multiprocessing

import time
import scipy as sp
import numpy.ma as ma
from core import algebra
from map import beam
from utils import file_tools as ft

class MapPair(ft.ClassPersistence):
    r"""Pair of maps that are processed together and cross correlated.

    Parameters
    ----------
    map1, map2: algebra_vector
        Input Maps.
    noise_inv1, noise_inv2: algebra_vector
        Input Noise inverses.
    freq: tuple of ints
        The frequency indices to use.

    Attributes
    ----------
    map1, map2: algebra_vector
        Input Maps.
    noise_inv1, noise_inv2: algebra_vector
        Input Noise inverses.
    map1_name, map2_name: str
        The names of the maps from file_middles.
    map1_code, map2_code: str
        A single letter representing which section the map is.
    freq: tuple of ints
        The frequency indices to use.
    modes1, modes2: 2D array
        The modes for the 2 maps in a pair (from 1 to number of modes
        to be subtracted).
    left_modes, right_modes: 3D array
        What was subtracted from the maps during cleaning.

    """

    def __init__(self, *args, **kwargs):
        # variable names that define the object (for loading and saving)
        self.varlist = ['map1', 'map2', 'noise_inv1', 'noise_inv2', 'freq']
        self.varlist.extend(['counts', 'modes1', 'modes2'])
        self.varlist.extend(['left_modes', 'right_modes'])
        self.varlist.extend(['map1_name', 'map2_name'])

        if ((len(args) == 0) and ("shelve_filename" in kwargs)):
            self.shelve_init(*args, **kwargs)
        else:
            self.standard_init(*args, **kwargs)

    def shelve_init(self, *args, **kwargs):
        shelve_filename = kwargs['shelve_filename']
        self.load_variables(shelve_filename)

    def standard_init(self, *args, **kwargs):
        r"""
        arguments: map1, map2, noise_inv1, noise_inv2, freq
        """
        (self.map1, self.map2) = (args[0], args[1])
        (self.noise_inv1, self.noise_inv2) = (args[2], args[3])
        self.freq = args[4]

        # Give infinite noise to unconsidered frequencies (This doesn't affect
        # anything but the output maps).
        noise_dimensions = self.noise_inv1.shape[0]
        for noise_index in range(noise_dimensions):
            if not noise_index in self.freq:
                self.noise_inv1[noise_index, ...] = 0
                self.noise_inv2[noise_index, ...] = 0

        # Set attributes.
        self.counts = 0
        self.modes1 = 0
        self.modes2 = 0
        self.left_modes = 0
        self.right_modes = 0
        # For saving, to keep track of each mapname.
        self.map1_name = ''
        self.map2_name = ''
        # Which section [A, B, C, D...] the maps is from.
        self.map1_code = ''
        self.map2_code = ''

    def set_names(self, name1, name2):
        r"""Set the map names and codes.

        Set map1_name to name1 and map2_name to name2.
        Also the map codes. Note it is hardcoded for 4 maps right now.

        Parameters
        ----------
        name1, name2: str
            The names of the maps.
        """

        # Note that "I" would not be a good idea since it represents
        # a polarization and all of the maps have it.
        sections = ["A", "B", "C", "D", "E", "F", "G", "H"]

        self.map1_name = name1
        self.map2_name = name2

        # Set map1's section.
        found1 = False
        for letter in sections:
            if ((not found1) and (letter in name1)):
                found1 = True
                self.map1_code = letter

        # Set map2's section.
        found2 = False
        for letter in sections:
            if ((not found2) and (letter in name2)):
                found2 = True
                self.map2_code = letter

        if ((not found1) or (not found2)):
            print "Maps section can only be named A, B, C, D, E, F, G, or H."
            raise

    def degrade_resolution(self):
        r"""Convolves the maps down to the lowest resolution.

        Also convolves the noise, making sure to deweight pixels near the edge
        as well.  Converts noise to factorizable form by averaging.
        """
        noise1 = self.noise_inv1
        noise2 = self.noise_inv2

        # Get the beam data.
        beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031])
        freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                             dtype=float)
        freq_data *= 1.0e6
        beam_diff = sp.sqrt(max(1.1 * beam_data) ** 2 - (beam_data) ** 2)
        common_resolution = beam.GaussianBeam(beam_diff, freq_data)
        # Convolve to a common resolution.
        self.map2 = common_resolution.apply(self.map2)
        self.map1 = common_resolution.apply(self.map1)

        # This block of code needs to be split off into a function and applied
        # twice (so we are sure to do the same thing to each).
        noise1[noise1 < 1.e-30] = 1.e-30
        noise1 = 1. / noise1
        noise1 = common_resolution.apply(noise1, cval=1.e30)
        noise1 = 1. / noise1
        noise1[noise1 < 1.e-20] = 0

        noise2[noise2 < 1.e-30] = 1.e-30
        noise2 = 1 / noise2
        noise2 = common_resolution.apply(noise2, cval=1.e30)
        noise2 = 1. / noise2
        noise2[noise2 < 1.e-20] = 0

        self.noise_inv1 = algebra.as_alg_like(noise1, self.noise_inv1)
        self.noise_inv2 = algebra.as_alg_like(noise2, self.noise_inv2)

    def make_noise_factorizable(self):
        r"""Convert noise weights such that the factor into a function a
        frequency times a function of pixel by taking means over the original
        weights.
        """

        # TODO: move me elsewhere!
        def make_factorizable(noise):
            r"""factorize the noise"""
            noise[noise < 1.e-30] = 1.e-30
            noise = 1. / noise
            noise = ma.array(noise)
            # Get the freqency averaged noise per pixel.  Propagate mask in any
            # frequency to all frequencies.
            for noise_index in range(noise.shape[0]):
                if sp.all(noise[noise_index, ...] > 1.e20):
                    noise[noise_index, ...] = ma.masked
            noise_fmean = ma.mean(noise, 0)
            noise_fmean[noise_fmean > 1.e20] = ma.masked
            # Get the pixel averaged noise in each frequency.
            noise[noise > 1.e20] = ma.masked
            noise /= noise_fmean
            noise_pmean = ma.mean(ma.mean(noise, 1), 1)
            # Combine.
            noise = noise_pmean[:, None, None] * noise_fmean[None, :, :]
            noise = (1. / noise).filled(0)

            return noise

        noise_inv1 = make_factorizable(self.noise_inv1)
        noise_inv2 = make_factorizable(self.noise_inv2)
        self.noise_inv1 = algebra.as_alg_like(noise_inv1, self.noise_inv1)
        self.noise_inv2 = algebra.as_alg_like(noise_inv2, self.noise_inv2)

    def subtract_weighted_mean(self):
        r"""Subtracts the weighted mean from each frequency slice."""
        means1 = sp.sum(sp.sum(self.noise_inv1 * self.map1, -1), -1)
        means1 /= sp.sum(sp.sum(self.noise_inv1, -1), -1)
        means1.shape += (1, 1)
        self.map1 -= means1
        means2 = sp.sum(sp.sum(self.noise_inv2 * self.map2, -1), -1)
        means2 /= sp.sum(sp.sum(self.noise_inv2, -1), -1)
        means2.shape += (1, 1)
        self.map2 -= means2

        # Zero out all the infinit noise pixels (0 weight).
        self.map1[self.noise_inv1 < 1.e-20] = 0
        self.map2[self.noise_inv2 < 1.e-20] = 0

    def subtract_frequency_modes(self, modes1, modes2=None):
        r"""Subtract frequency mode from the map.

        This does not save anything anymore. The outmaps (L and R modes)
        that were saved before are now stored as a variable in the class
        and the saving of everything (maps, noise_invs, and modes) is done
        later in it's own function.

        Parameters
        ---------
        modes1: list of 1D arrays.
            Arrays must be the same length as self.freq.  Modes to subtract out
            of the map one.
        modes2: list of 1D arrays.
            Modes to subtract out of map 2.  If `None` set to `modes1`.

        """

        if modes2 == None:
            modes2 = modes1

        map1 = self.map1
        map2 = self.map2
        freq = self.freq

        # First map.
        outmap_left = sp.empty((len(modes1), ) + map1.shape[1:])
        outmap_left = algebra.make_vect(outmap_left,
                                     axis_names=('freq', 'ra', 'dec'))
        outmap_left.copy_axis_info(map1)
        for ira in range(map1.shape[1]):
            for jdec in range(map1.shape[2]):
                # if sp.any(map1.data.mask[ira, jdec, freq]):
                #    continue
                # else:
                for mode_index, mode_vector in enumerate(modes1):
                    # v.shape = freq.shape
                    mode_vector = mode_vector.reshape(freq.shape)
                    # amp = sp.sum(mode_vector*map1.data[ira, jdec, freq])
                    amp = sp.dot(mode_vector, map1[freq, ira, jdec])
                    map1[freq, ira, jdec] -= amp * mode_vector
                    outmap_left[mode_index, ira, jdec] = amp
        self.left_modes = outmap_left

        # Second map.
        outmap_right = sp.empty((len(modes2), ) + map2.shape[1:])
        outmap_right = algebra.make_vect(outmap_right,
                                     axis_names=('freq', 'ra', 'dec'))
        outmap_right.copy_axis_info(map2)
        for ira in range(map2.shape[1]):
            for jdec in range(map2.shape[2]):
                # if sp.any(map2.data.mask[ira, jdec, freq]):
                #    continue
                # else:
                for mode_index, mode_vector in enumerate(modes2):
                    # mode_vector.shape = freq.shape
                    mode_vector = mode_vector.reshape(freq.shape)
                    amp = sp.dot(mode_vector, map2[freq, ira, jdec])
                    map2[freq, ira, jdec] -= amp * mode_vector
                    outmap_right[mode_index, ira, jdec] = amp
        self.right_modes = outmap_right

    def correlate(self, lags=(), speedup=False, verbose=False):
        r"""Calculate the cross correlation function of the maps.

        The cross correlation function is a function of f1, f2 and angular lag.
        The angular lag bins are passed, all pairs of frequencies are
        calculated.

        Parameters
        ----------
        lags: array like
            Angular lags bins (upper side bin edges).
        speedup: boolean
            Speeds up the correlation. This works fine, yes? Should be the
            normal way if so.

        Returns
        -------
        corr: array
            The correlation between 2 maps.
        counts: array
            The weighting of the correlation based on the maps' weights.

        """

        map1 = self.map1
        map2 = self.map2
        noise1 = self.noise_inv1
        noise2 = self.noise_inv2

        freq1 = self.freq
        freq2 = self.freq

        map1_ra = map1.get_axis('ra')
        map2_ra = map2.get_axis('ra')
        map1_dec = map1.get_axis('dec')
        map2_dec = map2.get_axis('dec')

        #print lags
        #print map1_ra
        #print map2_ra
        #print map1_dec
        #print map2_dec
        #print map1.shape
        #print map2.shape
        #print noise1.shape
        #print noise2.shape
        #print freq1
        #print freq2

        input_map1 = map1[freq1, :, :]
        input_map2 = map2[freq2, :, :]
        input_noise1 = noise1[freq1, :, :]
        input_noise2 = noise2[freq2, :, :]

        # Noise weight
        input_map1 *= input_noise1
        input_map2 *= input_noise2

        nlags = len(lags)
        nfreq = len(freq1)
        corr = sp.zeros((nfreq, nfreq, nlags), dtype=float)
        counts = sp.zeros(corr.shape, dtype=float)
        # Noting that if DEC != 0, then a degree of RA is less than a degree.
        ra_fact = sp.cos(sp.pi * map1.info['dec_centre'] / 180.0)

        # Calculate the pairwise lags.
        dra = (map1_ra[:, None] - map2_ra[None, :]) * ra_fact
        ddec = map1_dec[:, None] - map2_dec[None, :]
        lag = dra[:, None, :, None] ** 2 + ddec[None, :, None, :] ** 2
        lag = sp.sqrt(lag)
        # Bin this up.
        lag_inds = sp.digitize(lag.flatten(), lags)

        if speedup:
            print "Starting Correlation (sparse version)"

            (nr1, nd1) = (len(map1_ra), len(map1_dec))
            (nr2, nd2) = (len(map2_ra), len(map2_dec))
            (r1ind, d1ind) = (sp.arange(nr1), sp.arange(nd1))
            (r2ind, d2ind) = (sp.arange(nr2), sp.arange(nd2))
            ra1_pairind = r1ind.repeat(nr2 * nd1 * nd2)
            ra2_pairind = sp.tile(r2ind.repeat(nd2), (1, nr1 * nd1)).flatten()
            dec1_pairind = sp.tile(d1ind.repeat(nr2 * nd2), (1, nr1)).flatten()
            dec2_pairind = sp.tile(d2ind, (1, nr1 * nr2 * nd1)).flatten()

            # precalculate the pair indices for a given lag
            # could also imagine calculating the map slices here
            posmaskdict = {}
            for klag in range(nlags):
                mask = (lag_inds == klag)
                posmaskdict[repr(klag)] = (ra1_pairind[mask], ra2_pairind[mask],
                                           dec1_pairind[mask], dec2_pairind[mask])

            for if1 in range(len(freq1)):
                for jf2 in range(len(freq2)):
                    start = time.time()

                    data1 = input_map1[if1, :, :]
                    data2 = input_map2[jf2, :, :]
                    weights1 = input_noise1[if1, :, :]
                    weights2 = input_noise2[jf2, :, :]

                    for klag in range(nlags):
                        # cached, but written out:
                        #dprod = data1[ra1_pairind[mask], dec1_pairind[mask]]* \
                        #        data2[ra2_pairind[mask], dec2_pairind[mask]]
                        #wprod = weights1[ra1_pairind[mask], dec1_pairind[mask]] * \
                        #        weights2[ra2_pairind[mask], dec2_pairind[mask]]
                        (r1m, r2m, d1m, d2m) = posmaskdict[repr(klag)]
                        dprod = data1[r1m, d1m] * data2[r2m, d2m]
                        wprod = weights1[r1m, d1m] * weights2[r2m, d2m]
                        corr[if1, jf2, klag] += sp.sum(dprod)
                        counts[if1, jf2, klag] += sp.sum(wprod)

                    if verbose:
                        print if1, jf2, (time.time() - start), counts[if1, jf2, :]
        else:
            print "Starting Correlation (full version)"
            for if1 in range(len(freq1)):
                for jf2 in range(len(freq2)):
                    start = time.time()
                    # Calculate the pairwise products.
                    data1 = input_map1[if1, :, :]
                    data2 = input_map2[jf2, :, :]
                    weights1 = input_noise1[if1, :, :]
                    weights2 = input_noise2[jf2, :, :]
                    dprod = data1[..., None, None] * data2[None, None, ...]
                    wprod = weights1[..., None, None] * \
                            weights2[None, None, ...]
                    for klag in range(nlags):
                        mask = (lag_inds == klag)
                        corr[if1, jf2, klag] += sp.sum(dprod.flatten()[mask])
                        counts[if1, jf2, klag] += sp.sum(wprod.flatten()[mask])

                    if verbose:
                        print if1, jf2, (time.time() - start), counts[if1, jf2, :]

        mask = (counts < 1e-20)
        counts[mask] = 1
        corr /= counts
        corr[mask] = 0
        counts[mask] = 0

        return corr, counts
