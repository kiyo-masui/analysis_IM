r"""Program that calculates the correlation function across frequency slices.
"""
import gc
import scipy as sp
import numpy as np
import numpy.ma as ma
from core import algebra
from map import beam
from map import physical_gridding as pg
import kiyopy.utils
from quadratic_products import pwrspec_estimator as pe
from foreground_clean import find_modes
from utils import batch_handler as bh
# TODO: move single map operations to a separate class


class MapPair(object):
    r"""Pair of maps that are processed together and cross correlated.
    """

    def __init__(self, map1, map2, noise_inv1, noise_inv2, freq,
                 input_filenames=False):
        r"""
        arguments: map1, map2, noise_inv1, noise_inv2, freq
        """
        if input_filenames:
            self.map1 = algebra.make_vect(algebra.load(map1))
            self.map2 = algebra.make_vect(algebra.load(map2))
            if noise_inv1:
                self.noise_inv1 = algebra.make_vect(algebra.load(noise_inv1))
            else:
                print "WARNING: map1 has unity weight; no file given"
                self.noise_inv1 = algebra.ones_like(self.map1)

            if noise_inv2:
                self.noise_inv2 = algebra.make_vect(algebra.load(noise_inv2))
            else:
                print "WARNING: map2 has unity weight; no file given"
                self.noise_inv2 = algebra.ones_like(self.map2)

        else:
            self.map1 = map1
            self.map2 = map2
            self.noise_inv1 = noise_inv1
            self.noise_inv2 = noise_inv2

        self.freq = freq

        # maps in physical coordinates (derived)
        self.phys_map1 = None
        self.phys_map2 = None
        self.phys_noise_inv1 = None
        self.phys_noise_inv2 = None

        # give infinite noise to masked bands
        self.sanitize()

        # Set attributes.
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

    def sanitize(self):
        r"""set weights to zero in funny regions"""
        print "sanitizing the input arrays and weights"
        (self.map1, self.noise_inv1) = self.sanitize_single(self.map1,
                                                       self.noise_inv1)

        (self.map2, self.noise_inv2) = self.sanitize_single(self.map2,
                                                       self.noise_inv2)

    def sanitize_single(self, map, weightmap):
        weightmap[np.isnan(weightmap)] = 0.
        weightmap[np.isinf(weightmap)] = 0.
        weightmap[np.isnan(map)] = 0.
        weightmap[np.isinf(map)] = 0.
        weightmap[weightmap < 1.e-20] = 0.

        # This is not needed for sane maps, but some sims that to mean
        # subtraction before this stage could have nans by accident.
        # (np.nan * 0. is still np.nan)
        # not setting map[weightmap < 1.e-20] = 0 because these could be
        # "real" ares of the map which are mixed in through the common res.
        # convolution.
        map[np.isnan(map)] = 0.
        map[np.isinf(map)] = 0.

        weight_dimensions = weightmap.shape[0]
        for freq_index in range(weight_dimensions):
            if not freq_index in self.freq:
                weightmap[freq_index, ...] = 0

        #map[weightmap < 1.e-20] = 0.

        return (map, weightmap)

    def make_physical(self, refinement=1, pad=5, order=2, clear_obs=True):
        r"""Project the maps and weights into physical coordinates
        refinement allows the physical bins to be e.g. =2 times finer
        pad puts a number of padded pixels on all sides of the physical vol.
        the clear_obs flag deletes the obs. coord maps
        """

        self.phys_map1 = bh.repackage_kiyo(pg.physical_grid(self.map1,
                                           refinement=refinement,
                                           pad=pad, order=order))

        self.phys_map2 = bh.repackage_kiyo(pg.physical_grid(self.map2,
                                           refinement=refinement,
                                           pad=pad, order=order))

        self.phys_noise_inv1 = bh.repackage_kiyo(pg.physical_grid(
                                                 self.noise_inv1,
                                                 refinement=refinement,
                                                 pad=pad, order=order))

        self.phys_noise_inv2 = bh.repackage_kiyo(pg.physical_grid(
                                                 self.noise_inv2,
                                                 refinement=refinement,
                                                 pad=pad, order=order))

        del self.map1, self.map2, self.noise_inv1, self.noise_inv2
        gc.collect()

        return

    def degrade_resolution(self):
        r"""Convolves the maps down to the lowest resolution.

        Also convolves the noise, making sure to deweight pixels near the edge
        as well.  Converts noise to factorizable form by averaging.
        """
        print "degrading the resolution to a common beam"
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
        noise1[noise1 < 1.e-20] = 0.

        noise2[noise2 < 1.e-30] = 1.e-30
        noise2 = 1 / noise2
        noise2 = common_resolution.apply(noise2, cval=1.e30)
        noise2 = 1. / noise2
        noise2[noise2 < 1.e-20] = 0.

        self.noise_inv1 = algebra.as_alg_like(noise1, self.noise_inv1)
        self.noise_inv2 = algebra.as_alg_like(noise2, self.noise_inv2)

    def make_noise_factorizable(self, weight_prior=2):
        r"""Convert noise weights such that the factor into a function a
        frequency times a function of pixel by taking means over the original
        weights.

        weight_prior used to be 10^-30 before prior applied
        """
        print "making the noise factorizable"

        def make_factorizable(noise):
            r"""factorize the noise"""
            noise[noise < weight_prior] = 1.e-30
            noise = 1. / noise
            noise = ma.array(noise)
            # Get the freqency averaged noise per pixel.  Propagate mask in any
            # frequency to all frequencies.
            for noise_index in range(ma.shape(noise)[0]):
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
        print "subtracting the weighted mean from each slice"
        means1 = sp.sum(sp.sum(self.noise_inv1 * self.map1, -1), -1)
        means1 /= sp.sum(sp.sum(self.noise_inv1, -1), -1)
        means1.shape += (1, 1)
        self.map1 -= means1
        means2 = sp.sum(sp.sum(self.noise_inv2 * self.map2, -1), -1)
        means2 /= sp.sum(sp.sum(self.noise_inv2, -1), -1)
        means2.shape += (1, 1)
        self.map2 -= means2

        # Zero out all the infinite noise pixels (0 weight).
        self.map1[self.noise_inv1 < 1.e-20] = 0.
        self.map2[self.noise_inv2 < 1.e-20] = 0.

    def apply_map_weights(self):
        self.map1 = self.map1 * self.noise_inv1
        self.map2 = self.map2 * self.noise_inv2

    def unapply_map_weights(self):
        self.map1 = self.map1 / self.noise_inv1
        self.map2 = self.map2 / self.noise_inv2

    def subtract_frequency_modes(self, modes1, modes2=None,
                                 weighted=False, defer=False):
        r"""Subtract frequency modes from the map.
        """

        if modes2 == None:
            modes2 = modes1

        # First map.
        outmap_left = sp.empty((len(modes1), ) + self.map1.shape[1:])
        outmap_left = algebra.make_vect(outmap_left,
                                     axis_names=('freq', 'ra', 'dec'))
        outmap_left.copy_axis_info(self.map1)

        if defer:
            fitted = np.zeros_like(self.map1[self.freq, :, :])

        for mode_index, mode_vector in enumerate(modes1):
            mode_vector = mode_vector.reshape(self.freq.shape)

            if weighted:
                amp = sp.tensordot(mode_vector, self.map1[self.freq, :, :] *
                                self.noise_inv1[self.freq, :, :], axes=(0,0))
                amp /= sp.tensordot(mode_vector, mode_vector[:, None, None] *
                                self.noise_inv1[self.freq, :, :], axes=(0,0))
            else:
                amp = sp.tensordot(mode_vector,
                                   self.map1[self.freq, :, :], axes=(0,0))
                #amp /= sp.dot(mode_vector, mode_vector)

            if defer:
                fitted += mode_vector[:, None, None] * amp[None, :, :]
            else:
                fitted = mode_vector[:, None, None] * amp[None, :, :]
                self.map1[self.freq, :, :] -= fitted

            outmap_left[mode_index, :, :] = amp

        if defer:
            self.map1 -= fitted

        self.left_modes = outmap_left

        # Second map.
        outmap_right = sp.empty((len(modes2), ) + self.map2.shape[1:])
        outmap_right = algebra.make_vect(outmap_right,
                                     axis_names=('freq', 'ra', 'dec'))
        outmap_right.copy_axis_info(self.map2)

        if defer:
            fitted = np.zeros_like(self.map2[self.freq, :, :])

        for mode_index, mode_vector in enumerate(modes2):
            mode_vector = mode_vector.reshape(self.freq.shape)

            if weighted:
                amp = sp.tensordot(mode_vector, self.map2[self.freq, :, :] *
                                self.noise_inv2[self.freq, :, :], axes=(0,0))
                amp /= sp.tensordot(mode_vector, mode_vector[:, None, None] *
                                self.noise_inv2[self.freq, :, :], axes=(0,0))
            else:
                amp = sp.tensordot(mode_vector,
                                   self.map2[self.freq, :, :], axes=(0,0))
                #amp /= sp.dot(mode_vector, mode_vector)

            if defer:
                fitted += mode_vector[:, None, None] * amp[None, :, :]
            else:
                fitted = mode_vector[:, None, None] * amp[None, :, :]
                self.map2[self.freq, :, :] -= fitted

            outmap_right[mode_index, :, :] = amp

        if defer:
            self.map2 -= fitted

        self.right_modes = outmap_right

    def pwrspec_summary(self, window=None, unitless=True, bins=None,
                    truncate=False, nbins=40, logbins=True,
                    refinement=2, pad=5, order=2, return_3d=False):
        r"""calculate the 1D power spectrum
        """
        self.make_physical(refinement=refinement, pad=pad, order=order)

        xspec = pe.calculate_xspec(self.phys_map1, self.phys_map2,
                                   self.phys_noise_inv1, self.phys_noise_inv2,
                                   window=window, unitless=unitless,
                                   bins=bins, truncate=truncate, nbins=nbins,
                                   logbins=logbins, return_3d=return_3d)

        return xspec

    def freq_covariance(self):
        return find_modes.freq_covariance(self.map1, self.map2,
                                          self.noise_inv1,
                                          self.noise_inv2,
                                          self.freq, self.freq)
