r"""Program that calculates the correlation function across frequency slices.
"""

import time
import sys
import scipy as sp
import numpy as np
import numpy.ma as ma
from core import algebra
from map import beam
from map import physical_gridding as pg
from kiyopy import parse_ini
import kiyopy.utils
from correlate import corr_estimation
from correlate import pwrspec_estimation as pe
from utils import data_paths
from utils import batch_handler as bh
# TODO: move single map operations to a separate class

params_init = {
               'output_root': "./data_test/",
               'output_filetag': 'test',
               'map1': "map1.npy",
               'map2': "map2.npy",
               'noise_inv1': "noise_inv1.npy",
               'noise_inv2': "noise_inv2.npy",
               'freq_list': (),
               # Angular lags at which to calculate the correlation.  Upper
               # edge bins in degrees.
               'lags': tuple(sp.arange(0.002, 0.2, 0.12))
               }
prefix = 'fs_'


class MapPair(object):
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

    def __init__(self, map1, map2, noise_inv1, noise_inv2, freq, avoid_db=False):
        r"""
        arguments: map1, map2, noise_inv1, noise_inv2, freq
        """
        # just take filenames to these objects to avoid independent db lookups
        # which get expensive to do for e.g. 10000 quadratic products
        if avoid_db:
            self.map1 = algebra.make_vect(algebra.load(map1))
            self.map2 = algebra.make_vect(algebra.load(map2))
            self.noise_inv1 = algebra.make_vect(algebra.load(noise_inv1))
            self.noise_inv2 = algebra.make_vect(algebra.load(noise_inv2))
        else:
            self.datapath_db = data_paths.DataPath()
            self.map1 = self.datapath_db.fetch_multi(map1)
            self.map2 = self.datapath_db.fetch_multi(map2)
            self.noise_inv1 = self.datapath_db.fetch_multi(noise_inv1)
            self.noise_inv2 = self.datapath_db.fetch_multi(noise_inv2)

        self.freq = freq

        # set the physical-dimension maps to None
        self.phys_map1 = None
        self.phys_map2 = None
        self.phys_noise_inv1 = None
        self.phys_noise_inv2 = None

        # Give infinite noise to unconsidered frequencies
        self.sanitize()

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

    def make_physical(self, refinement=1, pad=5, order=2):
        r"""Project the maps and weights into physical coordinates
        refinement allows the physical bins to be e.g. =2 times finer
        pad puts a number of padded pixels on all sides of the physical vol.
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

    def make_noise_factorizable(self):
        r"""Convert noise weights such that the factor into a function a
        frequency times a function of pixel by taking means over the original
        weights.
        """
        print "making the noise factorizable"

        def make_factorizable(noise):
            r"""factorize the noise"""
            noise[noise < 1.e-30] = 1.e-30
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

    def subtract_frequency_modes(self, modes1, modes2=None):
        r"""Subtract frequency mode from the map.

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

    # TODO: add documentation
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
        # TODO: possibly revert to old correlation function calc?
        #return corr_estimation.corr_est(self.map1, self.map2,
        #                                self.noise_inv1, self.noise_inv2,
        #                                self.freq, self.freq,
        #                                lags=lags, speedup=speedup,
        #                                verbose=verbose)

        return corr_estimation.freq_covariance(self.map1, self.map2,
                                               self.noise_inv1,
                                               self.noise_inv2,
                                               self.freq, self.freq)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        CorrelateSingle(str(sys.argv[1])).execute()
    else:
        print 'Need one argument: parameter file name.'
