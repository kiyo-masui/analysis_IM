"""Dirty map making module.

Module converts data in the time domain into noise weighted data in the map
domain, i.e. it creats the dirty map.  Module also contains many utilities like
the pointing operator (`Pointing`) and the time domain noise operator
(`Noise`).
"""

import math

import scipy as sp
import numpy.ma as ma
import scipy.fftpack as fft
from scipy import linalg
from scipy import interpolate

import core.algebra as al
import tools
from noise import noise_power

# XXX
import matplotlib.pyplot as plt

# Constant that represents a very high noise level.  Setting the noise of a
# mode to this number deweights that mode.
T_infinity = 10000.0  # Kelvin**2

class Pointing(object):
    """Class represents the pointing operator.

    The pointing operator converts from the map domain to the time domain in
    its native form and from the time domain to the map domain in its
    transposed form.
    
    Parameters
    ----------
    axis_names: tuple of strings
        The names of the axes in the map domain e.g. ("ra", "dec")
    coords : tuple of 1D arrays
        Tuple must be same length as `axis_names`.  The coordinates as a 
        function of time for each of the map axes.
    map : al.vect object
        The map that we will be gridding onto.  Map axes must include
        `axis_names`.  No modification to the map is made, only the axis
        information is used.
    scheme : string
        Gridding scheme to use.  Choices are 'nearest'.
    """

    def __init__(self, axis_names, coords, map, scheme='nearest'):
        # Sanity check some of the inputs.
        if len(axis_names) != len(coords):
            msg = "Got %d pointing axis names, but got %d coordinate arrays."
            raise ValueError(msg % (len(axis_names), len(coords)))
        n_pointings = coords[0].shape[0]
        for coordinate_array in coords:
            if coordinate_array.shape != (n_pointings,):
                msg = "Coordinate arrays must all be 1D and same length."
                raise ValueError(msg)
        # Save the data we need for the pointing.
        self._axis_names = axis_names
        self._coords = tuple(coords) # A tuple of 1D arrays.
        self._scheme = scheme
        # Get all coordinate information from the map.
        # Figure out which axes we are doing the pointing for, get the shape of
        # those axes, the coordinate centres and the pixel size.
        map_axis_indices = ()
        map_axes_shape = ()
        #map_axes_centre = ()
        #map_axes_delta = ()
        list_map_axis_names = list(map.axes)
        for axis_name in axis_names:
            axis_index = list_map_axis_names.index(axis_name)
            map_axis_indices += (axis_index,)
            map_axes_shape += (map.shape[axis_index],)
            #map_axes_centre += (map.info[axis_name + "_centre"],)
            #map_axes_delta += (map.info[axis_name + "_delta"],)
        self._map_axes = map_axis_indices
        self._map_shape = map_axes_shape
        #self.map_delta = map_axes_delta
        #self.map_centre = map_axes_centre
        # Store a reference to the map.  We are storing a copy of the map that
        # has its own copy of all the metadata (including the numpy array
        # metadata, such as .shape) so we don't have to worry about changes to
        # the object whose reference we store.  However, we don't want to copy
        # the acctual data, since this class only uses the metadata.
        self._map = al.as_alg_like(map[...], map)

    def apply_to_time_axis(self, time_stream, map_out=None):
        """Use this operator to convert a 'time' axis to a coordinate axis.
        
        This functions implements a fast matrix multiplication.  For input
        `map`, the following operations should be equivalent, with the later
        much more efficient.
        
        >>> a = al.partial_dot(self.get_matrix.transpose(), map)
        >>> b = self.apply_to_time_axis(map)
        >>> sp.allclose(a, b)
        True
        
        """
        
        msg = "Fast multiply not implemented."
        raise NotImplementedError(msg)

    def get_matrix(self):
        """Gets the matrix representation of the pointing operator."""

        n_pointings = len(self._coords[0])
        n_coords = len(self._axis_names)
        if self._scheme == "nearest":
            dtype = int
        else:
            dtype = float
        # Initialize the output matrix.
        matrix = sp.zeros((n_pointings,) + self._map_shape, dtype=dtype)
        matrix = al.make_mat(matrix, axis_names=("time",) + self._axis_names,
                             row_axes=(0,), col_axes=range(1, n_coords + 1))
        # Loop over the time stream and get the weights for each pointing.
        for ii in xrange(n_pointings):
            coordinate = ()
            for jj in xrange(n_coords):
                coordinate += (self._coords[jj][ii],)
            pixels, weights = self._map.slice_interpolate_weights(
                self._map_axes, coordinate, self._scheme)
            index = (ii,)
            for jj in xrange(n_coords):
                index += (pixels[:,jj],)
            matrix[index] = weights
        return matrix


class Noise(object):
    """Object that represents the noise matrix for time stream data.
    
    The noise matrix is represented as separate components each with different
    symetries.  This is so the huge matrix does not have to be stored.

    Parameters
    ----------
    time_strem_data : al.vect object
        The data for which we want to represent the noise matrix.  Only meta
        data is used, not the acctual data.
    time : 1D array
        The time axis of the data.
    """

    # ---- Initialization methods. ----

    def __init__(self, time_stream_data, time):
        if len(time_stream_data.shape) != 2:
            raise ValueError("Only 2D data suported (freq, time).")
        self.n_chan = time_stream_data.shape[0]
        self.n_time = time_stream_data.shape[1]
        self.info = dict(time_stream_data.info)
        self._finalized = False
        self.time = time

    def _assert_not_finalized(self):
        """Make sure the noise matrix is not finalized and can be modified."""
        if self._finalized:
            raise AssertionError("Noise model closed for modification.")
    
    def _assert_finalized(self):
        """Make sure the noise matrix is finalized and can not be modified."""
        if not self._finalized:
            raise AssertionError("Noise model still being modified.")

    def initialize_diagonal(self):
        """Create the diagonal part of the noise matrix if it doesn't exist."""
        
        self._assert_not_finalized()
        # TODO: Need to copy the axis info from self.info.
        if hasattr(self, "diagonal"):
            return
        diagonal = sp.zeros((self.n_chan, self.n_time), dtype=float)
        diagonal = al.make_mat(diagonal, axis_names=("freq", "time"), 
                               row_axes=(0, 1), col_axes=(0, 1))
        self.diagonal = diagonal

    def initialize_allfreq(self):
        """Create the part of the noise matrix that contributes to a
        frequencies equally.
        """

        self._assert_not_finalized()
        # TODO: Need to copy the axis info from self.info.
        if hasattr(self, "allfreq"):
            return
        allfreq = sp.zeros((self.n_time, self.n_time), 
                           dtype=float)
        allfreq = al.make_mat(allfreq, axis_names=("time", "time"), 
                               row_axes=(0,), col_axes=(1,))
        self.allfreq = allfreq

    # ---- Methods that build up the noise matrix. ----

    def add_thermal(self, thermal_levels):
        """Add a thermal component to the noise.
        
        This modifies the diagonal part of the noise matrix.

        Parameters
        ----------
        thermal_levels : 1D array
            The noise level of each channel (frequency)
        """
        
        self.initialize_diagonal()
        self.diagonal += thermal_levels[:, None]

    def add_mask(self, mask_inds):
        """Add a mask to the noise.
        
        This modifies the diagonal part of the noise matrix.

        Parameters
        ----------
        mask_inds : 2 element tuple of integer arrays.
            The locations of data points to be masked.
        """
        
        self.initialize_diagonal()
        self.diagonal[mask_inds] += T_infinity

    def deweight_time_mean(self):
        """Deweights time mean in each channel.

        This modifies the part of the noise matrix that is the same for each
        channel.
        """
        
        self.initialize_allfreq()
        nt = self.n_time
        new_noise = sp.zeros((nt, nt), dtype=float)
        # Formula from Sept 7, 2011 of Kiyo's notes.
        new_noise += (T_infinity - 1.0)/nt
        new_noise.flat[::nt + 1] += 1.0
        self.allfreq += new_noise

    def deweight_time_slope(self):
        """Deweights time slope in each channel.
        
        This modifies the part of the noise matrix that is the same for each
        channel.
        """

        raise NotImplementedError()
    
    def add_correlated_over_f(self, amp, index, f0):
        """Add 1/f noise that is perfectly correlated between frequencies.

        This modifies the the correlated mode part of the noise matrix.

        Parameters
        ----------
        amp : float
            The amplitude of the noise at `f0`.
        index : float
            The spectral index of the spectrum (normaly near -1).
        f0 : float
            The pivot frequency.
        """
        
        time = self.time
        self._assert_not_finalized()
        if hasattr(self, "correlated_modes"):
            raise RuntimeError("The correlated modes already exist and would "
                               "be overwritten.")
        # Only have one channel mode, the mean mode.
        freq_modes = sp.ones((self.n_chan, 1), dtype=float)
        freq_modes /= sp.sqrt(self.n_chan)
        self.freq_modes = al.make_mat(freq_modes, axis_names=("freq", "mode"),
                                      row_axes=(0,), col_axes=(1,))
        # The covariance matrix for these modes (assumed orthoganal and
        # independant).
        mode_noise = sp.empty((1, self.n_time, self.n_time),
                              dtype=float)
        self.mode_noise = al.make_mat(mode_noise, axis_names=("mode", "time",
            "time"), row_axes=(0,1), col_axes=(0,2))
        # Build the matrix.
        time_deltas = abs(time[:, None] - time)
        # Smallest time step is the minimum of the first diagonal.
        dt = sp.amin(time_deltas.flat[1:len(time)+1])
        n_lags = sp.amax(time_deltas) // dt + 2
        correlation_function = noise_power.calculate_overf_correlation(amp, 
            index, f0, dt, n_lags)
        correlation_function_interpolator = \
            interpolate.interp1d(sp.arange(n_lags)*dt, correlation_function)
        self.mode_noise[0,...] = correlation_function_interpolator(time_deltas)
        # Multiply by number of channels since this is the mean, not the sum
        # (Sept 15, 2011 in Kiyo's notes).
        self.mode_noise *= self.n_chan

    def finalize(self):
        """Tell the class that you are done building the matrix.
        """
        
        # For now assume that all 3 noise components exist.  This should really
        # be adaptive so we don't have to have every component.
        if (not hasattr(self, 'mode_noise') or not hasattr(self, 'allfreq') or 
            not hasattr(self, 'diagonal')):
            raise RuntimeError("Not all noise components have been set.")
        # Calculate the inverses of all matricies.
        self.mode_noise_inv = al.empty_like(self.mode_noise)
        for ii in xrange(self.mode_noise.shape[0]):
            self.mode_noise_inv[ii, ...] = linalg.inv(self.mode_noise[ii, ...])
        allfreq_inv = linalg.inv(self.allfreq)
        self.allfreq_inv = al.as_alg_like(allfreq_inv, self.allfreq)
        diagonal_inv = self.diagonal**-1
        self.diagonal_inv = al.as_alg_like(diagonal_inv, self.diagonal)
        # Calculate the term in the bracket in the matrix inversion lemma.
        self.mode_noise_term = al.empty_like(self.mode_noise_inv)
        # TODO: Probably improve performance by interchanging the loops.  For
        # now it doesn't matter since there is only 1 mode.
        for ii in xrange(self.mode_noise.shape[0]):
            this_mode_noise_inv = self.mode_noise_inv[ii,...]
            this_mode = sp.copy(self.freq_modes[:,ii])
            tmp_mat = sp.zeros_like(this_mode_noise_inv)
            for jj in xrange(self.n_chan):
                tmp_mat += this_mode[jj]**2 * self.get_diag_allfreq_inverse(jj)
            tmp_mat += this_mode_noise_inv
            self.mode_noise_term[ii,...] = linalg.inv(tmp_mat, True)
        # Set flag so no more modifications to the matricies can occure.
        self.finalized = True
    
    # ---- Methods for using the Noise Matrix. ----

    def get_diag_allfreq_inverse(self, find):
        """Get the inverse of the combined diagonal and all frequency noise.
        
        Get the noise at only one frequency slice given by index `find`.
        """
        
        if find == "full":
            # Get the full matrix, not just one frequency slice.
            out = sp.empty((self.n_chan, self.n_time, self.n_time),
                           dtype=float)
            out = al.make_mat(out, axis_names=("freq", "time", "time"),
                              row_axes=(0, 1), col_axes=(0, 2))
            out.copy_axis_info(self.info)
            for ii in xrange(self.n_chan):
                out[ii, ...] = self.get_diag_allfreq_inverse(ii)
        else:
            # Get the data we need.
            diag = self.diagonal[find,...]
            allfreq = self.allfreq
            # Add the matricies.
            out = sp.copy(allfreq)
            out.flat[::self.n_time + 1] += diag
            # Invert it..
            out = linalg.inv(out, True)
            # Set up the output as an algebra object.
            out = al.as_alg_like(out, allfreq)
        return out

    def get_inverse(self):
        """Get the full noise inverse.

        This function is more for testing than accually being usefull (since in
        production we will only use part of the inverse at a time).
        """
        
        n_chan = self.n_chan
        n_time = self.n_time
        # Allowcate memory
        out = sp.zeros((n_chan, n_time, n_chan, n_time), dtype=float)
        # Loop over all frequency pairs.
        for ii in xrange(n_chan):
            mat1 = self.get_diag_allfreq_inverse(ii)
            for jj in xrange(n_chan):
                mat2 = self.get_diag_allfreq_inverse(jj)
                # Loop over modes.
                for kk in xrange(self.mode_noise_term.shape[0]):
                    tmp_mat = sp.copy(self.mode_noise_term[kk,...])
                    tmp_mat *= -self.freq_modes[ii,kk]*self.freq_modes[jj,kk]
                    if kk == 0:
                        mat_block = tmp_mat
                    else:
                        mat_block += tmp_mat
                mat_block = sp.dot(mat1, mat_block)
                mat_block = sp.dot(mat_block, mat2)
                if ii == jj:
                    mat_block += mat1
                out[ii,:,jj,:] = mat_block
        out = al.make_mat(out, axis_names=('freq', 'time', 'freq', 'time'),
                          row_axes=(0, 1), col_axes=(2, 3))
        return out

    def noise_weight_time_stream(self, data):
        """Noise weight a time stream data vector.
        """
        
        # Get the inverse of the 'independant frequency' part of the noise.
        block_freq_inverse = self.get_diag_allfreq_inverse("full")
        # Noise weight by the 'independant frequency' term.
        block_weighted = al.dot(block_freq_inverse, data)
        # Calculate the second term that couples frequencies.
        second_term = al.partial_dot(self.freq_modes.transpose(),
                                     block_weighted)
        second_term = al.dot(self.mode_noise_term, second_term)
        second_term = al.partial_dot(self.freq_modes, second_term)
        second_term = al.dot(block_freq_inverse, second_term)
        # Combine the terms.
        out = block_weighted - second_term

        return out

        




