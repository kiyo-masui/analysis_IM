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
import kiyopy.custom_exceptions as ce

# XXX
import matplotlib.pyplot as plt

# Constant that represents a very high noise level.  Setting the noise of a
# mode to this number deweights that mode.
T_infinity = 10000.0  # Kelvin**2

prefix ='dm_'
params_init = {               
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./testoutput_",
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               # In pixels.
               'map_shape' : (5, 5),
               'pixel_spacing' : 0.5, # degrees
               # How to treat the data.
               'polarizations' : ('I',)
               }

class DirtyMapMaker(object):
    """Dirty map maker.
    """

    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback



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
        
        # Store the full pointing matrix in sparse form.
        n_pointings = len(coords[0])
        n_coords = len(axis_names)
        if self._scheme == "nearest":
            self.dtype = int
        else:
            self.dtype = float
        # Loop over the time stream and get the weights for each pointing.
        memory_allocated = False
        for ii in xrange(n_pointings):
            coordinate = ()
            for jj in xrange(n_coords):
                coordinate += (self._coords[jj][ii],)
            try:
                pixels, weights = map.slice_interpolate_weights(
                    self._map_axes, coordinate, scheme)
            except ce.DataError:
                # This pointing is outside the map bounds.  When this happens,
                # point to the 0th pixel but with 0 weight.
                continue
            # On first iteration need to allocate memory for the sparse matrix
            # storage.
            if not memory_allocated:
                n_points_template = pixels.shape[0]
                self._pixel_inds = sp.zeros((n_pointings, n_coords,
                                              n_points_template), dtype=int)
                self._weights = sp.zeros((n_pointings, n_points_template),
                                         dtype=self.dtype)
                memory_allocated = True
            self._pixel_inds[ii,:,:] = pixels.transpose()
            self._weights[ii,:] = weights
    
    def get_sparse(self):
        """Return the arrays representing the pointing matrix in sparse form.

        Returns
        -------
        pixel_inds : array of ints
            Shape is (n_pointings, n_coordinates, n_pixels_per_pointing).
            These give which pixel entries are that are non zero for each
            pointing.  If a pointing is off the map, that row will be zeros.
        weights : array of floats
            Shape is (n_pointings, n_pixels_per_pointing).
            The entry of the pointing matrix corresponding with the matching
            pixel index.  If a pointing is off the map, that row will be zeros.
        """
        return self._pixel_inds, self._weights

    def apply_to_time_axis(self, time_stream, map_out=None):
        """Use this operator to convert a 'time' axis to a coordinate axis.
        
        This functions implements a fast matrix multiplication.  For input
        `map`, the following operations should be equivalent, with the later
        much more efficient.
        
        >>> a = al.partial_dot(self.get_matrix, map)
        >>> b = self.apply_to_time_axis(map)
        >>> sp.allclose(a, b)
        True
        
        """
        
        msg = "Fast multiply not implemented."
        raise NotImplementedError(msg)

    def get_matrix(self):
        """Gets the matrix representation of the pointing operator."""

        n_pointings = self._pixel_inds.shape[0]
        n_coords = self._pixel_inds.shape[1]
        n_pixels_per_pointing = self._pixel_inds.shape[2]
        # Initialize the output matrix.
        matrix = sp.zeros((n_pointings,) + self._map_shape, dtype=self.dtype)
        matrix = al.make_mat(matrix, axis_names=("time",) + self._axis_names,
                             row_axes=(0,), col_axes=range(1, n_coords + 1))
        # Loop over the time stream and get the weights for each pointing.
        for ii in xrange(n_pointings):
            for jj in xrange( n_pixels_per_pointing):
                matrix[(ii,) + tuple(self._pixel_inds[ii,:,jj])] = \
                        self._weights[ii, jj]
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

    # Internal nomiclature: The noise matrix is divided into three parts:  The
    # 'diagonal' contains one weight for every data point and represents a
    # fully diagonal matrix.  It has contributions from thermal noise and
    # deweights masked points.
    # The 'frequency_modes' part is the only part that
    # couples frequencies.  There are only a small number of modes along the
    # frequency axes, but each mode has a full time-time covariance matrix.
    # The modes are assumed to be uncorrelated.  The number of frequency modes
    # is named 'm'
    # The 'time_modes' part deweights certain modes along the time axis.  They
    # are uncorrelated between frequencies and the same for each frequency.
    # The number time modes is named 'q'.
    # The 'update' term will generally refer to the second term in the binomial
    # inverse identity, and the 'update_modes' refers to rotation matrices in
    # the identity.

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

    def add_time_modes(self, n_new_modes):
        """Initialize time modes to be deweighted at all frequencies.
        """

        self._assert_not_finalized()
        # TODO: Need to copy the axis info from self.info.
        if hasattr(self, "time_modes"):
            current_q = self.time_modes.shape[0]
            new_q = current_q + n_modes
            old_time_modes = self.time_modes
            time_modes = sp.zeros((new_q, self.n_time), 
                                    dtype=float)
            time_modes[:current_q,:] = old_allfreq
        else :
            current_q = 0
            time_modes = sp.zeros((n_new_modes, self.n_time), 
                               dtype=float)
        time_modes = al.make_mat(time_modes, axis_names=("time_mode", "time"), 
                                 row_axes=(0,), col_axes=(1,))
        self.time_modes = time_modes
        return current_q

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
        
        start = self.add_time_modes(1)
        self.time_modes[start,:] = 1.0/math.sqrt(self.n_time)

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
        freq_modes = sp.ones((1, self.n_chan), dtype=float)
        freq_modes /= sp.sqrt(self.n_chan)
        self.freq_modes = al.make_mat(freq_modes, 
                                      axis_names=("freq_mode", "freq"),
                                      row_axes=(0,), col_axes=(1,))
        # The covariance matrix for these modes (assumed orthoganal and
        # independant).
        freq_mode_noise = sp.empty((1, self.n_time, self.n_time),
                              dtype=float)
        self.freq_mode_noise = al.make_mat(freq_mode_noise, 
            axis_names=("freq_mode", "time", "time"), 
            row_axes=(0,1), col_axes=(0,2))
        # Build the matrix.
        time_deltas = abs(time[:, None] - time)
        # Smallest time step is the minimum of the first diagonal.
        dt = sp.amin(time_deltas.flat[1:len(time)+1])
        n_lags = sp.amax(time_deltas) // dt + 2
        correlation_function = noise_power.calculate_overf_correlation(amp, 
            index, f0, dt, n_lags)
        corr_func_interpolator = \
            interpolate.interp1d(sp.arange(n_lags)*dt, correlation_function)
        self.freq_mode_noise[0,...] = corr_func_interpolator(time_deltas)
        # Multiply by number of channels since this is the mean, not the sum
        # (Sept 15, 2011 in Kiyo's notes).
        self.freq_mode_noise *= self.n_chan

    def finalize(self):
        """Tell the class that you are done building the matrix.
        """
        
        n_time = self.n_time
        n_chan = self.n_chan
        # For now assume that all 3 noise components exist.  This should really
        # be adaptive so we don't have to have every component.
        if (not hasattr(self, 'freq_modes') or not hasattr(self, 'time_modes')
            or not hasattr(self, 'diagonal')):
            raise RuntimeError("Not all noise components have been set.")
        # Calculate the inverses of all matricies.
        freq_mode_inv = al.empty_like(self.freq_mode_noise)
        for ii in xrange(self.freq_mode_noise.shape[0]):
            freq_mode_inv[ii, ...] = \
                    linalg.inv(self.freq_mode_noise[ii, ...])
        self.freq_mode_inv = freq_mode_inv
        diagonal_inv = self.diagonal**-1
        self.diagonal_inv = al.as_alg_like(diagonal_inv, self.diagonal)
        # Calculate the term in the bracket in the matrix inversion lemma.
        # Get the size of the update term.
        # First, the rank of the correlated frequency part.
        m = self.freq_modes.shape[0]
        n_update =  m * n_time
        # Next, the rank of the all frequencies part.
        q = self.time_modes.shape[0]
        n_update += q * n_chan
        # Build the update matrix in blocks.
        freq_mode_update = sp.zeros((m, n_time, m, n_time), dtype=float)
        freq_mode_update = al.make_mat(freq_mode_update, 
            axis_names=('freq_mode', 'time', 'freq_mode', 'time'),
            row_axes=(0, 1), col_axes=(2, 3))
        cross_update = sp.zeros((m, n_time, q, n_chan), dtype=float)
        cross_update = al.make_mat(cross_update, 
            axis_names=('freq_mode', 'time', 'time_mode', 'freq'),
            row_axes=(0, 1), col_axes=(2, 3))
        time_mode_update = sp.zeros((q, n_chan, q, n_chan), dtype=float)
        time_mode_update = al.make_mat(time_mode_update, 
            axis_names=('time_mode', 'freq', 'time_mode', 'freq'),
            row_axes=(0, 1), col_axes=(2, 3))
        # Build the matrices.
        # Add the update mode noise in thier proper space.
        for ii in range(m):
            freq_mode_update[ii,:,ii,:] = self.freq_mode_inv[ii,:,:]
        time_mode_update.flat[::q * n_chan + 1] = 1.0/T_infinity
        # Now transform the diagonal noise to this funny space and add it to
        # the update term.
        # Do this one pair of modes at a time to make things less complicated.
        for ii in xrange(m):
            freq_mode1 = self.freq_modes.index_axis(0, ii)
            for jj in xrange(m):
                freq_mode2 = self.freq_modes.index_axis(0, jj)
                tmp_mat = al.partial_dot(freq_mode1, diagonal_inv)
                freq_mode_update[ii,:,jj,:].flat[::n_time + 1] += \
                        al.partial_dot(tmp_mat, freq_mode2)
        for ii in xrange(m):
            freq_mode = self.freq_modes.index_axis(0, ii)
            for jj in xrange(q):
                time_mode = self.time_modes.index_axis(0,jj)
                tmp_mat = al.partial_dot(freq_mode, diagonal_inv)
                cross_update[ii,:,jj,:] += al.partial_dot(tmp_mat, time_mode)
        for ii in xrange(q):
            time_mode1 = self.time_modes.index_axis(0, ii)
            for jj in xrange(q):
                time_mode2 = self.time_modes.index_axis(0,jj)
                tmp_mat = al.partial_dot(time_mode1, diagonal_inv)
                time_mode_update[ii,:,jj,:].flat[::n_chan + 1] += \
                        al.partial_dot(tmp_mat, time_mode2)
        # Put all the update terms in one big matrix and invert it.
        update_matrix = sp.empty((n_update, n_update), dtype=float)
        # Top left.
        update_matrix[:m * n_time,:m * n_time].flat[...] = \
            freq_mode_update.flat
        # Bottom right.
        update_matrix[m * n_time:,m * n_time:].flat[...] = \
            time_mode_update.flat
        # Top right.
        update_matrix[:m * n_time,m * n_time:].flat[...] = \
            cross_update.flat
        # Bottom left.
        tmp_mat = sp.swapaxes(cross_update, 0, 2)
        tmp_mat = sp.swapaxes(tmp_mat, 1, 3)
        update_matrix[m * n_time:,:m * n_time].flat[...] = \
            tmp_mat.flat
        update_matrix_inv = linalg.inv(update_matrix)
        # Copy the update terms back to thier own matrices and store them.
        freq_mode_update.flat[...] = \
                update_matrix_inv[:m * n_time,:m * n_time].flat
        self.freq_mode_update = freq_mode_update
        time_mode_update.flat[...] = \
                update_matrix_inv[m * n_time:,m * n_time:].flat
        self.time_mode_update = time_mode_update
        cross_update.flat[...] = \
                update_matrix_inv[:m * n_time,m * n_time:].flat
        self.cross_update = cross_update
        # Set flag so no more modifications to the matricies can occure.
        self._finalized = True
        
    # ---- Methods for using the Noise Matrix. ----

    def get_inverse(self):
        """Get the full noise inverse.

        This function is more for testing than accually being usefull (since in
        production we will only use part of the inverse at a time).
        """
        
        self._assert_finalized()
        n_chan = self.n_chan
        n_time = self.n_time
        freq_modes = self.freq_modes
        time_modes = self.time_modes
        # Get the size of the update term.
        # First, the rank of the correlated frequency part.
        m = self.freq_modes.shape[0]
        n_update =  m * n_time
        # Next, the rank of the all frequencies part.
        q = self.time_modes.shape[0]
        n_update += q * n_chan
        # Allowcate memory.
        out = sp.zeros((n_chan, n_time, n_chan, n_time), dtype=float)
        out = al.make_mat(out, axis_names=('freq', 'time', 'freq', 'time'),
                          row_axes=(0, 1), col_axes=(2, 3))
        # Loop over the frequency indeces to reduce workspace memory and
        # for simplicity.
        for ii in xrange(n_chan):
            this_freq_modes1 = freq_modes.index_axis(1, ii)
            for jj in xrange(n_chan):
                this_freq_modes2 = freq_modes.index_axis(1, jj)
                # The freq_mode-freq_mode part of the update term.
                tmp_mat = al.partial_dot(this_freq_modes1, 
                                         self.freq_mode_update)
                out[ii,:,jj,:] -= al.partial_dot(tmp_mat, this_freq_modes2)
                # The off diagonal blocks.
                this_cross1 = self.cross_update.index_axis(3, ii)
                this_cross2 = self.cross_update.index_axis(3, jj)
                tmp_mat = al.partial_dot(this_cross2, time_modes)
                out[ii,:,jj,:] -= al.partial_dot(this_freq_modes1, tmp_mat)
                tmp_mat = al.partial_dot(this_cross1, time_modes)
                tmp_mat = al.partial_dot(this_freq_modes2, tmp_mat)
                out[ii,:,jj,:] -= tmp_mat.transpose()
                # Finally the time_mode-time_mode part.
                this_time_update = self.time_mode_update.index_axis(1, ii)
                this_time_update = this_time_update.index_axis(2, jj)
                tmp_mat = al.partial_dot(time_modes.mat_transpose(), 
                                         this_time_update)
                out[ii,:,jj,:] -= al.partial_dot(tmp_mat, time_modes)
                # Multiply by thermal.
                out[ii,:,jj,:] *= self.diagonal_inv[ii,:,None]
                out[ii,:,jj,:] *= self.diagonal_inv[jj,None,:]
        # Add the thermal term.
        out.flat[::n_chan * n_time + 1] += self.diagonal_inv.flat[...]
        return out

    def noise_weight_time_stream(self, data):
        """Noise weight a timeleft_part_update_term.axes stream data vector.
        """
        
        time_modes = self.time_modes
        freq_modes = self.freq_modes
        freq_mode_update = self.freq_mode_update
        time_mode_update = self.time_mode_update
        cross_update = self.cross_update
        diagonal_inv = self.diagonal_inv
        # Noise weight by the diagonal part of the noise.
        diag_weighted = al.dot(diagonal_inv, data)
        # Calculate the update term carrying the freq modes and the time modes
        # through separately.
        # Transform to the update space.
        tmp_update_term_freq = al.partial_dot(freq_modes, diag_weighted)
        update_term_time = al.partial_dot(time_modes, diag_weighted)
        # Multiply by the update matrix.
        update_term_freq = (al.partial_dot(freq_mode_update, 
                                           tmp_update_term_freq)
                            + al.partial_dot(cross_update, update_term_time))
        update_term_time = (al.partial_dot(time_mode_update, update_term_time)
                            + al.partial_dot(cross_update.mat_transpose(),
                                             tmp_update_term_freq))
        # Transform back.
        update_term_freq = al.partial_dot(freq_modes.mat_transpose(),
                                          update_term_freq)
        update_term_time = al.partial_dot(time_modes.mat_transpose(),
                                          update_term_time)
        # Combine.
        update_term = update_term_freq + update_term_time.transpose()
        # Final multiply by the diagonal component.
        update_term = al.partial_dot(diagonal_inv, update_term)
        # Combine the terms.
        out = diag_weighted - update_term
        return out

    def set_pointing(self, P):
        """Sets the pointing for the for the noise.

        Calling this function lets the Noise class know how to transform to map
        space.
        """

        if len(P._map_shape) != 2:
            raise NotImplementedError("Only 2 D pointings implemented.")
        self.pointing_inds, self.pointing_weights = P.get_sparse()

    def update_map_noise(self, f_ind, pixel_inds, map_noise_inv):
        """Convert noise to map space.

        We fill build the matrix one row at a time for performance reasons.
        This is implemented without using the algebra interface for the arrays,
        again for performance.
        """
        
        # Unpack pixel inds (assumes 2D pointing).
        x_ind = pixel_inds[0]
        y_ind = pixel_inds[1]
        # Get some variables from self.
        time_modes = self.time_modes
        freq_modes = self.freq_modes
        pointing_inds = self.pointing_inds
        pointing_weights = self.pointing_weights
        n_time = self.n_time
        n_pix_per_pointing = pointing_inds.shape[-1]
        # Numbers of frequency modes and time modes.
        m = freq_modes.shape[0]
        q = time_modes.shape[0]
        # Isolate some the parts of variouse arrays that are relevant to the
        # current frequency row.
        this_chan_time_mode_update = self.time_mode_update[:,f_ind,...]
        this_chan_freq_modes = freq_modes[:,f_ind,...]
        this_chan_cross_update = self.cross_update[:,:,:,f_ind]
        # For cross update at this channel, we really want the transpose.
        this_chan_cross_update = sp.rollaxis(this_chan_cross_update, -1, 0)
        this_chan_cross_update = sp.copy(this_chan_cross_update)
        # Figure out which pointing indecies apply to this row of the matrix.
        row_pointings = []
        row_weights = []
        for ii in xrange(n_time):
            for jj in xrange(n_pix_per_pointing):
                if (pointing_inds[ii,0,jj] == x_ind 
                    and pointing_inds[ii,1,jj] == y_ind):
                    row_pointings.append(ii)
                    row_weights.append(pointing_weights[ii,jj])
        # If this data set never points to the pixel for this row, do nothing.
        if not row_pointings:
            return
        # Now loop over the pointings that we just found.
        for ii in xrange(len(row_pointings)):
            time_ind = row_pointings[ii]
            weight = row_weights[ii]
            # Isolate the relevant parts of some of the matrices and vectors.
            this_time_time_modes = time_modes[:,time_ind]
            this_time_freq_mode_update = self.freq_mode_update[:,time_ind,...]
            this_time_cross_update = self.cross_update[:,time_ind,...]
            # Now build up the central matrix.




