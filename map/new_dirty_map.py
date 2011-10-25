"""Dirty map making module.

Module converts data in the time domain into noise weighted data in the map
domain, i.e. it creats the dirty map.  Module also contains many utilities like
the pointing operator (`Pointing`) and the time domain noise operator
(`Noise`).
"""

import math
import threading
from Queue import Queue

import scipy as sp
import numpy.ma as ma
import scipy.fftpack as fft
from scipy import linalg
from scipy import interpolate
import numpy as np

import core.algebra as al
import tools
from noise import noise_power
import kiyopy.custom_exceptions as ce
from kiyopy import parse_ini
import _mapmaker as _mapmaker_c

# XXX
import matplotlib.pyplot as plt

# Constant that represents a very high noise level.  Setting the noise of a
# mode to this number deweights that mode.
# 100 K**2 seems like a good number.  Bigger than any actual noise mode we would
# encounter, but not so big that it screws up numerical stability.
T_infinity = 100.0  # Kelvin**2
#T_infinity = 10.  # Kelvin**2

prefix ='dm_'
params_init = {               
               # IO:
               'input_root' : './',
               # The unique part of every fname
               'file_middles' : ("testfile_GBTfits",),
               'input_end' : ".fits",
               'output_root' : "./testoutput_",
               # What data to include from each file.
               'scans' : (),
               'bands' : (0,),
               'polarizations' : ('I',),
               # Map parameters (Ra (deg), Dec (deg)).
               'field_centre' : (325.0, 0.0),
               # In pixels.
               'map_shape' : (5, 5),
               'pixel_spacing' : 0.5, # degrees
               # How to treat the data.
               # How much data to include at a time (scan by scan or file by
               # file)
               'time_block' : 'file',
               # What kind of frequency correlations to use.  Options are
               # 'None', 'mean' and 'measured'.
               'frequency_correlations' : 'None',
               # Ignored unless 'frequency_correlations' is measured.
               'number_frequency_modes' : 1,
               # What number to use as the thermal nose.  Options are
               # 'channel_var', and 'thermal'.
               'thermal_weight' : 'thermal',
               'deweight_time_mean' : True,
               'deweight_time_slope' : False,
               }

class DirtyMap(object):
    """Dirty map maker.
    """

    def __init__(self, parameter_file_or_dict=None, feedback=2):
        # Read in the parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                 prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def preprocess_data(self):
        """Reads data and does some preprocessing and data selection."""
        pass

    def get_noise_parameter(self, parameter_name):
        """Reads a desired noise parameter for the current data."""
        pass

    def execute(self, n_processes):
        """Driver method."""
        
        self.n_processes = n_processes
        params = self.params
        n_pols = len(params["polarizations"])
        n_bands = len(params["bands"])
        for ii in range(n_pols):
            for jj in range(n_bands):
                self.pol = ii
                self.band = jj
                # Initialization of the outputs.
                # XXX
                self.n_chan = n_chan
                self.n_ra = n_ra
                self.n_dec = n_dec
                self.map = map
                self.cov_inv = cov_inv
                # XXX
                # Do work.
                self.make_map()
                # IO.

    def make_map(self):
        """Makes map for current polarization and band.
        
        This worker function has been split off for testing reasons.
        """
        
        params = self.params
        map = self.map
        cov_inv = self.cov_inv
        n_files = len(params["file_middles"])
        # Initialize lists for the data and the noise.
        data_list = []
        noise_list = []
        pointing_list = []
        # Loop over the files get the data and build the noise.
        # TODO change this to an iterator over blocks (scans or files).
        for ii in xrange(n_files):
            # Setting the file number tells the IO methods which data to read.
            self.file_number = ii
            # Get all the information we need from the data files.
            # ### Divid by 0 error here.
            time_stream, ra, dec, az, el, time, mask_inds = \
                    self.preprocess_data()
            P = Pointing(("ra", "dec"), (ra, dec), map, "linear")
            # Now build up our noise model for this piece of data.
            N = Noise(time_stream, time)
            N.add_mask(mask_inds)
            # The thermal part.
            if params["thermal_weight"] is "thermal":
                thermal_noise = self.get_noise_parameter("thermal")
                N.add_thermal(thermal_noise)
            else:
                raise ValueError("Invalid 'thermal_weight'.")
            # Frequency correlations.
            if params['frequency_correlations'] is 'mean':
                mean_overf = self.get_noise_parameter("mean_over_f")
                N.add_correlated_over_f(*mean_overf)
            else:
                raise ValueError("Invalid frequency correlations.")
            # Things to do along the time axis.
            if params['deweight_time_mean']:
                N.deweight_time_mean()
            if params['deweight_time_slope']:
                N.deweight_time_slope()
            N.finalize()
            # Make the dirty map.
            weighted_time_stream = N.weight_time_stream(time_stream)
            map += P.apply_to_time_axis(weighted_time_stream)
            # Store all these for later.
            data_list.append(time_stream)
            pointing_list.append(P)
            noise_list.append(N)
        # Now we have lists with all the data and thier noise.  Accumulate it
        # into a dirty map and its covariance.
        n_time_blocks = len(data_list)
        # Initialize the queue of work to be done.
        index_queue = Queue()
        # Define a function that takes work off a queue and does it.  Variables
        # local to this function prefixed with 'thread_'.
        def thread_work():
            while True:
                thread_inds = index_queue.get()
                # None will be the flag that there is no more work to do.
                if thread_inds is None:
                    return
                thread_f_ind = thread_inds[0]
                thread_ra_ind = thread_inds[1]
                thread_cov_inv_row = sp.zeros((self.n_dec, self.n_chan,
                                               self.n_ra, self.n_dec),
                                              dtype=float)
                for thread_kk in xrange(n_time_blocks):
                    thread_P = pointing_list[thread_kk]
                    thread_N = noise_list[thread_kk]
                    thread_P.noise_to_map_domain(thread_N, thread_f_ind, 
                                    thread_ra_ind, thread_cov_inv_row)
                cov_inv[thread_f_ind,thread_ra_ind,...] = thread_cov_inv_row
        # Start the worker threads.
        thread_list = []
        for ii in range(self.n_processes):
            T = threading.Thread(target=thread_work)
            T.start()
            thread_list.append(T)
        # Now put work on the queue for the threads to do.
        for ii in xrange(self.n_chan):
            for jj in xrange(self.n_ra):
                index_queue.put((ii, jj))
        # At the end of the queue, tell the threads that they are done.
        for ii in range(self.n_processes):
            index_queue.put(None)
        # Wait for the threads.
        for T in thread_list:
            T.join()
        if not index_queue.empty():
            raise RuntimeError("A thread had an error.")
        # Now go through and make sure that the noise isn't singular by adding
        # a bit of information to untouched pixels.
        diag_slice = slice(None, None, 
                           self.n_ra * self.n_dec * self.n_chan + 1)
        cov_diag = cov_inv.flat[diag_slice]
        untouched_inds = cov_diag < 1.0e-8 / T_infinity
        tmp_add = sp.zeros(cov_diag.size, dtype=float)
        tmp_add[untouched_inds] = 1.0 / T_infinity
        cov_inv.flat[diag_slice] += tmp_add


#### Classes ####

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
        self.dtype = np.float
        # Loop over the time stream and get the weights for each pointing.
        memory_allocated = False
        for ii in xrange(n_pointings):
            coordinate = ()
            for jj in xrange(n_coords):
                coordinate += (self._coords[jj][ii],)
            pixels, weights = map.slice_interpolate_weights(
                self._map_axes, coordinate, scheme)
            # On first iteration need to allocate memory for the sparse matrix
            # storage.
            if not memory_allocated:
                n_points_template = pixels.shape[0]
                self._pixel_inds = sp.zeros((n_pointings, n_coords,
                                              n_points_template), dtype=np.int)
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
        
        This functions implements a fast matrix multiplication. It is roughly
        equivalent to using `algebra.partial_dot` except in axis placement.
        This function "replaces" the time axis with the map axes which is
        different from `partial_dot`'s behaviour.  This function is much more
        efficient than using `partial_dot`.
        
        For input
        `map`, the following operations should be equivalent, with the later
        much more efficient.
        
        # XXX This example is wrong and broken.
        >>> a = al.partial_dot(map)
        >>> b = self.apply_to_time_axis(map)
        >>> sp.allclose(a, b)
        True
        """
        
        if not isinstance(time_stream, al.vect):
            raise TypeError("Input data must be an algebra.vect object.")
        # Find the time axis for the input.
        for ii in range(time_stream.ndim):
            if time_stream.axes[ii] == 'time':
                time_axis = ii
                break
        else :
            raise ValueError("Input data vect doesn't have a time axis.")
        # Get some dimensions.
        n_pointings = self._pixel_inds.shape[0]
        n_pixels_template = self._pixel_inds.shape[2]
        n_axes = len(self._axis_names)
        if time_stream.shape[time_axis] != n_pointings:
            msg = ("Time stream data and pointing have different number of"
                   " time points.")
            raise ValueError(msg)
        # Get the shape and axis names of the output.
        out_shape = (time_stream.shape[:time_axis] + self._map_shape
                     + time_stream.shape[time_axis + 1:])
        out_axes = (time_stream.axes[:time_axis] + self._axis_names
                     + time_stream.axes[time_axis + 1:])
        # Allowcate output memory if not passed.
        if map_out is None:
            map_out = sp.zeros(out_shape, dtype=float)
            map_out = al.make_vect(map_out, axis_names=out_axes)
        else :
            if map_out.shape != out_shape:
                raise ValueError("Output array is the wrong shape.")
        # Initialize tuples that will index the input and the output.
        data_index = [slice(None),] * time_stream.ndim + [None]
        out_index = [slice(None),] * map_out.ndim
        # Loop over the time axis and do the dot.
        for ii in xrange(n_pointings):
            data_index[time_axis] = ii
            for kk in xrange(n_axes):
                out_index[time_axis + kk] = self._pixel_inds[ii,kk,:]
            map_out[tuple(out_index)] += (self._weights[ii,:]
                                   * time_stream[tuple(data_index)])
        return map_out

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

    def noise_to_map_domain(self, Noise, f_ind, ra_ind, map_noise_inv):
        """Convert noise to map space.
        
        For performace and IO reasons this is done with a call to this function
        for each frequency row and each ra row.  All dec rows and all columns
        are handled in this function simultaniousely.

        This function is designed to be thread safe in that if it is called
        from two separate threads but with different `f_ind` or `ra_ind`, there
        should be no race conditions.
        """
        
        Noise._assert_finalized()
        _mapmaker_c.update_map_noise_chan_ra_row(Noise.diagonal_inv,
                Noise.freq_modes, Noise.time_modes, Noise.freq_mode_update,
                Noise.time_mode_update, Noise.cross_update, self._pixel_inds,
                self._weights, f_ind, ra_ind, map_noise_inv)
        

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
            new_q = current_q + n_new_modes
            old_time_modes = self.time_modes
            time_modes = sp.zeros((new_q, self.n_time), 
                                    dtype=float)
            time_modes[:current_q,:] = old_time_modes
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
        # Purposely not normalized to better represent the 'mean mode'.
        # XXX To make this normalized, need to make a time_modes_noise matrix
        # and put the extra factors in it (like in time_modes_noise).
        self.time_modes[start,:] = 1.0

    def deweight_time_slope(self):
        """Deweights time slope in each channel.
        
        This modifies the part of the noise matrix that is the same for each
        channel.
        """

        start = self.add_time_modes(1)
        mode = self.time - sp.mean(self.time)
        # Purposely normalized to sqrt(n).
        mode *= sp.sqrt(self.n_time) / sp.sqrt(sp.sum(mode**2)) 
        self.time_modes[start,:] = mode

    def add_correlated_over_f(self, amp, index, f0):
        """Add 1/f noise that is perfectly correlated between frequencies.

        This modifies the the correlated mode part of the noise matrix. It adds
        a single mode with equal amplitude in all frequencies.

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
        # TODO: Som sort of condition number check here?
        #if hasattr(self, 'flag'):
        #    e, v = linalg.eigh(update_matrix)
        #    print "Update term condition number:", max(e)/min(e)
        #    e, v = linalg.eigh(update_matrix_inv)
        #    print "Update term  inverse condition number:", max(e)/min(e)
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
                # XXX
                #if hasattr(self, 'flag'):
                #    e, v = linalg.eigh(out[ii,:,jj,:])
                #    print "Intermediate 1 eig range:", max(-e), min(-e)
                ###
                # Multiply by thermal.
                out[ii,:,jj,:] *= self.diagonal_inv[ii,:,None]
                out[ii,:,jj,:] *= self.diagonal_inv[jj,None,:]
                # XXX
                #if hasattr(self, 'flag'):
                #    e, v = linalg.eigh(out[ii,:,jj,:])
                #    print "Intermediate 2 eig range:", max(-e), min(-e),
                #    print max(-e)/min(-e)
                #    print -e
                ###
        # Add the thermal term.
        out.flat[::n_chan * n_time + 1] += self.diagonal_inv.flat[...]
        # XXX
        #if hasattr(self, 'flag'):
        #    a = out.flat[::n_chan * n_time + 1] / self.diagonal_inv.flat[...]
        #    print "subtraction precision loss:", max(a), min(a)
        return out

    def weight_time_stream(self, data):
        """Noise weight a timeleft_part_update_term.axes stream data vector.
        """
        
        self._assert_finalized()
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


#### Utilities ####

def trim_time_stream(data, coords, lower_bounds, upper_bounds):
    """Discards data that is out of the map bounds.

    Parameters
    ----------
    data : algebra.vect
        The time stream data.  The vect object axes must be ('freq', 'time').
    coords : tuple (length `n`) of arrays
        The pointing information for `n` map coordinates as a function of time.
    lower_bounds : tuple (length `n`) of floats.
        Lower map boundaries for each coordinate.
    upper_bounds : tuple (length `n`) of floats.
        Upper map boundaries for each coordinate.

    Returns
    -------
    trimmed_data : algebra.vect
        New data vector with ashortened time axes to exclude off map
        pointings.
    inds : 1d array of ints
        Which tim stream array indices were retained (this array can be used 
        index other arrays).
    """
    
    n_time = data.shape[1]
    # Loop through the pointing information and see which times are on the map.
    inds = []
    for ii in range(n_time):
        for jj in range(len(coords)):
            if (coords[jj][ii] < lower_bounds[jj]
                or coords[jj][ii] > upper_bounds[jj]):
                break
        else:
            inds.append(ii)
    inds = sp.asarray(inds, dtype=int)
    # Shorten the time axis and cast as a vect object.
    trimmed_data = sp.ascontiguousarray(data[:,inds])
    trimmed_data = al.make_vect(trimmed_data, axis_names=('freq', 'time'))
    return trimmed_data, inds
