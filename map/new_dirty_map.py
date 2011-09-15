"""Dirty map making module."""

import math

import scipy as sp
import numpy.ma as ma
import scipy.fftpack as fft

import core.algebra as al
import tools

T_infinity = 10000.0  # Kelvin**2

class Pointing(object):
    """Class represents the pointing operator."""

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
    """

    def __init__(self, time_stream_data):
        if len(time_stream_data.shape) != 2:
            raise ValueError("Only 2D data suported (freq, time).")
        self.data_shape = time_stream_data.shape

    def initialize_diagonal(self):
        """Create the diagonal part of the noise matrix if it doesn't exist."""

        if hasattr(self, "diagonal"):
            return
        diagonal = sp.zeros(self.data_shape, dtype=float)
        diagonal = al.make_mat(diagonal, axis_names=("freq", "time"), 
                               row_axes=(0, 1), col_axes=(0, 1))
        self.diagonal = diagonal

    def initialize_allfreq(self):
        """Create the part of the noise matrix that contributes to a
        frequencies equally.
        """

        if hasattr(self, "allfreq"):
            return
        allfreq = sp.zeros((self.data_shape[1], self.data_shape[1]), 
                           dtype=float)
        allfreq = al.make_mat(allfreq, axis_names=("time", "time"), 
                               row_axes=(0,), col_axes=(1,))
        self.allfreq = allfreq


    def add_thermal(self, thermal_levels):
        """Add a thermal component to the noise.
        """
        
        self.initialize_diagonal()
        self.diagonal += thermal_levels[:, None]

    def add_mask(self, mask_inds):
        """Add a mask to the noise.
        """
        
        self.initialize_diagonal()
        self.diagonal[mask_inds] += T_infinity

    def deweight_time_mean(self):
        """Deweights time mean in each channel.
        """
        
        self.initialize_allfreq()
        nt = self.data_shape[1]
        new_noise = sp.zeros((nt, nt), dtype=float)
        new_noise += T_infinity/nt
        new_noise.flat[::nt + 1] += 1.0
        self.allfreq += new_noise

    def deweight_time_slope(self):
        """Deweights time slope in each channel.
        """

        raise NotImplementedError()
    
    def add_correlated_over_f(self, amp, index, f0):
        """Add 1/f noise that is perfectly correlated between frequencies.
        """

        if hasattr(self, "correlated_modes"):
            raise RuntimeError("The correlated modes already exist and would "
                               "be overwritten.")
        pass




