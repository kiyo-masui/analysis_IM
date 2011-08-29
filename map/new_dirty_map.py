"""Dirty map making module."""

import scipy as sp
import numpy.ma as ma

import core.algebra as al
import tools

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
        """Use this operator to convert a 'time' axis to a coordinate axis."""

        pass

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



class Noise_independant_channels(object):

    def __init__(self):
        pass

    def apply_inverse(self):
        pass

    def transform_map_space(self, Pointing, matrix_out=None):
        pass


class Noise(Noise_independant_channels):

    def __init__(self):
        pass


