"""Provides vector and matrix interfaces tailored to the needs of IM analysis.

At the heart of this module is the need to organize multidimensional data (such
as a map which has 3 axes: ra, dec, frequency) as a 1D vector for linear
algebra operations.  However we do not want to lose the multidimensional
organization of the data.  In addition there are efficiencies to be gained by
organizing data in a multidimensional way.  For example, some matrices will
be block diagonal in that they will not couple different frequency bins.  We
would like to exploit this.  To address this, classes are provided that store
data in numpy `ndarray` format but know how to reorganize the data into a
matrix or vector format.  Some basic matrix operations are also provided.

Also fully supported is writing the data to and from disk in a standard,
portable format that can easily be read in another language.  This is based on
numpy's NPY format.  Also fully supported is memory mapping; the ability to
manipulate an array stored on disk as if it was in memory.  This will be very
important as data sets become to large to store in memory and when operations
need to be performed in parallel using SCALAPACK.

Examples
--------

Matrices.

>>> import algebra as al
>>> mat_arr = np.arange(120)
>>> mat_arr.shape = (4, 5, 6)
>>> mat_arr = al.make_mat(mat_arr, row_axes=(0,), col_axes=(1, 2),
...                       axis_names=('mode', 'ra', 'dec'))
>>> mat_arr.row_shape()
(4,)
>>> mat_arr.col_shape()
(5, 6)
>>> mat_arr.mat_shape()
(4, 30)
>>> mat_arr.col_names()
('ra', 'dec')

Vectors.

>>> vect_arr = np.arange(30)
>>> vect_arr.shape = (5, 6)
>>> vect_arr = al.make_vect(vect_arr, axis_names=('ra', 'dec')
>>> vect_arr.mat_shape()
(30,)

Algebra.

>>> new = al.dot(mat_arr, vect_arr)
>>> new.shape
(4,)
>>> isinstance(new, al.vect_array)
True
>>> al.axes
('mode',)
"""

import os
import sys
import warnings
from functools import wraps
import operator

import scipy as sp
import numpy as np
import numpy.lib.format as npfor
from numpy.lib.utils import safe_eval

#import utils.cubic_conv_interpolation as cci
#from utils import cubic_conv_interpolation as cci

# TODO:
# when submitted as batch on Sunnyvale, the PYTHONPATH seems to get clobbered
# on the compute nodes. For now, just add these as-needed, but need better
# solution
try:
    import kiyopy.custom_exceptions as ce
except ImportError:
    import sys, site
    print "Your python path does not contain kiyopy (clobbered?)"
    site.addsitedir('/home/eswitzer/local/lib/')
    site.addsitedir('/home/eswitzer/local/lib/python2.6/site-packages/')
    import kiyopy.custom_exceptions as ce

# ---- Slight modification to the NPY format. ---------------------------------

def write_array_header_1_0(fp, d):
    """ Write the header for an array using the 1.0 format.

    This version of write array header has been modified to align the start of
    the array data with the 4096 bytes, corresponding to the page size of most
    systems.  This is so the npy files can be easily memmaped.

    Parameters
    ----------
    fp : filelike object
    d : dict
        This has the appropriate entries for writing its string representation
        to the header of the file.
    """
    import struct
    header = ["{"]
    for key, value in sorted(d.items()):
        # Need to use repr here, since we eval these when reading
        header.append("'%s': %s, " % (key, repr(value)))
    header.append("}")
    header = "".join(header)
    # Pad the header with spaces and a final newline such that the magic
    # string, the header-length short and the header are aligned on a 16-byte
    # boundary.  Hopefully, some system, possibly memory-mapping, can take
    # advantage of our premature optimization.
    # 1 for the newline
    current_header_len = npfor.MAGIC_LEN + 2 + len(header) + 1  
    topad = 4096 - (current_header_len % 4096)
    header = '%s%s\n' % (header, ' '*topad)
    if len(header) >= (256*256):
        raise ValueError("header does not fit inside %s bytes" % (256*256))
    header_len_str = struct.pack('<H', len(header))
    fp.write(header_len_str)
    fp.write(header)

def _replace_write_header(f) :
    """Wrap functions such that np.lib.format.write_array_header_1_0 is
    replaced by the local version, but only for the one function call."""
    @wraps(f)
    def wrapper(*args, **kwds):
        # Replace the header writer in the format module.
        tmp_write_header = npfor.write_array_header_1_0
        npfor.write_array_header_1_0 = write_array_header_1_0
        # Evaluate the function.
        try :
            result = f(*args, **kwds)
        finally :
            # Restore the header.
            npfor.write_array_header_1_0 = tmp_write_header
        return result
    return wrapper

# ---- Array classes for holding vectors and matricies. -----------------------

class info_array(sp.ndarray) :
    """A standard numpy ndarray object with a dictionary for holding extra info.

    This class should work exactly the same as a numpy ndarray object but has an
    attribute named info, which is a dictionary.  This class performs basic
    meta data handling for the higher lever classes that subclass this one:
    mat_array and vect_array. 

    Parameters
    ----------
    input_array : array like
        Array to converted to an info_array.  The info_array will be a
        view to the input array.
    info : dictionary
        Dictionary to be set as the `info` attribute (default is None, which
        implies create a new empty dictionary).

    Attributes
    ----------
    info : dictionary
        Holds any meta data associated with this array.  All items should be
        easily represented as a string so they can be written to and read from
        file.  Setting this attribute is generally not safe, but modifying it
        is.
    
    See Also
    --------
    info_memmap : Analogous class to this one with data stored on disk.
    vect_array : Vector object based on this class.
    mat_array : Matrix object based on this class.

    Notes
    -----
    All new from template array creation operations make a copy of the metadata
    not a reference.  To get a reference you need to explicitly make a call to
    the `view` function.

    See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html for more
    information
    """

    def __new__(cls, input_array, info=None):
        # Input array is an already formed ndarray instance.
        # We first cast to be our class type.
        obj = sp.asarray(input_array).view(cls)
        # Add the new attribute to the created instance.
        if not info is None :
            obj.info = info
        # Finally, we must return the newly created object.
        return obj

    def __array_finalize__(self, obj):
        # `info` is a reference to the origional only for an explicit call to
        # self.view() (view casting).  Otherwise we copy to protect the data.
        self.info = dict(getattr(obj, 'info', {}))

    def view(self, *args) :
        """Return a numpy view of self.
        
        This is mostly the same as the numpy version of this method, but it
        also makes the view's `info` attribute a reference to `self`s (where
        applicable).

        See Also
        --------
        np.ndarray.view
        """
        # Create the normal view.
        out = sp.ndarray.view(self, *args)
        # If it's info_array, replace the copy of the info attribute with a
        # reference (they will share metadata).
        if isinstance(out, info_array) :
            out.info = self.info
        return out

class info_memmap(sp.memmap) :
    """A standard numpy memmap object with a dictionary for holding extra info.

    This class should work exactly the same as a numpy memmap object but has an
    attribute named info, which is a dictionary. This class performs basic
    meta data handling for the higher lever classes that subclass this one:
    mat_memmap and vect_memmap.  This array is written to file at the same time
    that the memmap is flushed.
    
    Parameters
    ----------
    marray : numpy.memmap
        Array to be converted to an info_memmap.  The info_memmap will be a
        view to the input array.
    info : dictionary
        Dictionary to be set as the `info` attribute (default is None, which
        implies create a new empty dictionary).
    metafile : str 
        filename to write the metadata to.  In some
        versions of numpy, the metadata will be written to file even if the
        memmap is in read only mode.  To avoid this pass metafile=None, which
        prevents the metadata from being stored on disk at all.

    Attributes
    ----------
    info : dictionary
        Holds any meta data associated with this array.  All items should be
        easily represented as a string so they can be written to and read from
        file. Setting this attribute is generally not safe, but modifying it
        is.
    metafile : str
        filename where the metadata is written to.  `info` is written to this
        file whenever the `flush` method is called (which includes deletion of
        the object).  This can happen even if the memmap was opened in 'r'
        mode.  Set to None if you wish to protect the data on file.

    See Also
    --------
    info_array : Similar class with data stored in memory.
    vect_memmap : Vector object based on this class.
    mat_memmap : Matrix object based on this class.
    open_memmap : Open a file on disk as an info_memmap.

    Notes
    -----
    All new from template array creation operations make a copy of the metadata
    not a reference.  To get a reference you need to explicitly make a call to
    the `view` function.  Also, the metafile is set to None on new from
    template oberations.  
    
    See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html for more
    information
    """

    def __new__(cls, marray, info=None, metafile=None):
        # Input array is an already formed ndarray instance.
        # We first cast to be our class type.
        if not isinstance(marray, sp.memmap) :
            raise TypeError("info_memmaps can only be initialized off of "
                            "numpy memmaps.")
        obj = marray.view(cls)
        # Add the new attribute to the created instance.
        if info is None :
            info = {}
        obj.info = info
        obj.metafile = metafile
        # Finally, we must return the newly created object.
        return obj

    def __array_finalize__(self, obj):
        sp.memmap.__array_finalize__(self, obj)
        # Info is a reference to the origional for views.
        self.info = dict(getattr(obj, 'info', {}))
        # Do not copy the metafile attribute, new arrays will clobber the data.
        # This attribute is only copied on an explicit view() call.
        self.metafile = None

    def view(self, *args) :
        """Return a numpy view of self.
        
        This is mostly the same as the numpy version of this method, but it
        also makes the view's `info` attribute a reference to `self`s (where
        applicable).

        See Also
        --------
        np.ndarray.view
        """
        # Create the normal view.
        out = sp.memmap.view(self, *args)
        # If it's info_array, replace the copy of the info attribute with a
        # reference (they will share metadata).
        if isinstance(out, info_memmap) :
            out.info = self.info
            out.metafile = self.metafile
        return out

    def flush(self) :
        """Flush changes to disk.

        This method saves the info dictionary to metafile and then calls the
        flush method from the numpy memmap.
        """
        # Write the info dictionary to disk.
        self._info_flush()
        # Now flush the actual memmap.
        sp.memmap.flush(self)

    def _info_flush(self) :
        """Write the info array to disk only."""
        # Prior to numpy 1.5, we can't get the mode, so just assume we are
        # allowed to write
        mode = getattr(self, 'mode', 'r+')
        if ('+' in mode or 'w' in mode) and not self.metafile is None :
            # Convert info dictionary to a pretty string.
            infostring = repr(self.info)
            try:
                safe_eval(infostring)
            except SyntaxError :
                raise ce.DataError("Array info not representable as a string.")
            # Save the meta data.
            info_fid = open(self.metafile, 'w')
            try :
                info_fid.write(infostring)
            finally :
                info_fid.close()

    def __del__(self) :
        self._info_flush()
        sp.memmap.__del__(self)

    def __deepcopy__(self, copy) :
        """Not implemented, raises an exception."""
        raise NotImeplementedError("Deep copy won't work.")

def assert_info(array) :
    """Check if passed array is an info_array or info_memmap.

    Raises a ValueError if check fails.
    """
    if not (isinstance(array, info_array) or isinstance(array, info_memmap)) :
        raise TypeError("Array is not an algebra.info_array or "
                         "algebra.info_memmap.")


# ---- Functions for reading and writing the above to and from file. ----------

@_replace_write_header
def open_memmap(filename, mode='r+', dtype=None, shape=None,
                fortran_order=False, version=(1,0), metafile=None) :
    """Open a file and memory map it to an info_memmap object.

    This is similar to the numpy.lib.format.openmemmap() function but also
    deals with the meta data dictionary, which is read and written from a
    meta data file.

    The only extra argument over the numpy version is the meta data file name
    `metafile`.
    
    Parameters
    ----------
    metafile : str
        File name for which the `info` attribute of the returned info_memmap
        will be read from and written to. Default is None, where the it is
        assumed to be `filename` + ".meta".
        
    Returns
    -------
    marray : info_memmap
        The `info` is intialized as an empty dictionary if `mode` is 'w' or if
        the file corresponding to `metafile` does not exist.  The `metafile`
        attribute of marray is set to the `metafile` parameter unless `mode` is
        'r' or 'c' in which case it is set to None.
    """
    
    # Restrict to version (1,0) because we've only written write_header for
    # this version.
    if version != (1,0) :
        raise ValueError("Only version (1,0) is safe from this function.")
    # Memory map the data part.
    marray = npfor.open_memmap(filename, mode, dtype, shape, fortran_order,
                               version)
    # Get the file name for the meta data.
    if metafile is None :
        metafile = filename + '.meta'
    # Read the meta data if need be.
    if ('r' in mode or mode is 'c') and os.path.isfile(metafile) :
        info_fid = open(metafile, 'r')
        try :
            infostring = info_fid.readline() 
        finally : 
            info_fid.close()
        info = safe_eval(infostring)
    else :
        info = {}
    # In read mode don't pass a metafile to protect the meta data.
    if mode is 'r' or mode is 'c' :
        metafile = None
    marray = info_memmap(marray, info, metafile)
        
    return marray

def load(file, metafile=None) :
    """Open a .npy file and load it into memory as an info_aray.
    
    Similar to the numpy.load function.  Does not support memory
    mapping (use open_memmap).
    
    Parameters
    ----------
    file : file handle or str
        .npy file or file name to read the array from.
    metafile : str
        File name for which the `info` attribute of the returned info_array
        will be read from. Default is None, where the it is
        assumed to be the file name associated with `file` with ".meta"
        appended. If the file does not exist, the info attribute is initialized
        to an empty dictionary.

    Returns
    -------
    iarray : info_array object
    """
    
    # Load the data from .npy format.
    array = sp.load(file)
    # Figure out what the filename for the meta data should be.
    if metafile is None :
        try :
            fname = file.name
        except AttributeError :
            fname = file
        metafile = fname + ".meta"
    # Read the meta data.
    if os.path.isfile(metafile) :
        info_fid = open(metafile, 'r')
        try :
            infostring = info_fid.readline()
        finally : 
            info_fid.close()
        info = safe_eval(infostring)
    else :
        info = {}
    # Construct the infor array.
    array = info_array(array, info)
    
    return array

@_replace_write_header
def save(file, iarray, metafile=None, version=(1,0)) :
    """Save a info array to a .npy file and a metadata file.
    
    Similar to the numpy.save function.

    Parameters
    ----------
    file : file handle or str
        File or file name to write the array to in .npy format.
    iarray : info_array object or array with similar interface
        Array to be written to file with meta data.
    metafile : str
        File name for the meta data.  The `info` attribute of `iarray` will be
        written here. Default is None, where the it is
        assumed to be the file name associated with `file` with ".meta"
        appended.
    """
    
    # Restrict to version (1,0) because we've only written write_header for
    # this version.
    if version != (1,0) :
        raise ValueError("Only version (1,0) is safe from this function.")
    # Make sure that the meta data will be representable as a string.
    infostring = repr(iarray.info)
    try:
        safe_eval(infostring)
    except SyntaxError :
        raise ce.DataError("Array info not representable as a string.")
    # Save the array in .npy format.
    if isinstance(file, basestring):
        fid = open(file, "wb")
    else:
        fid = file
    npfor.write_array(fid, iarray, version=version)
    # Figure out what the filename for the meta data should be.
    if metafile is None :
        try :
            fname = file.name
        except AttributeError :
            fname = file
        metafile = fname + ".meta"
    # Save the meta data.
    info_fid = open(metafile, 'w')
    try :
        info_fid.write(infostring)
    finally :
        info_fid.close()


# ---- Functions for manipulating above arrays as matrices and vectors. -------

#### Some helper functions that perform some quick checks. ####

def _set_type_axes(array, type, axis_names) :
    """Sets the array.info['type'] and array.info[axes] metadata and does some
    checks.  Used in vect and mat constructors.
    """
    
    assert_info(array)

    if axis_names is None :
        axes = (None,)*array.ndim
    else :
        _check_axis_names(array, axis_names)
        axes = tuple(axis_names)
    array.info['type'] = type
    array.axes = axes

def _check_axis_names(array, axis_names=None) :
    """Checks that axis names  sequence is valid for array."""

    if axis_names is None :
        axis_names = array.axes

    if len(axis_names) != array.ndim :
        raise ValueError("axis_names parameter must be a sequence of length "
                         "arr.ndim")
    else :
        for name in axis_names :
            if (not isinstance(name, str)) and (not name is None) :
                raise TypeError("Invalid axis name.")

def _check_rows_cols(arr, row_axes=None, col_axes=None) :
    """Check that rows and cols are valid for the matrix."""

    if row_axes is None and col_axes is None :
        row_axes = arr.info['rows']
        col_axes = arr.info['cols']
    # Both parameters had better be sequences of integers.
    for ind in row_axes :
        if not ind in range(arr.ndim) :
            raise ValueError("Invalid row axes.")
    for ind in col_axes :
        if not ind in range(arr.ndim) :
            raise ValueError("Invalid col axes")
    # Make sure each axis is spoken for.
    for ii in range(arr.ndim) :
        if (not ii in row_axes) and (not ii in col_axes) :
            raise ValueError("Every axis must be identified varying over "
                             "as the matrix row, column or both.")


#### Common class definitions ####

class alg_object(object) :
    """Base class for all vectors and matricies.
    
    This is not an actual class by itself, just defines some methods common to
    both `mat` objects and `vect` objects.
    """

    def __array_finalize__(self, obj) :
        self.info_base.__array_finalize__(self, obj)
        if (not obj is None) and (self.shape != obj.shape) :
            self.__class__ = self.info_base
    
    def set_axis_info(self, axis_name, centre, delta) :
        """Set meta data for calculating values of an axis.
        
        This provides the meta data required to calculate one of the axes.
        This data is stored in the `info` attribute, which is carried
        around by this class and easily written to disk.

        The information provided is subsequently used in `get_axis` to
        calculate the values along a given axis. 

        Parameters
        ----------
        axis_name : str
            Name of the axis for which you are setting the meta data.
            Must match one of the entries of the axes attribute.
        centre : float
            The value of the axis at the centre bin (indexed by n//2),
            where n = self.shape[i] and self.axes[i] is `axis_name`.
        delta : float
            The width of each bin.

        See Also
        --------
        copy_axis_info
        get_axis

        Examples
        --------
        >>> import algebra
        >>> a = algebra.make_vect(sp.zeros(5, 5), axis_names=('ra', 'dec'))
        >>> a.set_axis_info('ra', 2, 0.5)
        >>> a.get_axis('ra')
        array([1.0, 1.5, 2.0, 2.5, 3.0])
        """

        if not axis_name in self.axes :
            raise ValueError("axis_name not in self.axes.")

        self.info[axis_name + '_centre'] = float(centre)
        self.info[axis_name + '_delta'] = float(delta)

    def copy_axis_info(self, alg_obj) :
        """Set the axis info by copying from another alg_object.

        This transfers meta data that is set with `set_axis_info` another
        alg_object instance.
        
        Parameters
        ----------
        alg_obj : algebra.alg_object instance
            Object from which to copy axis meta data.  Meta data for all axis
            names that occur in both `alg_obj.axes` and `self.axes` is copied.

        See Also
        --------
        set_axis_info
        get_axis
        """
        
        if isinstance(alg_obj, dict):
            info = alg_obj
        else:
            info = alg_obj.info
        for axis in self.axes :
            if axis in info["axes"] :
                try :
                    centre = info[axis + '_centre']
                    delta = info[axis + '_delta']
                except KeyError :
                    continue
                self.info[axis + '_centre'] = centre
                self.info[axis + '_delta'] = delta

    def get_axis(self, axis_name) :
        """Calculate the array representing a named axis.  
        
        For a given axis name, calculate the 1D array that gives the value of
        that axis.  This requires that the relevant meta data be set by
        `set_axis_info`.

        Parameters
        ----------
        axis_name : str or int
            Name of the axis to be calculated.  `axis_name` must occur in the
            `axes` attribute.  If an int is passed, than it is convered to a
            string by indexing the `axes` attribute.

        Returns
        -------
        axis_array : np.ndarray
            The array corresponding to the values of the axis quantity along
            it's axis.

        See Also
        --------
        set_axis_info
        copy_axis_info
        """
        
        if isinstance(axis_name, int) :
            axis_name = self.axes[axis_name]
        len = self.shape[self.axes.index(axis_name)]
        return (self.info[axis_name + '_delta']*(sp.arange(len) - len//2) 
                + self.info[axis_name + '_centre'])

    def slice_interpolate_weights(self, axes, coord, kind='linear') :
        """Get the interpolation weights for a subset of the dimensions.

        This method gets the interpolation weights for interpolating the
        alg_object is some subset of it's dimensions.  This provides the
        freedom in the uninterpolated dimensions to either slice or otherwise
        index the array.

		Note: When using the cubic convolution interpolation, the points
		may out of bounds. That is how the algorithm works when interpolating
		right beside a boundary. If the value of those nodes are needed, use
		'get_value' from the cci helper module.
		Also, some weights will be negative as needed by the algorithm.

        Parameters
        ----------
        axes : int or sequence of ints (length N)
            Over which axes to interpolate.
        coord : float or sequence of floats (length N)
            The coordinate location to interpolate at.
        kind : string
            The interpolation algorithm.  Options are: 'linear' or 'nearest'.
			And now 'cubic' too!

        Returns
        -------
        points : array of ints shape (M x N)
            The indices for the N `axes` at the M interpolation data points
            that are used.
        weights : array of floats length M
            Weights for the interpolations data points.
        """

        if not hasattr(axes, '__iter__') :
            axes = (axes,)
        if not hasattr(coord, '__iter__') :
            coord = (coord,)
        n = len(axes)
        if n != len(coord) :
            message = "axes and coord parameters must be same length."
            raise ValueError(message)
        if kind in ('linear',) :
            # Any interpolation scheme that only depends on data points
            # directly surrounding. There are 2^n of them.
            m = 2**n
            points = sp.empty((m, n), dtype=int)
            weights = sp.empty((m,), dtype=float)
            # Find the indices of the surrounding points, as well as the 
            single_inds = sp.empty((n, 2), dtype=int)
            normalized_distance = sp.empty((n, 2), dtype=float)
            for ii in range(n) :
                axis_ind = axes[ii]
                value = coord[ii]
                # The spacing between points of the axis we are considering.
                delta = abs(self.info[self.axes[axis_ind] + "_delta"])
                # For each axis, find the indicies that surround the
                # interpolation location.
                axis = self.get_axis(axis_ind)
                if value > max(axis) or value < min(axis) :
                    message = ("Interpolation coordinate outside of "
                               "interpolation range.  axis: " + str(axis_ind)
                               + ", coord: " + str(value) + ", range: "
                               + str((min(axis), max(axis))))
                    raise ce.DataError(message)
                distances = abs(axis - value)
                min_ind = distances.argmin()
                single_inds[ii, 0] = min_ind
                normalized_distance[ii, 0] = distances[min_ind]/delta
                distances[min_ind] = distances.max() + 1
                min_ind = distances.argmin()
                single_inds[ii, 1] = min_ind
                normalized_distance[ii, 1] = distances[min_ind]/delta
            # Now that we have all the distances, figure out all the weights.
            for ii in range(m) :
                temp_ii = ii
                weight = 1.0
                for jj in range(n) :
                    points[ii, jj] = single_inds[jj, temp_ii%2]
                    weight *= (1.0 - normalized_distance[jj, temp_ii%2])
                    temp_ii = temp_ii//2
                weights[ii] = weight
        elif kind == 'nearest':
            # Only one grid point to consider and each axis is independant.
            m = 1
            weights = sp.ones(1, dtype=int)
            points = sp.empty((1, n), dtype=int)
            # Loop over the axes we are interpolating.
            for ii in range(n):
                axis_name = self.axes[axes[ii]]
                axis_centre = self.info[axis_name + "_centre"]
                axis_delta = self.info[axis_name + "_delta"]
                index = (coord[ii] - axis_centre)/axis_delta
                index += self.shape[axes[ii]]//2
                if index < 0 or index > self.shape[axes[ii]] - 1:
                    message = ("Interpolation coordinate outside of "
                               "interpolation range.  axis: " + str(axes[ii])
                               + ", coord: " + str(coord[ii]) + ".")
                    raise ce.DataError(message)
                points[0, ii] = round(index)
        elif kind == 'cubic':
            import utils.cubic_conv_interpolation as cci
            # Make sure the given point is an array.
            pnt = sp.array(coord)
            # Get the array containing the first value in each dimension.
            # And get the arrays of deltas for each dimension.
            x0 = []
            step_sizes = []
            for ax in axes:
                ax_name = self.axes[ax]
                x0.append(self.get_axis(ax_name)[0])
                ax_delta_name = ax_name + "_delta"
                step_sizes.append(self.info[ax_delta_name])
            x0 = sp.array(x0)
            step_sizes = sp.array(step_sizes)
            # Get the maximum possible index in each dimension in axes.
            max_inds = np.array(self.shape) - 1
            max_needed_inds = []
            for ax in axes:
                max_needed_inds.append(max_inds[ax])
            max_inds = np.array(max_needed_inds)
            # If there are less than four elements along some dimension
            # then raise an error since cubic conv won't work.
            too_small_axes = []
            for i in range(len(max_inds)):
                if max_inds[i] < 3:
                    too_small_axes.append(axes[i])
            if not (len(too_small_axes) == 0):
                msg = "Need at least 4 points for cubic interpolation " + \
                      "on axis (axes): " + str(too_small_axes)
                raise ce.DataError(msg)
            # Get the nodes needed and their weights.
            points, weights = cci.interpolate_weights(axes, pnt, x0, \
                                  step_sizes, max_inds)
        else :
            message = "Unsupported interpolation algorithm: " + kind
            raise ValueError(message)
        return points, weights
        
    def slice_interpolate(self, axes, coord, kind='linear') :
        """Interpolate along a subset of dimensions.

        This method interpolate the array object along some subset of it's
        dimensions.  The array is sliced long the uninterpolated dimensions.

        Parameters
        ----------
        axes : int or sequency of ints (length N)
            Over which axes to interpolate.
        coord : float or sequence of floats (length N)
            The coordinate location to interpolate at.
        kind : string
            The interpolation algorithm.  Options are: 'linear', 'nearest'.
        
        Returns
        -------
        slice : numpy array
            The array has a shape corresponding to the uninterpolated
            dimensions of `self`.  That is it has dimensions the same as `self`
            along dimension included in `axes` which are interpolated.
        """
        
        if not hasattr(axes, '__iter__') :
            axes = (axes,)
        if not hasattr(coord, '__iter__') :
            coord = (coord,)
        n = len(axes)
        if n != len(coord) :
            message = "axes and coord parameters must be same length."
            raise ValueError(message)
        # Get the contributing points and their weights.
        points, weights = self.slice_interpolate_weights(axes, coord, kind)
        # Sum up the points
        q = points.shape[0]
        out = 0.0
        for ii in range(q) :
            index = [slice(None)] * self.ndim
            for jj, axis in enumerate(axes) :
                index[axis] = points[ii, jj]
            index = tuple(index)
            out += weights[ii] * self[index]
        return out


#### Vector class definitions ####

class vect(alg_object) :
    """Multidimentional array interpreted as a vector.
    
    This class gets most of its functionality from the numpy ndarray class.
    In addition it provides support for orgainizing it's data as a vector.
    This class comes in two flavours: `vect_array` and `vect_memmap`
    depending on whether the array is stored in memory or on disk.  The raw
    `vect` class is not a valid class by itself.

    The vector representation of the array is the flattend array.

    One of the features that `vect`s implement is named axes.  This allows
    maps to carry axis information with them, among other things.  For
    `vect`s axis names must be unique (This is not true of `mat`s).

    Parameters
    ----------
    input_array : info_array (for vect_array) or info_memmap (for vect_memmap)
        Array to be converted to a vect.
    axis_names : tuple of strings, optional
        The sequence contains the name of each axis.  This sequence will be
        stored in the `axes` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the axis names.

    Attributes
    ----------
    axes : tuple of strings
        The names of each of the axes of the array.

    See Also
    --------
    make_vect : Cast any array as a vector.
    info_array, info_memap : Base classes that handle meta data.

    Notes
    -----
    Since much of the functionality provided by this class is only valid
    for a certain shape of the array, shape changing operations in
    general return an `info_array` or `info_memmap` as appropriate (except
    explicit assignment to vect.shape).

    The `axes` attribute is actually stored in the `info_array`'s info
    dictionary.  This is just an implimentation detail.
    """

    __array_priority__ = 2.0
    
    # self.info_base is set in the class factory.
    def __new__(cls, input_array, axis_names=None) :
        
        if not isinstance(input_array, cls.info_base) :
            raise ValueError("Array to convert must be instance of " +
                             str(cls.info_base))
        obj = input_array.view(cls)
        if obj.info.has_key('type') :
            if not axis_names is None :
                warnings.warn("Initialization argument ignored. Requisite "
                              "metadata for vector already exists. "
                              "Clear info dictionary if you want opposite "
                              "behaviour.")
            if not obj.info['type'] == 'vect' :
                raise ValueError("Meta data present is incompatible.")
            _check_axis_names(obj)
        else :
            _set_type_axes(obj, 'vect', axis_names)
        return obj

    def __setattr__(self, name, value) :
        if name == 'axes' :
            _check_axis_names(self, value)
            self.info['axes'] = value
        else :
            self.info_base.__setattr__(self, name, value)

    def __getattr__(self, name) :
        if name == 'axes' :
            return self.info['axes']
        else :
            # Since numpy uses __get_attribute__ not __getattr__, we should
            # raise an Attribute error.
            raise AttributeError("Attribute " + name + " not found.")

    def flat_view(self) :
        """Returns a view of the vector that has been flattened.

        The view will be cast as a scipy.ndarray and its shape will be
        (self.size, ).  It is a view, so writing to it will write back to the
        original vector object.

        Returns
        -------
        flat_view : np.ndarray
            A view of `self` as an ndarray, flattened to 1D.
        """
        
        flat = self.view(sp.ndarray)
        flat.shape = (self.size, )
        return flat

    def mat_shape(self) :
        """Get the shape of the represented matrix (vector)."""
        _check_axis_names(self)
        return (self.size,)


def _vect_class_factory(base_class) :
    """Internal class factory for making a vector class that inherits from
    either info_array or info_memmap."""

    if (not base_class is info_array) and (not base_class is info_memmap) :
        raise TypeError("Vectors inherit from info arrays or info memmaps.")

    class vect_class(vect, base_class) :
        __doc__ = vect.__doc__
        info_base = base_class

    return vect_class

vect_array = _vect_class_factory(info_array)
vect_array.__name__ = 'vect_array'
vect_memmap = _vect_class_factory(info_memmap)
vect_memmap.__name__ = 'vect_memmap'

def make_vect(array, axis_names=None) :
    """Do whatever it takes to make a vect out of an array.
    
    Convert any class that can be converted to a vect (array, info_array,
    memmap, info_memmap) to the appropriate vect object (vect_array,
    vect_memmap).

    This convenience function just simplifies the constructor hierarchy.  
    Normally to get a vect out of an array, you would need to construct an
    intermediate info_array object.  This bypasses that step.
    
    Parameters
    ----------
    array : array_like
        Array to be converted to vect object (if possible).
    axis_names : tuple of strings, optional
        The sequence contains the name of each axis.  This sequence will be
        stored in the `axes` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the axis names.

    Returns
    -------
    vect_arr : vect_array or vect_memmap
        A view of `array` converted to a vect object.
    
    """

    if isinstance(array, sp.memmap) :
        if not isinstance(array, info_memmap) :
            array = info_memmap(array)
        return vect_memmap(array, axis_names)
    elif isinstance(array, sp.ndarray) :
        if not isinstance(array, info_array) :
            array = info_array(array)
        return vect_array(array, axis_names)
    else :
        raise TypeError("Object cannot be converted to a vector.")

#### Matrix class definitions ####

class mat(alg_object) :
    """Multidimentional array interpreted as a matrix.
    
    This class gets most of its functionality from the numpy ndarray class.
    In addition it provides support for organizing it's data as a vector.
    This class comes in two flavours: `mat_array` and `mat_memmap`
    depending on whether the array is stored in memory or on disk.  The raw
    `mat` class is not a valid class by itself.

    To make the association between a multidimentional array and a matrix,
    each axis of the array must be identified as varying over either the
    rows or columns of a matrix.  For instance the shape of the array could
    be (3, 5, 7).  We could identify the first axis as a row axis and the
    second two as column axes in which case the matrix would have 3 rows
    and 35 colums.  We generally make the rows axes left of the columns and
    many algorithms assume this.  It is also possible for an axis to be
    identified as both a row and a column, in which case the matrix is
    block diagonal over that axis.  Generally, the block diagonal axes are
    the leftmost.

    Like `vect`s, mats have named axes, however 2 axes may have the same
    name as long as one is identified as a row axis and the other as a col
    axis.

    Parameters
    ----------
    input_array : info_array (for mat_array) or info_memmap (for
                  mat_memmap)
        Array to be converted to a vect.
    row_axes : tuple of ints
        Sequence contains the axis numbers of the array to identify as
        varying over the matrix rows. This sequence is stored in the
        `rows` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the rows.
    col_axis : tuple of ints
        Sequence contains the axis numbers of the array to identify as
        varying over the matrix columns. This sequence is stored in the
        `cols` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the cols.
    axis_names : tuple of strings, optional
        The sequence contains the name of each axis.  This sequence will be
        stored in the `axes` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the axis names.

    Attributes
    ----------
    axes : tuple of strings
        The names of each of the axes of the array.
    rows : tuple of ints
        Which of the array's axes to identify as varying over the matrix
        rows.
    cols : tuple of ints
        Which of the array's axes to identify as varying over the matrix
        columns.

    Notes
    -----
    Since much of the functionality provided by this class is only valid
    for a certain shape of the array, shape changing operations in
    general return an `info_array` or `info_memmap` as appropriate (except
    explicit assignment to mat.shape).

    The `axes`, `rows` and `cols` attributes are actually stored in the 
    `info_array`'s info dictionary.  This is just an implementation detail.

    See Also
    --------
    vect_array, vect_memmap : Vector classes.
    make_mat : Function that casts any array as a matrix.
    info_array, info_memmap : Base classes that handle meta data.
    """

    __array_priority__ = 3.0
    
    def __new__(cls, input_array, row_axes=None, col_axes=None, 
                axis_names=None) :
        
        if not isinstance(input_array, cls.info_base) :
            raise ValueError("Array to convert must be instance of " +
                             str(cls.info_base))
        
        obj = input_array.view(cls)

        if obj.info.has_key('type') :
            if ((not axis_names is None) or (not row_axes is None) or 
                (not col_axes is None)) :
                warnings.warn("Initialization argument ignored. Requisite "
                              "metadata for matrix already exists. "
                              "Clear info dictionary if you want opposite "
                              "behaviour.")
            if not obj.info['type'] == 'mat' :
                raise ValueError("Meta data present is incompatible.")
            _check_axis_names(obj)
            _check_rows_cols(obj)
        else :
            if ((row_axes is None) and (col_axes is None)
                and (input_array.ndim==2)) :
                row_axes = (0,)
                col_axes = (1,)
            else :
                _check_rows_cols(input_array, row_axes, col_axes)
            _set_type_axes(obj, 'mat', axis_names)
            obj.rows = row_axes
            obj.cols = col_axes
        return obj

    def __setattr__(self, name, value) :
        if name == 'axes' :
            _check_axis_names(self, value)
            self.info['axes'] = value
        elif name == 'rows' :
            for ind in value :
                if not ind in range(self.ndim) :
                    raise ValueError("Invalid row axes.")
            self.info['rows'] = tuple(value)
        elif name == 'cols' :
            for ind in value :
                if not ind in range(self.ndim) :
                    raise ValueError("Invalid col axes.")
            self.info['cols'] = tuple(value)
        else :
            self.info_base.__setattr__(self, name, value)

    def __getattr__(self, name) :
        if name == 'axes' :
            return self.info['axes']
        elif name == 'rows' :
            return self.info['rows']
        elif name == 'cols' :
            return self.info['cols']
        else :
            # Since numpy uses __get_attribute__ not __getattr__, we should
            # raise an Attribute error.
            raise AttributeError("Attribute " + name + " not found.")

    def check_rows_cols(self) :
        """Check that rows and cols are valid for the matrix.
        
        Raises an exception if the rows or columns are invalid.
        """
        _check_rows_cols(self)

    def assert_axes_ordered(self) :
        """Enforces a specific ordering to the matrix row and column axis
        associations.
        """

        rows = self.info['rows']
        cols = self.info['cols']
        r_ind = len(rows) - 1
        c_ind = len(cols) - 1
        in_rows = False
        in_cols = True

        for axis in range(self.ndim-1, -1, -1) :
            if cols[c_ind] == axis and in_cols and c_ind >= 0:
                c_ind -= 1
                continue
            elif in_cols :
                in_cols = False
                in_rows = True
            if rows[r_ind] == axis and in_rows and r_ind >= 0:
                r_ind -= 1
                continue
            elif in_rows :
                in_rows = False
            if rows[r_ind] == axis and cols[c_ind] == axis :
                r_ind -= 1
                c_ind -= 1
            else :
                raise NotImplementedError("Matrix row and column array axis"
                                          "associations not ordered correctly.")
    
    def get_num_blocks(self, return_block_shape=False, 
                       return_n_axes_diag=False) :
        """Get the number of blocks in a block diagonal matrix."""
        
        shape = self.mat_shape()
        # Current algorithm assumes specific format.
        self.assert_axes_ordered()
        
        diag_axes = [ii for ii in range(self.ndim) if ii in self.rows 
                     and ii in self.cols]
        num_blocks = sp.prod([self.shape[ii] for ii in diag_axes])
        if return_block_shape and return_n_axes_diag :
            return num_blocks, (shape[0]/num_blocks, shape[1]/num_blocks), \
                   len(diag_axes)
        elif return_block_shape :
            return num_blocks, (shape[0]/num_blocks, shape[1]/num_blocks)
        elif return_n_axes_diag :
            return num_blocks, len(diag_axes)
        else :
            return num_blocks
    
    def row_names(self) :
        """Return the axis names that correspond to rows."""

        names = ()
        for axis_ind in self.rows :
            names = names + (self.axes[axis_ind],)

        return names

    def col_names(self) : 
        """Return the axis names that correspond to columns."""

        names = ()
        for axis_ind in self.cols :
            names = names + (self.axes[axis_ind],)

        return names

    def row_shape(self) :
        """Return the shape of the array only including axes that correspond to
        rows."""

        shape = ()
        for axis_ind in self.rows :
            shape = shape + (self.shape[axis_ind],)

        return shape
    
    def col_shape(self) :
        """Return the shape of the array only including axes that correspond to
        rows."""

        shape = ()
        for axis_ind in self.cols :
            shape = shape + (self.shape[axis_ind],)

        return shape

    def mat_shape(self) :
        """Get the shape of the represented matrix."""

        self.check_rows_cols()
        _check_axis_names(self)
        nrows = 1
        for axis in self.rows :
            nrows *= self.shape[axis]
        ncols = 1
        for axis in self.cols :
            ncols *= self.shape[axis]
        return (nrows, ncols)

    def iter_blocks(self) :
        """Returns an iterator over the blocks of a matrix."""
        
        # Build the iterator class.
        class iterator(object) :
            
            def __init__(self, arr) :
                self.arr = arr
                self.n_blocks, self.block_shape, self.n_axes_diag = \
                    arr.get_num_blocks(True, True)
                self.ii = 0

            def __iter__(self) :
                return self

            def next(self) :
                if self.ii >= self.n_blocks :
                    raise StopIteration()
                else :
                    # Find the indices for this block.
                    array_index = ()
                    tmp_block_num = self.ii
                    self.ii += 1
                    for jj in range(self.n_axes_diag -1, -1, -1) :
                        array_index = ((tmp_block_num%self.arr.shape[jj],) 
                                       + array_index)
                        tmp_block_num = tmp_block_num//self.arr.shape[jj]
                    # Return the data.
                    return sp.reshape(self.arr[array_index], self.block_shape)
        
        # Initialize it and return it.
        return iterator(self)
    
    def _iter_row_col_index(self, row_or_col) :
        """Implementation of iter_col_index and iter_row_index."""
        #
        # Build the iterator class.
        class iterator(object) :

            def __init__(self, arr) :
                self.ndim = arr.ndim
                self.ii = 0
                # Axes that we will iteration over.
                if row_or_col == 'col' :
                    self.axes = list(arr.cols)
                    other_axes = arr.rows
                elif row_or_col == 'row' :
                    self.axes = list(arr.rows)
                    other_axes = arr.cols
                else :
                    raise RunTimeError()
                # Do not iterate over axes that are shared with rows.
                for oth in other_axes :
                    if oth in self.axes :
                        self.axes.remove(oth)
                # Get the shape of these axes, as well as the total size.
                self.shape = ()
                self.size = 1
                for axis in self.axes :
                    self.shape += (arr.shape[axis],)
                    self.size *= arr.shape[axis]

            def __iter__(self) :
                return self

            def next(self) :
                inds = ()
                ii = self.ii
                self.ii += 1
                if ii >= self.size :
                    raise StopIteration()
                # The sequence that will eventually be used to subscript the
                # array.
                array_index = [slice(sys.maxint)] * self.ndim
                # Get the indices.  Loop through the axes backward.
                for jj in range(len(self.axes) - 1, -1, -1) :
                    array_index[self.axes[jj]] = ii%self.shape[jj]
                    ii = ii//self.shape[jj]
                # Return the indices.
                return array_index

        # Initiallize and return iterator.
        return iterator(self)

    def iter_row_index(self) :
        """Returns an iterator over row axes of the mat.

        This iterates over all the axes that are assotiated only with rows
        of the mat.  Any axis that is identified as both a column and a row is
        not iterated over.  The iterator returns an tuple that can subscript
        the mat (not a view of the mat).  This is useful when you have an
        operation that has to be applied uniformly to all the columns of a mat.

        Examples
        --------
        >>> for index in mat.iter_row_index() :
        >>>     sub_arr = mat[index]
        >>>     sub_arr.shape == mat.col_shape()
        True
        """

        return self._iter_row_col_index('row')

    def iter_col_index(self) :
        """Returns an iterator over column axes of the mat.

        This iterates over all the axes that are assotiated only with columns
        of the mat.  Any axis that is identified as both a column and a row is
        not iterated over.  The iterator returns an tuple that can subscript
        the mat (not a view of the mat).  This is useful when you have an
        operation that has to be applied uniformly to all the columns of a mat.

        Examples
        --------
        >>> for index in mat.iter_col_index() :
        >>>     sub_arr = mat[index]
        >>>     sub_arr.shape == mat.row_shape()
        True
        """

        return self._iter_row_col_index('col')

    def mat_diag(self) :
        """Get the daigonal elements of the matrix, as a vect object."""

        # Current algorithm assumes specific format.
        self.assert_axes_ordered()
        # We expect a square matrix
        shape = self.mat_shape()
        if shape[0] != shape[1] :
            raise NotImplementedError("Only works for square mats.")
        # output memory
        out = sp.empty((shape[0],))
        # Figure out how many axes are in both row and col (and therefore block
        # diagonal).
        n_blocks, block_shape = self.get_num_blocks(True, False)
        # For square matricies, n_blocks*block_shape[ii] == shape[ii].
        block_size = shape[0]//n_blocks
        # Transfer over the diagonals.
        for ii, mat_block in enumerate(self.iter_blocks()) :
            out[ii*block_size:(ii+1)*block_size] = sp.diag(mat_block)

        # Now make this a vect object and transfer the relevant metadata.
        if self.row_shape() == self.col_shape() :
            out.shape = self.row_shape()
        out = make_vect(out)
        if self.row_names() == self.col_names() :
            out.axes = self.row_names()
            out.copy_axis_info(self)
        
        return out

    def expand(self) :
        """Calculates expanded matrix in 2 dimensional form.

        Takes an arbitrary matrix and returns the expanded version of it,
        as matrix with internal array dimensions of shape(mat).  If the
        original matrix has efficiency from any block diagonal structure, this
        is lost in the returned matrix.
        """

        # XXX Obviouse improvement: Check if this matrix is already full 
        # (ie not diagonal structure)
        # and if so, return a view.
        
        # Also verifies the validity of the matrix.
        shape = self.mat_shape()
        # Current algorithm assumes specific format.
        self.assert_axes_ordered()
        # Allocate memory.
        out_mat = sp.zeros(shape, dtype=self.dtype)
        out_mat = info_array(out_mat)
        out_mat = mat_array(out_mat)

        # Figure out how many axes are in both row and col (and therefore block
        # diagonal).
        n_blocks, block_shape = self.get_num_blocks(True, False)
        
        # Loop over the blocks and assign data.
        for ii, mat_block in enumerate(self.iter_blocks()) :
            # Figure out where this block starts.
            row_start = ii*block_shape[0]
            col_start = ii*block_shape[1]
            out_mat[row_start:row_start + block_shape[0], 
                    col_start:col_start + block_shape[1]] = mat_block
        return out_mat

    def mat_transpose(self):
        """Transpose the matrix.

        Returns an `mat` object with rows and columns exchanged.  The
        underlying array is not modified.

        Returns
        -------
        transposed_matrix : mat object
            A view of `self`, with different meta data such that the matrix is
            transposed.
        """
        
        # Make a copy of the info dictionary.
        info = dict(self.info)
        # Make a view of self.
        out = self.view()
        # Replace the info dictionary (which is a reference to self.info) with
        # the copy.
        out.info = info
        # Transpose the axes.
        out.cols = self.rows
        out.rows = self.cols
        return out

    def index_axis(self, axis, index):
        """Remove an axis by indexing with an integer.
        
        Returns a view of the `mat` with an axis removed by indexing it with
        the passed integer.  The output retains it's `mat` characteristics,
        which are updated to reflect the removal of the index.

        Parameters
        ----------
        axis : integer
            Axis to remove.
        index : integer
            Which entry to choose from `axis`.

        Returns
        -------
        out : mat or vect
            If the removal of the axis deplete either the rows or columns, a
            vect object is returned.  Otherwise a mat is returned.
        """
        
        # Suport negitive indicies.
        if axis < 0:
            axis = self.ndim + axis
        # Get copies of rows and colums ommitting the removed axis.  Also, any
        # indecies that are greater than axis, need to decremented.
        rows = list(self.rows)
        if axis in rows:
            rows.remove(axis)
        cols = list(self.cols)
        if axis in cols:
            cols.remove(axis)
        for ii in range(len(rows)):
            if rows[ii] > axis:
                rows[ii] -= 1
        for ii in range(len(cols)):
            if cols[ii] > axis:
                cols[ii] -= 1
        rows = tuple(rows)
        cols = tuple(cols)
        # New name attribute.
        names = self.axes[:axis] + self.axes[axis+1:]
        # Index the array.
        slice_index = [slice(None)] * self.ndim
        slice_index[axis] = index
        out = sp.asarray(self[slice_index])
        # Set the matrix meta data.
        if not rows or not cols:
            out = make_vect(out, axis_names=names)
        else:
            out = make_mat(out, axis_names=names, row_axes=rows, col_axes=cols)
        out.copy_axis_info(self)
        return out


def _mat_class_factory(base_class) :
    """Internal class factory for making a matrix class that inherits from
    either info_array or info_memmap."""

    if (not base_class is info_array) and (not base_class is info_memmap) :
        raise TypeError("Matrices inherit from info arrays or info memmaps.")

    class mat_class(mat, base_class) :
        __doc__ = mat.__doc__
        info_base = base_class
    return mat_class

mat_array = _mat_class_factory(info_array)
mat_array.__name__ = 'mat_array'
mat_memmap = _mat_class_factory(info_memmap)
mat_memmap.__name__ = 'mat_memmap'

def make_mat(array, row_axes=None, col_axes=None, axis_names=None) :
    """Do what ever it takes to make a mat out of an array.
    
    Convert any class that can be converted to a mat (array, info_array,
    memmap, info_memmap) to the appropriate mat object (mat_array or
    mat_memmap).
    
    This convieiance function just simplifies the constructor heirarchy.  
    Normally to get an mat out of an array, you would need to construct an
    intermediate info_array object.  This bypasses that step.
    
    Parameters
    ----------
    array : array_like
        Array to be converted to mat object (if possible).
    row_axes : tuple of ints
        Sequence contains the axis numbers of the array to identify as
        varying over the matrix rows. This sequence is stored in the
        `rows` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the rows.
    col_axis : tuple of ints
        Sequence contains the axis numbers of the array to identify as
        varying over the matrix columns. This sequence is stored in the
        `cols` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the cols.
    axis_names : tuple of strings, optional
        The sequence contains the name of each axis.  This sequence will be
        stored in the `axes` attribute.  This parameter is ignored if
        `input_array`'s info attribute already contains the axis names.

    Returns
    -------
    mat_arr : mat_array or mat_memmap
        A view of `array` converted to a mat object.
    """

    if isinstance(array, sp.memmap) :
        if not isinstance(array, info_memmap) :
            array = info_memmap(array)
        return mat_memmap(array, row_axes, col_axes, axis_names)
    elif isinstance(array, sp.ndarray) :
        if not isinstance(array, info_array) :
            array = info_array(array)
        return mat_array(array, row_axes, col_axes, axis_names)
    else :
        raise TypeError("Object cannot be converted to a matrix.")

#### Functions ####

def dot(arr1, arr2, check_inner_axes=True) :
    """Perform matrix multiplication."""

    shape1 = arr1.mat_shape()
    shape2 = arr2.mat_shape()

    if shape1[-1] != shape2[0] :
        raise ValueError("Matrix dimensions incompatible for matrix "
                         "multiplication.")
    # Matrix-vector product case.
    if len(shape1) == 2 and len(shape2) ==1 :
        # Strict axis checking has been requested, make sure that the axis
        # number, lengths and names of the input vector are equal to the 
        # column axis names of the input matrix.
        if check_inner_axes : 
            if not arr2.ndim == len(arr1.cols) :
                raise ce.DataError("Matrix column axis number are not the "
                                   "same as vector ndim and strict checking "
                                   "has been requested.")
            for ii, name in enumerate(arr2.axes) :
                if not arr1.shape[arr1.cols[ii]] == arr2.shape[ii] :
                    raise ce.DataError("Matrix column axis lens are not the "
                                       "same as vector axis lens and strict "
                                       "checking has been requested.")
                if not name == arr1.axes[arr1.cols[ii]] :
                    raise ce.DataError("Matrix column axis names are not the "
                                       "same as vector axes names and strict "
                                       "checking has been requested.")
        # Figure out what the output vector is going to look like.
        out_shape = [arr1.shape[ii] for ii in range(arr1.ndim)
                     if ii in arr1.info['rows']]
        out_names = [arr1.info['axes'][ii] for ii in range(arr1.ndim)
                     if ii in arr1.info['rows']]
        out_vect = sp.empty(out_shape)
        out_vect = make_vect(out_vect, out_names)
        n_blocks, block_shape = arr1.get_num_blocks(return_block_shape=True)
        # Make flattened veiws for the acctual matrix algebra.
        out_flat = out_vect.flat_view()
        arr2_flat = arr2.flat_view()

        for ii, block in enumerate(arr1.iter_blocks()) :
            out_flat[ii*block_shape[0]:(ii+1)*block_shape[0]] = \
                sp.dot(block, arr2_flat[ii*block_shape[1]:
                                        (ii+1)*block_shape[1]])
        return out_vect
    else :
        raise NotImplementedError("Matrix-matrix multiplication has not been "
                                  "Implemented yet.")

def partial_dot(left, right):
    """Perform matrix multiplication on some subset of the axes.
    
    This is similar to a numpy `tensordot` but it is aware of the matrix and
    vector nature of the inputs and returns appropriate objects.  It decides
    which axes to 'dot' based on the axis names.

    If a `vect` is passed, it is treated as `mat` with one row if it's the
    first arguments and a matrix with one column if it's the second.  If the
    output matrix has either only a single row or a single column, it is cast
    as a `vect`.

    This function can properly deal with block diagonal structure and axes
    sorted in any order.

    The axes in the output array as sorted such in order of block diagonal axes
    then row-only axes then col-only axes.

    Parameters
    ----------
    left : mat or vect
    right : mat or vect
    
    Returns
    -------
    out : mat or vect
        Tensor product of `left` and `right`, with any named axes appearing in
        both `left`'s columns and `right`'s rows contracted.
    """
    
    # Figure out what kind of object the inputs are.
    msg = "Inputs must be either mat or vect objects."
    if isinstance(left, mat):
        left_rows = list(left.rows)
        left_cols = list(left.cols)
    elif isinstance(left, vect):
        left_rows = []
        left_cols = range(left.ndim)
    else:
        raise TypeError(msg)
    if isinstance(right, mat):
        right_rows = list(right.rows)
        right_cols = list(right.cols)
    elif isinstance(right, vect):
        right_rows = range(right.ndim)
        right_cols = []
    else:
        raise TypeError(msg)
    # Find axes that are block diagonal make copies of the rows and cols that
    # ommits the block diagonal ones.
    left_cols_only = list(left_cols)
    left_rows_only = list(left_rows)
    right_cols_only = list(right_cols)
    right_rows_only = list(right_rows)
    left_diag = []
    left_diag_names = []
    left_col_only_names = []
    for axis in left_cols:
        if axis in left_rows:
            left_diag.append(axis)
            left_diag_names.append(left.axes[axis])
            left_rows_only.remove(axis)
            left_cols_only.remove(axis)
        else :
            left_col_only_names.append(left.axes[axis])
    right_diag = []
    right_diag_names = []
    right_row_only_names = []
    for axis in list(right_rows):
        if axis in right_cols:
            right_diag.append(axis)
            right_diag_names.append(right.axes[axis])
            right_rows_only.remove(axis)
            right_cols_only.remove(axis)
        else:
            right_row_only_names.append(right.axes[axis])
    # Divide all axes into groups based on what we are going to do with them.
    # To not be dotted.
    left_notdot = []
    # Block diagonal axis to not dot.
    left_notdot_diag = []
    # To be dotted with a normal axis.
    left_todot = []
    right_todot = []
    # To be dotted with a block diagonal axis.
    left_todot_with_diag = []
    right_todot_with_diag = []
    # Block diagonal axes to be dotted with a normal axis.
    left_todot_diag = []
    right_todot_diag = []
    # Block diagonal axes to be dotted with a block diagonal axis.
    left_todot_diag_diag = []
    right_todot_diag_diag = []
    for axis in left_cols:
        axis_name = left.axes[axis]
        if (axis_name in left_col_only_names 
            and axis_name not in right_row_only_names
            and axis_name not in right_diag_names):
            left_notdot.append(axis)
        elif (axis_name in left_diag_names 
              and axis_name not in right_row_only_names
              and axis_name not in right_diag_names):
            left_notdot_diag.append(axis)
        elif (axis_name in left_col_only_names 
              and axis_name in right_row_only_names):
            left_todot.append(axis)
            right_todot.append(right.axes.index(axis_name))
        elif (axis_name in left_diag_names 
              and axis_name in right_row_only_names):
            left_todot_diag.append(axis)
            right_todot_with_diag.append(right.axes.index(axis_name))
        elif (axis_name in left_col_only_names 
              and axis_name in right_diag_names):
            left_todot_with_diag.append(axis)
            right_todot_diag.append(right.axes.index(axis_name))
        elif axis_name in left_diag_names and axis_name in right_diag_names:
            left_todot_diag_diag.append(axis)
            right_todot_diag_diag.append(right.axes.index(axis_name))
    right_notdot = list(set(right_rows_only) - set(right_todot)
                        - set(right_todot_with_diag))
    right_notdot.sort()
    right_notdot_diag = list(set(right_diag) - set(right_todot_diag)
                             - set(right_todot_diag_diag))
    right_notdot_diag.sort()
    # Need shapes and names for all of these.
    left_notdot_shape = [left.shape[axis] for axis in left_notdot]
    left_notdot_names = [left.axes[axis] for axis in left_notdot]
    left_notdot_diag_shape = [left.shape[axis] for axis in left_notdot_diag]
    left_notdot_diag_names = [left.axes[axis] for axis in left_notdot_diag]
    left_todot_shape = [left.shape[axis] for axis in left_todot]
    left_todot_names = [left.axes[axis] for axis in left_todot]
    left_todot_with_diag_shape = [left.shape[axis] for axis in
                                  left_todot_with_diag]
    left_todot_with_diag_names = [left.axes[axis] for axis in
                                  left_todot_with_diag]
    left_todot_diag_shape = [left.shape[axis] for axis in left_todot_diag]
    left_todot_diag_names = [left.axes[axis] for axis in left_todot_diag]
    left_todot_diag_diag_shape = [left.shape[axis] for axis in
                                  left_todot_diag_diag]
    left_todot_diag_diag_names = [left.axes[axis] for axis in
                                  left_todot_diag_diag]
    left_rows_only_shape = [left.shape[axis] for axis in left_rows_only]
    left_rows_only_names = [left.axes[axis] for axis in left_rows_only]
    right_notdot_shape = [right.shape[axis] for axis in right_notdot]
    right_notdot_names = [right.axes[axis] for axis in right_notdot]
    right_notdot_diag_shape = [right.shape[axis] 
                               for axis in right_notdot_diag]
    right_notdot_diag_names = [right.axes[axis] for axis in right_notdot_diag]
    right_todot_shape = [right.shape[axis] for axis in right_todot]
    right_todot_names = [right.axes[axis] for axis in right_todot]
    right_todot_with_diag_shape = [right.shape[axis] for axis in
                                  right_todot_with_diag]
    right_todot_with_diag_names = [right.axes[axis] for axis in
                                  right_todot_with_diag]
    right_todot_diag_shape = [right.shape[axis] for axis in right_todot_diag]
    right_todot_diag_names = [right.axes[axis] for axis in right_todot_diag]
    right_todot_diag_diag_shape = [right.shape[axis] for axis in
                                  right_todot_diag_diag]
    right_todot_diag_diag_names = [right.axes[axis] for axis in
                                  right_todot_diag_diag]
    right_cols_only_shape = [right.shape[axis] for axis in right_cols_only]
    right_cols_only_names = [right.axes[axis] for axis in right_cols_only]
    # Figure out the shape, names, rows and cols of the output.
    out_shape = (left_notdot_diag_shape + left_todot_diag_diag_shape
                 + right_notdot_diag_shape + left_todot_diag_shape
                 + left_rows_only_shape + right_notdot_shape 
                 + left_notdot_shape + right_todot_diag_shape\
                 + right_cols_only_shape)
    out_names = (left_notdot_diag_names + left_todot_diag_diag_names
                 + right_notdot_diag_names + left_todot_diag_names
                 + left_rows_only_names + right_notdot_names
                 + left_notdot_names + right_todot_diag_names
                 + right_cols_only_names)
    # First add the block diagonal axes as both rows and columns.
    out_rows = range(len(left_notdot_diag) + len(left_todot_diag_diag)
                     + len(right_notdot_diag))
    out_cols = list(out_rows)
    # Now add the others.
    out_rows += range(len(out_rows), len(out_rows) + len(left_todot_diag)
                      + len(left_rows_only) + len(right_notdot))
    out_cols += range(len(out_rows), len(out_shape))
    # Output data type.
    # This is no good because it crashes for length 0 arrays.
    #out_dtype = (left.flat[[0]] * right.flat[[0]]).dtype
    # There are functions that do this in higher versions of numpy.
    out_dtype = sp.dtype(float)
    # Allowcate memory.
    out = sp.empty(out_shape, dtype=out_dtype)
    # All the block diagonal axes will be treated together. Get the global
    # shape of them.
    all_diag_shape = (left_notdot_diag_shape + left_todot_diag_diag_shape 
                      + right_notdot_diag_shape + left_todot_diag_shape 
                      + right_todot_diag_shape)
    n_diag_axes = len(all_diag_shape)
    all_diag_size = 1
    for s in all_diag_shape:
        all_diag_size *= s
    # Each of these block diagonal axes are associated with different axes in
    # the input and output arrays.  Figure out the associations for each of
    # them.  These arrays are the length of the number of diagonal axes and
    # each entry refers to an axis in that array.  If the axis does not apply
    # to that array, put None.
    out_diag_inds = []
    left_diag_inds = []
    right_diag_inds = []
    tmp_n_out_axes_passed = 0
    for ii in range(len(left_notdot_diag)):
        left_diag_inds.append(left_notdot_diag[ii])
        right_diag_inds.append(None)
        out_diag_inds.append(ii)
    tmp_n_out_axes_passed += len(left_notdot_diag)
    for ii in range(len(left_todot_diag_diag)):
        left_diag_inds.append(left_todot_diag_diag[ii])
        right_diag_inds.append(right_todot_diag_diag[ii])
        out_diag_inds.append(ii + tmp_n_out_axes_passed)
    tmp_n_out_axes_passed += len(left_todot_diag_diag)
    for ii in range(len(right_notdot_diag)):
        left_diag_inds.append(None)
        right_diag_inds.append(right_notdot_diag[ii])
        out_diag_inds.append(ii + tmp_n_out_axes_passed)
    tmp_n_out_axes_passed += len(right_notdot_diag)
    for ii in range(len(left_todot_diag)):
        left_diag_inds.append(left_todot_diag[ii])
        right_diag_inds.append(right_todot_with_diag[ii])
        out_diag_inds.append(ii + tmp_n_out_axes_passed)
    tmp_n_out_axes_passed += (len(left_todot_diag) + len(left_rows_only) 
                              + len(right_notdot) + len(left_notdot))
    for ii in range(len(right_todot_diag)):
        left_diag_inds.append(left_todot_with_diag[ii])
        right_diag_inds.append(right_todot_diag[ii])
        out_diag_inds.append(ii + tmp_n_out_axes_passed)
    # Once we index all the diagonal axes, the ones we want to dot will be in
    # the wrong place.  Find the location of the axes we'll need in the new
    # array.
    left_sliced_rows_only = list(left_rows_only)
    left_sliced_notdot = list(left_notdot)
    left_sliced_todot = list(left_todot)
    for diag_axis in left_diag_inds:
        if not diag_axis is None:
            for ii in range(len(left_rows_only)):
                if diag_axis < left_rows_only[ii]:
                    left_sliced_rows_only[ii] -= 1
            for ii in range(len(left_notdot)):
                if diag_axis < left_notdot[ii]:
                    left_sliced_notdot[ii] -= 1
            for ii in range(len(left_todot)):
                if diag_axis < left_todot[ii]:
                    left_sliced_todot[ii] -= 1
    right_sliced_cols_only = list(right_cols_only)
    right_sliced_notdot = list(right_notdot)
    right_sliced_todot = list(right_todot)
    for diag_axis in right_diag_inds:
        if not diag_axis is None:
            for ii in range(len(right_cols_only)):
                if diag_axis < right_cols_only[ii]:
                    right_sliced_cols_only[ii] -= 1
            for ii in range(len(right_notdot)):
                if diag_axis < right_notdot[ii]:
                    right_sliced_notdot[ii] -= 1
            for ii in range(len(right_todot)):
                if diag_axis < right_todot[ii]:
                    right_sliced_todot[ii] -= 1
    # Once we slice the arrays to get one block, we will permute the axes to be
    # in the proper order for dotting.
    left_sliced_permute = (left_sliced_rows_only + left_sliced_notdot
                           + left_sliced_todot)
    right_sliced_permute = (right_sliced_todot + right_sliced_notdot
                            + right_sliced_cols_only)
    # Then we'll reshape both into 2D arrays (matricies).
    left_sliced_reshape = (sp.prod(left_rows_only_shape + left_notdot_shape),
                           sp.prod(left_todot_shape))
    right_sliced_reshape = (sp.prod(right_todot_shape),
                            sp.prod(right_notdot_shape +
                                    right_cols_only_shape))
    # After the dot, we will neet to reshape back.
    out_sliced_reshape = tuple(left_rows_only_shape + left_notdot_shape
                               + right_notdot_shape + right_cols_only_shape)
    # And finally we will need to permute the out axes.
    out_sliced_permute = range(len(left_rows_only))
    out_sliced_permute += range(len(left_rows_only) + len(left_notdot_shape),
                                len(left_rows_only) + len(left_notdot_shape) 
                                + len(right_notdot))
    out_sliced_permute += range(len(left_rows_only), len(left_rows_only)
                                + len(left_notdot))
    out_sliced_permute += range(len(out_sliced_reshape)- len(right_cols_only),
                                len(out_sliced_reshape))
    # Create an index for each of left, right and out.
    left_slice = [slice(None)] * left.ndim
    right_slice = [slice(None)] * right.ndim
    out_slice = [slice(None)] * out.ndim
    # Flags for corner cases of 0D arrays.
    left_scalar_flag = False
    right_scalar_flag = False
    # Now we loop over all the block diagonal axes.
    for ii in xrange(all_diag_size):
        # Figure out exactly which blocks we are dealing with.
        tmp_ii = ii
        for kk in xrange(n_diag_axes - 1, -1, -1):
            this_index = tmp_ii % all_diag_shape[kk]
            out_slice[out_diag_inds[kk]] = this_index
            if not left_diag_inds[kk] is None:
                left_slice[left_diag_inds[kk]] = this_index
            if not right_diag_inds[kk] is None:
                right_slice[right_diag_inds[kk]] = this_index
            tmp_ii = tmp_ii // all_diag_shape[kk]
        this_left = left[tuple(left_slice)]
        this_right = right[tuple(right_slice)]
        # Permute and reshape the axes so the fast dot function can be used.
        # Corner case, if this_left or this_right are scalars (0D arrays).
        # TODO: Really these flags can be set outside the loop with a few more
        # lines of code.
        if ii == 0:
            if this_left.ndim == 0:
                left_scalar_flag = True
            if this_left.ndim == 0 and this_right.ndim == 0:
                right_scalar_flag = True
        if not left_scalar_flag:
            this_left = this_left.transpose(left_sliced_permute)
            this_left = sp.reshape(this_left, left_sliced_reshape)
        if not right_scalar_flag:
            this_right = this_right.transpose(right_sliced_permute)
            this_right = sp.reshape(this_right, right_sliced_reshape)
        # Dot them.
        if left_scalar_flag or right_scalar_flag:
            this_out = this_left * this_right
        else:
            this_out = sp.dot(this_left, this_right)
        # Reshape, permute and copy to the output.
        if not (left_scalar_flag and right_scalar_flag):
            this_out.shape = out_sliced_reshape
            this_out = this_out.transpose(out_sliced_permute)
        out[tuple(out_slice)] = this_out
    # XXX There is a bug where this crashes for certain cases if these lines
    # are before the loop.
    if not out_rows or not out_cols:
        out = make_vect(out, out_names)
    else:
        out = make_mat(out, axis_names=out_names, row_axes=out_rows,
                       col_axes=out_cols)
    return out

def empty_like(obj) :
    """Create a new algebra object with uninitialized data but otherwise the
    same as the passed object."""

    out = sp.empty_like(obj)
    return as_alg_like(out, obj)

def zeros_like(obj) :
    """Create a new algebra object full of zeros but otherwise the same 
    as the passed object."""

    out = sp.zeros_like(obj)
    return as_alg_like(out, obj)

def ones_like(obj) :
    """Create a new algebra object full of zeros but otherwise the same 
    as the passed object."""

    out = sp.ones_like(obj)
    return as_alg_like(out, obj)

def as_alg_like(array, obj):
    """Cast an array as an algebra object similar to the passed object.
    
    Parameters
    ----------
    array : numpy array
        Array to be cast
    obj : alg_object
        Algebra object from which propertise should be copied.
    """
    
    if not isinstance(obj, alg_object):
        raise TypeError("Object to mimic must be an `alg_object`.")
    out = array
    out = info_array(out)
    out.info = dict(obj.info)
    if isinstance(obj, vect) :
        out = make_vect(out)
    elif isinstance(obj, mat) :
        out = make_mat(out)
    else :
        raise TypeError("Expected `obj` to be an algebra mat or vect.")
    
    return out


# TODO: These need scipy standard documentation.
def array_summary(array, testname, axes, meetall=False, identify_entries=True):
    """helper function for summarizing arrays
    meetall: prints those entries for which all values in the slice meet the
    criteria (normal behavior is print all entries where _any_ value in the
    slice meets the criteria
    identify_entries: prints entries meeting the criteria
    """
    total_matching = array.sum()
    if total_matching != 0:
        print testname + "s:"
        match_count = np.apply_over_axes(np.sum, array, axes)
        print match_count.flatten()
        if identify_entries:
            if meetall:
                arrayshape = array.shape
                subarray_size = reduce(operator.mul,
                                       [arrayshape[i] for i in axes])
                print "with all " + testname + "s: " + \
                      repr(np.where(match_count.flatten() == subarray_size))
            else:
                print "has " + testname + ": " + \
                      repr(np.where(match_count.flatten() != 0))
            print "total " + testname + "s: " + repr(total_matching)
        print "-" * 80
    else:
        print "There are no " + testname + " entries"

def compressed_array_summary(array, name, axes=[1, 2], extras=False):
    """print various summaries of arrays compressed along specified axes"""

    print "-" * 80
    print "array property summary for " + name + ":"
    array_summary(np.isnan(array), "nan", axes)
    array_summary(np.isinf(array), "inf", axes)
    array_summary((array == 0.), "zero", axes, meetall=True)
    array_summary((array < 0.), "negative", axes, identify_entries=False)

    if extras:
        sum_nu = np.apply_over_axes(np.sum, array, axes)
        min_nu = np.apply_over_axes(np.min, array, axes)
        max_nu = np.apply_over_axes(np.max, array, axes)
        print sum_nu.flatten()
        print min_nu.flatten()
        print max_nu.flatten()
    print ""


def cartesian(arrays, out=None):
    """
    SO: using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
    Generate a cartesian product of input arrays. ~5x faster than itertools

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:, 0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1: ], out=out[0: m, 1: ])
        for j in xrange(1, arrays[0].size):
            out[j * m: (j + 1) * m, 1: ] = out[0: m, 1: ]

    return out


def roll_zeropad(a, shift, axis=None):
    """
    SO: python-numpy-roll-with-padding
    Roll array elements along a given axis.

    Elements off the end of the array are treated as zeros.

    Parameters
    ----------
    a : array_like
        Input array.
    shift : int
        The number of places by which elements are shifted.
    axis : int, optional
        The axis along which elements are shifted.  By default, the array
        is flattened before shifting, after which the original
        shape is restored.

    Returns
    -------
    res : ndarray
        Output array, with the same shape as `a`.

    See Also
    --------
    roll     : Elements that roll off one end come back on the other.
    rollaxis : Roll the specified axis backwards, until it lies in a
               given position.

    Examples
    --------
    >>> x = np.arange(10)
    >>> roll_zeropad(x, 2)
    array([0, 0, 0, 1, 2, 3, 4, 5, 6, 7])
    >>> roll_zeropad(x, -2)
    array([2, 3, 4, 5, 6, 7, 8, 9, 0, 0])

    >>> x2 = np.reshape(x, (2,5))
    >>> x2
    array([[0, 1, 2, 3, 4],
           [5, 6, 7, 8, 9]])
    >>> roll_zeropad(x2, 1)
    array([[0, 0, 1, 2, 3],
           [4, 5, 6, 7, 8]])
    >>> roll_zeropad(x2, -2)
    array([[2, 3, 4, 5, 6],
           [7, 8, 9, 0, 0]])
    >>> roll_zeropad(x2, 1, axis=0)
    array([[0, 0, 0, 0, 0],
           [0, 1, 2, 3, 4]])
    >>> roll_zeropad(x2, -1, axis=0)
    array([[5, 6, 7, 8, 9],
           [0, 0, 0, 0, 0]])
    >>> roll_zeropad(x2, 1, axis=1)
    array([[0, 0, 1, 2, 3],
           [0, 5, 6, 7, 8]])
    >>> roll_zeropad(x2, -2, axis=1)
    array([[2, 3, 4, 0, 0],
           [7, 8, 9, 0, 0]])

    >>> roll_zeropad(x2, 50)
    array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0]])
    >>> roll_zeropad(x2, -50)
    array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0]])
    >>> roll_zeropad(x2, 0)
    array([[0, 1, 2, 3, 4],
           [5, 6, 7, 8, 9]])

    """
    a = np.asanyarray(a)
    if shift == 0: return a
    if axis is None:
        n = a.size
        reshape = True
    else:
        n = a.shape[axis]
        reshape = False
    if np.abs(shift) > n:
        res = np.zeros_like(a)
    elif shift < 0:
        shift += n
        zeros = np.zeros_like(a.take(np.arange(n-shift), axis))
        res = np.concatenate((a.take(np.arange(n-shift,n), axis), zeros), axis)
    else:
        zeros = np.zeros_like(a.take(np.arange(n-shift,n), axis))
        res = np.concatenate((zeros, a.take(np.arange(n-shift), axis)), axis)
    if reshape:
        return res.reshape(a.shape)
    else:
        return res
