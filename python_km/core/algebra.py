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

import scipy as sp
import numpy.lib.format as npfor
from numpy.lib.utils import safe_eval

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

    Notes
    -----
    All new from template array creation operations make a copy of the metadata
    not a reference.  To get a reference you need to explicitly make a call to
    the `view` function.

    See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html for more
    information
    
    See Also
    --------
    info_memmap : Analogous class to this one with data stored on disk.
    vect_array : Vector object based on this class.
    mat_array " Matrix object based on this class.
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
    
    def __deepcopy__(self, copy):
        """Not implemented, raises an exception."""
        raise NotImeplementedError("Deep copy won't work.")

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

    Notes
    -----
    All new from template array creation operations make a copy of the metadata
    not a reference.  To get a reference you need to explicitly make a call to
    the `view` function.  Also, the metafile is set to None on new from
    template oberations.  
    
    See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html for more
    information

    See Also
    --------
    info_array : Similar class with data stored in memory.
    vect_memmap : Vector object based on this class.
    mat_memmap : Matrix object based on this class.
    open_memmap : Open a file on disk as an info_memmap.
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
        raise NotImeplementedError("Deep copy w on't work.")

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
    both `mat` objects and `vect` objects."""
    
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

        for axis in self.axes :
            if axis in alg_obj.axes :
                try :
                    centre = alg_obj.info[axis + '_centre']
                    delta = alg_obj.info[axis + '_delta']
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
        axis_name : str
            Name of the axis to be calculated.  `axis_name` must occur in the
            `axes` attribute.

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
        
        len = self.shape[self.axes.index(axis_name)]
        return (self.info[axis_name + '_delta']*(sp.arange(len) - len//2) 
                + self.info[axis_name + '_centre'])

    def __array_finalize__(self, obj) :
        self.info_base.__array_finalize__(self, obj)
        if (not obj is None) and (self.shape != obj.shape) :
            self.__class__ = self.info_base


#### Vector class definitions ####

class vect(alg_object) :
    """Base class for vectors. 
    
    This is only half a class.  A complete vector class is created by making
    a class that inherits from both this class and from info_array or
    info_memmap, for example the classes vect_array and vect_memmap.
    """
    
    def flat_view(self) :
        """Returns a view of the vector that has been flattend.

        The view will be cast as a scipy.ndarray and its shape will be
        (self.size, ).  It is a view, so writing to it will write back to the
        origional vector object.

        Returns
        -------
        flat_view : np.ndarray
            A veiw of `self` as an ndarray, flattened to 1D.
        """
        
        flat = self.view(sp.ndarray)
        flat.shape = (self.size, )
        return flat

    def mat_shape(self) :
        """Get the shape of the represented matrix (vector)."""
        _check_axis_names(self)
        return (self.size,)

    __array_priority__ = 2.0

def _vect_class_factory(base_class) :
    """Internal class factory for making a vector class that inherits from
    either info_array or info_memmap."""

    if (not base_class is info_array) and (not base_class is info_memmap) :
        raise TypeError("Vectors inherit from info arrays or info memmaps.")

    class vect_class(vect, base_class) :
        """Multidimentional array interpreted as a vector.
        
        This class gets most of its functionality from the numpy ndarray class.
        In addition it provides support for orgainizing it's data as a vector.
        This class comes in two flavours: `vect_array` and `vect_memmap`
        depending on whether the array is stored in memory or on disk.

        The vector representation of the array is the flattend array.

        One of the features that `vect`s implement is named axes.  This allows
        maps to carry axis information with them, amonge other things.  For
        `vect`s axis names must be unique (This is not true of `mat`s).

        Parameters
        ----------
        input_array : info_array (for vect_array) or info_memmap (for
                      vect_memmap)
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
        The `axes` attribute is acctually stored in the `info_array`'s info
        dictionary.  This is just an implimentation detail.
        """
        
        info_base = base_class

        # Only define methods here that have to refer to base_class.
        # Everything else can go into vect.
        def __new__(cls, input_array, axis_names=None) :
            
            if not isinstance(input_array, base_class) :
                raise ValueError("Array to convert must be instance of " +
                                 str(base_class))
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
                base_class.__setattr__(self, name, value)

        def __getattr__(self, name) :
            if name == 'axes' :
                return self.info['axes']
            else :
                # Since numpy uses __get_attribute__ not __getattr__, we should
                # raise an Attribute error.
                raise AttributeError("Attribute " + name + " not found.")

    return vect_class

vect_array = _vect_class_factory(info_array)
vect_array.__name__ = 'vect_array'
vect_memmap = _vect_class_factory(info_memmap)
vect_memmap.__name__ = 'vect_memmap'
    
def make_vect(array, axis_names=None) :
    """Do what ever it takes to make an vect out of an array.
    
    Convert any class that can be converted to a vect (array, info_array,
    memmap, info_memmap) to the appropriate vect object (vect_array,
    vect_memmap).

    This convieiance function just simplifies the constructor heirarchy.  
    Normally to get an vect out of an array, you would need to construct an
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
    """Base class for matricies.
    
    This is only half a class.  A complete matris class is created by making
    a class that inherits from both this class and from info_array or
    info_memmap, for example the classes mat_array and mat_memmap.
    """
    
    __array_priority__ = 3.0

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
        >>> for index in mat.iter_row_axes :
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
        >>> for index in mat.iter_col_axes :
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

def _mat_class_factory(base_class) :
    """Internal class factory for making a matrix class that inherits from
    either info_array or info_memmap."""

    if (not base_class is info_array) and (not base_class is info_memmap) :
        raise TypeError("Matrices inherit from info arrays or info memmaps.")

    class mat_class(mat, base_class) :
        """Multidimentional array interpreted as a matrix.
        
        This class gets most of its functionality from the numpy ndarray class.
        In addition it provides support for orgainizing it's data as a vector.
        This class comes in two flavours: `mat_array` and `mat_memmap`
        depending on whether the array is stored in memory or on disk.

        To make the assotiation between a multidimentional array and a matrix,
        each axis of the array must be identified as varying over either the
        rows or colums of a matrix.  For instance the shape of the array could
        be (3, 5, 7).  We could identify the first axis as a row axis and the
        second two as column axis in which case the matrix would have 3 rows
        and 35 colums.  We generally make the rows axes left of the columns and
        many algorithms assume this.  It is also possible for an axis to be
        identified as both a row and a column, in which case the matrix is
        block diagonal over that axis.  Generally the block diagonal axes are
        the left most.

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
        The `axes`, `rows` and `cols` attributes are acctually stored in the 
        `info_array`'s info dictionary.  This is just an implimentation detail.

        See Also
        --------
        vect_array, vect_memmap : Vector classes.
        make_mat : Funciton that casts any array as a matrix.
        info_array, info_memmap : Base classes that handle meta data.
        """

        info_base = base_class
        
        def __new__(cls, input_array, row_axes=None, col_axes=None, 
                    axis_names=None) :
            
            if not isinstance(input_array, base_class) :
                raise ValueError("Array to convert must be instance of " +
                                 str(base_class))
            
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
                if row_axes is None and col_axes is None and input_array.ndim==2 :
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
                base_class.__setattr__(self, name, value)

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
        raise ValueError("Matrix dimensions incompatible for matrix ",
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

def zeros_like(obj) :
    """Create a new algebra object full of zeros but otherwise the same 
    as the passed object."""

    out = sp.zeros_like(obj)

    if isinstance(obj, vect) :
        out = make_vect(out)
        out.info = dict(obj.info)
    elif isinstance(obj, mat) :
        out = info_array(out)
        out.info = dict(obj.info)
        out = make_mat(out)
    else :
        raise TypeError("Expected an algebra mat or vect.")
    
    return out
