"""This module provides vector and matrix interfaces tylored to the needs of
intensity mapping analysis.

The considerations for intensity mapping analysis can be summarized as follows:
    - Elements in data vectors often have a significant amount of assotiated
      meta data.  For instance in map space, a data point is a pixel and it has
      pointing and frequency information that needs to be stored with it.
    - Matricies can generally refer to data vectors for meta data.
    - Data vectors can normally be held in the memory of a single machine, but
      matricies cannot.
    - Matricies are occationally Block diagonal and we would like to take
      advantage of this.
    - We will evenetually need to do operations massivly in parralelle with
      SCALLAPACK.

To satisfy these needs, this module provides:
    - A data structure for vectors and infrastructure to read and write them to
      disk.  This is implemented using numpy record arrays, with an added
      header. Memory mapping is fully supported.
    - Infrastructure for using numpy array objects as matricies, with memory
      mapping a being a consideration throughout.
    - To take advantage of block diagonal matricies, all matricies are 3 
      dimentional arrays, with the 0th dimention indicating the block.  The 0th
      dimension may be of lenth 1 (for a dense array).
    - Infrastructure for maniputlating these formats and performing standard
      operations.
"""

import os

import scipy as sp
import numpy.lib.format as npfor
from numpy.lib.utils import safe_eval

import kiyopy.custom_exceptions as ce


# ---- Array classes for holding vectors and matricies. -----------------------

class info_array(sp.ndarray) :
    """A standard numpy ndarray object with a dictionary for holding extra info.

    This class should work exactly the same as a numpy ndarray object but has an
    attribut named info, which is a dictionary.

    Initialize this class by passing it a ndarray and optionally a dictionary
    which will be used as info.
    """

    def __new__(cls, input_array, info=None):
        # Input array is an already formed ndarray instance.
        # We first cast to be our class type.
        obj = sp.asarray(input_array).view(cls)
        # Add the new attribute to the created instance.
        if info is None :
            info = {}
        obj.info = info
        # Finally, we must return the newly created object.
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.info = getattr(obj, 'info', {})

class info_memmap(sp.memmap) :
    """A standard numpy memmap object with a dictionary for holding extra info.

    This class should work exactly the same as a numpy memmap object but has an
    attribut named info, which is a dictionary.
    
    Initialization arguments :
    marray : A numpy memmap.
    info=None : A dictionary to be used as the info metadata dictionary. 
    metafile=None : String filename to write the metadata to.  In some
        versions of numpy, the metadata will be written to file even if the
        memmap is in read only mode.  To avoid this pass metafile=None.
    """

    def __new__(cls, input_array, info=None, metafile=None):
        # Input array is an already formed ndarray instance.
        # We first cast to be our class type.
        if not isinstance(input_array, sp.memmap) :
            raise TypeError("info_memmaps can only be initialized off of "
                            "numpy memmaps.")
        obj = input_array.view(cls)
        # Add the new attribute to the created instance.
        if info is None :
            info = {}
        obj.info = info
        obj.metafile = metafile
        # Finally, we must return the newly created object.
        return obj

    def __array_finalize__(self, obj):
        sp.memmap.__array_finalize__(self, obj)
        if obj is None: return
        self.info = getattr(obj, 'info', {})
        self.metafile = getattr(obj, 'metafile', None)

    def flush(self) :
        """Flush changes to disk.

        This method saves the info dictionary to metafile and then calls the
        flush method from the numpy memmap.
        """
        
        # Prior to numpy 1.5, we can't get teh mode, so just assume we are
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
        # Now flush the actual memmap.
        sp.memmap.flush(self)

    def __del__(self) :
        self.flush()
        sp.memmap.__del__(self)

def assert_info(array) :
    """Check if passed array is an info_array or info_memmap.

    Raises a ValueError if check fails.
    """
    if not (isinstance(array, info_array) or isinstance(array, info_memmap)) :
        raise TypeError("Array is not an algebra.info_array or "
                         "algebra.info_memmap.")


# ---- Function for reading and writing the above to and from file. -----------

def open_memmap(filename, mode='r+', dtype=None, shape=None,
                fortran_order=False, version=(1,0), metafile=None) :
    """Open a file and memory map it to an info_memmap object.

    This is similar to the numpy.lib.format.openmemmap() function but also
    deals with the meta data dictionary, which is read and written from a
    meta data file.
    """
    
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
    """Open a .npy file and load it into memory as a info_aray.
    
    Similar to the numpy.load funciton.  Does not support memory
    mapping.
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

def save(file, info_array, metafile=None) :
    """Save a info array to a .npy file and a metadata file.
    """
    # Make sure that the meta data will be representable as a string.
    infostring = repr(info_array.info)
    try:
        safe_eval(infostring)
    except SyntaxError :
        raise ce.DataError("Array info not representable as a string.")
    # Save the array in .npy format.
    sp.save(file, info_array)
    # Figure out what the filename for the meta data should be.
    if metafile is None :
        try :
            fname = file.name
        except AttributeError :
            fname = file
        metafile = fname + ".meta"
    #Save the meta data.
    info_fid = open(metafile, 'w')
    try :
        info_fid.write(infostring)
    finally :
        info_fid.close()


# ---- Functions for manipulating above arrays as matircies and vectors. ------

def _set_type_axes(array, type, axis_names) :
    """Sets the array.info['type'] and array.info[axes] metadata and does some
    checks.  Used in make_vect and make_mat.
    """
    
    assert_info(array)

    if axis_names is None :
        axes = (None,)*array.ndim
    elif len(axis_names) != array.ndim :
        raise ValueError("axis_names parameter must be a sequence of length "
                         "vect.ndims")
    else :
        for name in axis_names :
            if (not isinstance(name, str)) and (not name is None) :
                raise TypeError("Invalid axis name.")
        axes = tuple(axis_names)
    array.info['type'] = type
    array.info['axes'] = axes

def make_vect(vect, axis_names=None) :
    """Identify an info array as a vector.

    When passed an info_array or info_memmap, set the nessisary metadata to
    identify it as a vector.  Optionally pass names for the axes in the array.
    The axis names should be passed as a tuple with length vect.ndims.  The
    entries should be either string names or None.
    """
    
    _set_type_axes(vect, 'vect', axis_names)

def make_mat(mat, row_axes=None, col_axes=None, axis_names=None) :
    """Identify an info array as a matrix.

    When passed an info_array or info_memmap, set the nessisary metadata to
    identify it as a Matrix.
    Other arguments:
    row_axes and col_axes: Unless the array is exactly 2D, the row_axes and
                col_axes parameters must be specified.  They are sequences of
                integers telling which axes to identify as the matrix rows and
                which to identify as matrix columns.  All axes of the array
                must be accounted.  An axis may be identified as both in which
                case the matrix is block diagonal in some space.
    axis_names: name for the axes in the array. The axis names should be passed
                as a tuple with length vect.ndims.  The entries should be either
                string names or None.
    """

    # Check the row and col arguments.
    if row_axes is None and col_axes is None and mat.ndim==2 :
        row_axes = (0,)
        col_axes = (1,)
    else :
        # Both parameters had better be sequences of integers.
        for ind in row_axes :
            if not ind in range(mat.ndim) :
                raise ValueError("Invalid row axes.")
        for ind in col_axes :
            if not ind in range(mat.ndim) :
                raise ValueError("Invalid col axes")
        # Make sure each axis is spoken for.
        for ii in range(mat.ndim) :
            if (not ii in row_axes) and (not ii in col_axes) :
                raise ValueError("Every axis must be identified varying over "
                                 "as the matrix row, column or both.")

    _set_type_axes(mat, 'mat', axis_names)
    mat.info['rows'] = tuple(row_axes)
    mat.info['cols'] = tuple(col_axes)

def mat_shape(mat) :
    """Returns the external dimension of an array like object.

    The passed array should be a 3D array.  The zeroth dimension gives the 
    block of the assumed block diagonal matrix.  The external matrix shape is
    the shape seen by the matrix interface and is a 2 element tuple.  It is
    given by (mat.shape[0]*mat.shape[1], mat.shape[0]*mat.shape[2]).

    This function also checks that mat is a valid matrix for this module and
    this function may be used to check this.
    """
    assert(0)

    if not isinstance(mat, sp.ndarray) :
        raise TypeError("Argument is not an array.")
    if mat.ndim != 3 :
        raise TypeError("Arrays must be 3 dimensional to be compatible with "
                        "this module as a matrix.")

    return (mat.shape[0]*mat.shape[1], mat.shape[0]*mat.shape[2])

    
