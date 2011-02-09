"""This module provides vector and matrix interfaces tailored to the needs of
intensity mapping analysis.

The considerations for intensity mapping analysis can be summarized as follows:
    - Elements in data vectors often have a significant amount of associated
      meta data.  For instance in map space, a data point is a pixel and it has
      pointing and frequency information that needs to be stored with it.
    - Matricies can generally refer to data vectors for meta data.
    - Data vectors can normally be held in the memory of a single machine, but
      matrices cannot.
    - Matricies are occationally Block diagonal and we would like to take
      advantage of this.
    - We will evenetually need to do operations massively in parallel with
      SCALLAPACK.

To satisfy these needs, this module provides:
    - A data structure for vectors and infrastructure to read and write them to
      disk.  This is implemented using numpy record arrays, with an added
      header. Memory mapping is fully supported.
    - Infrastructure for using numpy array objects as matrices, with memory
      mapping a being a consideration throughout.
    - Infrastructure for manipulating these formats and performing standard
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
    attribute named info, which is a dictionary.

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
    attribute named info, which is a dictionary.
    
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


# ---- Functions for reading and writing the above to and from file. -----------

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
    
    Similar to the numpy.load function.  Does not support memory
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
    # Save the meta data.
    info_fid = open(metafile, 'w')
    try :
        info_fid.write(infostring)
    finally :
        info_fid.close()


# ---- Functions for manipulating above arrays as matrices and vectors. ------

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
    array.info['axes'] = axes

def _check_axis_names(array, axis_names=None) :
    """Checks that axis names sequence is valid for array."""

    if axis_names is None :
        axis_names = array.info['axes']

    if len(axis_names) != array.ndim :
        raise ValueError("axis_names parameter must be a sequence of length "
                         "vect.ndims")
    else :
        for name in axis_names :
            if (not isinstance(name, str)) and (not name is None) :
                raise TypeError("Invalid axis name.")

class vect(object) :
    """Base class for vectors.
    
    This is only half a class.  A complete vector class is created by making a
    class that inherits from both this class and from info_array or
    info_memmap, for example the classes vect_array and vect_memmap.
    """
    pass
    
def _vect_class_factory(base_class) :
    """Internal class factory for making a vector class that inherits from
    either info_array or info_memmap."""

    if (not base_class is info_array) and (not base_class is info_memmap) :
        raise TypeError("Vectors inherit from info arrays or info memmaps.")

    class vect_class(base_class, vect) :
        """Vector class for this module.
        
        When passed an info_array or info_memmap, set the necessary metadata
        to identify it as a vector.  Optionally pass names for the axes in the
        array. The axis names should be passed as a tuple with length 
        vect.ndims.  The entries should be either string names or None.
        """
        
        def __new__(cls, input_array, axis_names=None) :
            
            if not isinstance(input_array, base_class) :
                raise ValueError("Array to convert must be instance of " +
                                 str(base_class))
            obj = input_array.view(cls)
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
            elif name == 'mat_shape' :
                _check_axis_names(self)
                return (self.size,)
            else :
                raise AttributeError("Attribute " + name + " not found.")

    return vect_class

vect_array = _vect_class_factory(info_array)
vect_memmap = _vect_class_factory(info_memmap)
    
def make_vect(array, axis_names=None) :
    """Do what ever it takes to make an vect out of an array."""

    if isinstance(array, sp.memmap) :
        if not isinstance(array, info_memmap) :
            array = info_memmap(array)
        return vect_memmap(array, axis_names)
    elif sp.isarray(array, sp.ndarray) :
        if not isinstance(array, info_array) :
            array = info_array(array)
        return vect_array(array, asix_names)
    else :
        raise TypeError("Object cannot be converted to a vector.")


class mat(object) :
    """Base class for matricies.
    
    This is only half a class.  A complete matris class is created by making a
    class that inherits from both this class and from info_array or
    info_memmap, for example the classes mat_array and mat_memmap.
    """
    
    def check_rows_cols(self) :
        """Check that rows and cols are valid for the matrix."""
        _check_rows_cols(self)
        
def _mat_class_factory(base_class) :
    """Internal class factory for making a matrix class that inherits from
    either info_array or info_memmap."""

    if (not base_class is info_array) and (not base_class is info_memmap) :
        raise TypeError("Matrices inherit from info arrays or info memmaps.")

    class mat_class(base_class, mat) :
        """Matrix class for this module."""
        # TODO: Write a long doc string.
                
        def __new__(cls, input_array, row_axes=None, col_axes=None, 
                    axis_names=None) :
            
            if not isinstance(input_array, base_class) :
                raise ValueError("Array to convert must be instance of " +
                                 str(base_class))

            if row_axes is None and col_axes is None and input_array.ndim==2 :
                row_axes = (0,)
                col_axes = (1,)
            else :
                _check_rows_cols(input_array, row_axes, col_axes)
            
            obj = input_array.view(cls)
            
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
                        raise ValueError("Invalid row axes.")
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
            elif name == 'mat_shape' :
                self.check_rows_cols()
                _check_axis_names(self)
                nrows = 1
                for axis in self.info['rows'] :
                    nrows *= self.shape[axis]
                ncols = 1
                for axis in self.info['cols'] :
                    ncols *= self.shape[axis]
                return (nrows, ncols)
            else :
                raise AttributeError("Attribute " + name + "not found.")

    return mat_class

mat_array = _mat_class_factory(info_array)
mat_memmap = _mat_class_factory(info_memmap)








def make_mat(mat, row_axes=None, col_axes=None, axis_names=None) :
    """Identify an info array as a matrix.

    When passed an info_array or info_memmap, set the necessary metadata to
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
        _check_rows_cols(mat, row_axes, col_axes)
    _set_type_axes(mat, 'mat', axis_names)
    mat.info['rows'] = tuple(row_axes)
    mat.info['cols'] = tuple(col_axes)

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

def get_shape(arr) :
    """Returns the external dimensions of a matrix or vector.

    Given an info array that has been identified as a matrix pr vector
    (using make_mat or make_vect),
    this function returns the dimensions of that matrix.  It also performs
    several checks to ensure that the matrix is valid.  Calling this function
    is the preferred way to check the validity of a matrix.
    """

    assert_info(arr)
    _check_axis_names(arr)
    
    if arr.info['type'] == 'vect' :
        return (sp.size(arr),)
    elif arr.info['type'] == 'mat' :
        _check_rows_cols(arr)
        nrows = 1
        for axis in arr.info['rows'] :
            nrows *= arr.shape[axis]
        ncols = 1
        for axis in arr.info['cols'] :
            ncols *= arr.shape[axis]
        return (nrows, ncols)
    else :
        raise ValueError("Unrecognized algebra type.")

def assert_axes_ordered(mat) :
    """Enforces a specific ordering to the matrix row and column axis
    associations.
    """

    rows = mat.info['rows']
    cols = mat.info['cols']
    r_ind = len(rows) - 1
    c_ind = len(cols) - 1
    in_rows = False
    in_cols = True

    for axis in range(mat.ndim-1, -1, -1) :
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

def get_num_blocks(mat, return_block_shape=False, return_n_axes_diag=False) :
    """Get the number of blocks in a block diagonal matrix."""
    
    shape = get_shape(mat)
    # Current algorithm assumes specific format.
    assert_axes_ordered(mat)
    
    diag_axes = [ii for ii in range(mat.ndim) if ii in mat.info['rows'] and 
                 ii in mat.info['cols']]
    num_blocks = sp.prod([mat.shape[ii] for ii in diag_axes])
    if return_block_shape and return_n_axes_diag :
        return num_blocks, (shape[0]/num_blocks, shape[1]/num_blocks), \
               len(diag_axes)
    elif return_block_shape :
        return num_blocks, (shape[0]/num_blocks, shape[1]/num_blocks)
    elif return_n_axes_diag :
        return num_blocks, len(diag_axes)
    else :
        return num_blocks

class iterate_blocks(object) :
    """An iterator that iterates over the blocks of an matrix."""
    
    def __init__(self, mat) :
        self.mat = mat
        self.n_blocks, self.block_shape, self.n_axes_diag = \
            get_num_blocks(mat, True, True)
        self.ii = 0

    def __iter__(self) :
        return self

    def next(self) :
        if self.ii >= self.n_blocks :
            raise StopIteration
        else :
            # Find the indices for this block.
            array_index = ()
            tmp_block_num = self.ii
            self.ii += 1
            for jj in range(self.n_axes_diag -1, -1, -1) :
                array_index = (tmp_block_num%self.mat.shape[jj],) + array_index
                tmp_block_num = tmp_block_num//self.mat.shape[jj]
            # Return the data.
            return sp.reshape(self.mat[array_index], self.block_shape)
        

def expand_mat(mat) :
    """Calculates expanded matrix in 2 dimensional form.

    Takes an arbitrary matrix and returns the expanded version of it, as matrix
    with internal array dimensions of shape(mat).  If the original matrix has
    efficiency from any block diagonal structure, this is lost in the returned
    matrix.
    """
    
    # Also verifies the validity of the matrix.
    shape = get_shape(mat)
    # Current algorithm assumes specific format.
    assert_axes_ordered(mat)
    # Allocate memory.
    out_mat = sp.zeros(shape, dtype=mat.dtype)
    out_mat = info_array(out_mat)
    make_mat(out_mat)

    # Figure out how many axes are in both row and col (and therefore block
    # diagonal).
    n_blocks, block_shape = get_num_blocks(mat, True, False)
    
    # Loop over the blocks and assign data.
    for ii, mat_block in enumerate(iterate_blocks(mat)) :
        # Figure out where this block starts.
        row_start = ii*block_shape[0]
        col_start = ii*block_shape[1]
        out_mat[row_start:row_start + block_shape[0], 
                col_start:col_start + block_shape[1]] = mat_block
    return out_mat

def dot(arr1, arr2) :
    """Perform matrix multiplication."""

    shape1 = get_shape(arr1)
    shape2 = get_shape(arr2)

    if shape1[-1] != shape2[0] :
        raise ValueError("Matrix dimensions incompatible for matrix ",
                         "multiplication.")
    
    # Matrix-vector product case.
    if len(shape1) == 2 and len(shape2) ==1 :
        out_shape = [arr1.shape[ii] for ii in range(arr1.ndim)
                     if ii in arr1.info['rows']]
        out_names = [arr1.info['axes'][ii] for ii in range(arr1.ndim)
                     if ii in arr1.info['rows']]
        out_vect = sp.empty(out_shape)
        out_vect = info_array(out_vect)
        out_vect = vect_array(out_vect, out_names)

        return out_vect
