"""This module provides efficient storage for matricies that have block
diagonal structure."""

import scipy as sp
import numpy.lib.format as format

import kiyopy.custom_exceptions as ce

class FoldedMatrix(object) :
    """This class effeciently stores a matrix with block diagonal structure.

    The key idea of a Folded matrix is to store a 2D array (i.e. a matrix), that
    has some sort of block diagonal structure, as a 3D array so that the empty
    blocks do not occupy memory.  This class fully supports memory mapping such
    that it can be used for arrays that are to big to hold in memory.  The
    class is not meant to be a drop in replacement for numpy arrays.  When you
    want to use this as a numpy array, you will have to pull either all or part
    of the data out of this object.

    The external interface of this object looks somewhat like an (n, m) array
    where internally it is stored as (z, x, y) array with z*x=n and z*y=m. z is
    the number of blocks and (x, y) is the shape of each block.  Compatibility
    with these conditions is checked when creating a new instance.

    Usage :
    FoldedMatrix((n,m), (z,x,y))
        Create a new FoldedMatrix object, allowcating new memory.
    FoldedMatrix((n,m), (z,x,y), 'file_name')
        Same as above but memory map it to file_name (see open).
    FoldedMatrix((n,m), array)
        Create a FoldedMatrix object from a numpy array.  The array must have
        ndim=3 and have a shape compatible with (n,m).
    Mat = open('file_name'[, mode='r+']) 
        Open a saved memory mapped FoldedMatrix object. 
    """

    def __init__(self, array, file_name=None) :
        """See class doc string."""
        
        # If user passed an array, use that as storage.
        if isinstance(array, sp.ndarray) :
            if array.ndim == 3 :
                self.data = sp.array(array, copy=False)
            elif array.ndim == 2 :
                self.data = sp.array(array, copy=False).reshape((1,) +
                                                                array.shape)
        # Otherwise allowcate memory.
        elif len(array) == 3 or len(array) == 2 :
            if len(array) == 2 :
                array = (1,) + tuple(array)
            if file_name is None :
                self.data = sp.empty(array, dtype=float)
            elif isinstance(file_name, str) :
                self.data = format.open_memmap(file_name, mode='w+', 
                                               shape=array, dtype=float)
            else :
                raise TypeError("Argument 'file_name' must be a string.")
        else :
            raise ValueError("Second argument must be either a 2D or 3D array" 
                             "or a sequence of length 2 or 3.")
        # Set the outer dimensions: the shape of the expanded block diagonal
        # array.
        self.shape = (self.data.shape[0]*self.data.shape[1], 
                      self.data.shape[0]*self.data.shape[2])

