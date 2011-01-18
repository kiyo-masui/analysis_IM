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

    Note that the structure of the matrix does not have to be restricted to the
    primary diagonal. More on this idea later.
    """

    def __init__(self, mat_dims, array, file_name=None) :
        """See class doc string."""
        
        if len(mat_dims) == 2 :
            self.shape = tuple(mat_dims)
        else :
            raise ValueError("Argument mat_dims must be sequence of length 2")
        if isinstance(array, sp.ndarray) and array.ndim==3 :
            self.data = array
        elif len(array) == 3 :
            if file_name is None :
                self.data = sp.empty(array, dtype=float)
            elif isinstance(file_name, str) :
                self.data = format.open_memmap(file_name, mode='w+', 
                                               shape=array, dtype=float)
            else :
                raise TypeError("Argument 'file_name' must be a string.")
        else :
            raise ValueError("Second argument must be either of 3D array or a "
                             "sequence of length 3.")
        
