"""GPU aceclerated Cholesky for clean map maker."""

import numpy as np
import scipy as sp
from constants import T_infinity

# Cython import
cimport numpy as np
cimport cython
# Doesn't seem to work, cdef extern instead.
#from libc.stdlib cimport malloc, free
np.import_array()    # Exposes the numpy c-api.

# We will do everything in double precision.
DTYPE = np.float
DTYPE_num = np.NPY_FLOAT64
ctypedef np.float_t DTYPE_t

cdef extern void cholesky_(int *n, double **A, int *lda)
cdef extern void test_cholesky_(int *n)
cdef extern void tinv_(int *n, double **A, int *lda)

cdef extern from "stdlib.h":
    void free(void* ptr)
    void* malloc(size_t size)

@cython.boundscheck(False)
@cython.wraparound(False)
def test_cholesky(n):

    cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] matrix
    cdef int n1 = n
    cdef int ii, jj

    matrix = large_empty((n, n))
    for ii in range(n):
        for jj in range(ii, n):
            matrix[ii, jj] = 0.5
        matrix[ii, ii] += ii + 1.
    import time
    st = time.time()
    cholesky_(&n1, <DTYPE_t **> matrix.data, &n1)
    print "Time was: ", time.time() - st

def test_cholesky_f(n):
    
    cdef int n1 = n
    test_cholesky_(&n1)

cdef class memory_container:
    """Creates and holds a memory buffer.
    
    Deallocates the memory when instance goes out of scope, even if other
    objects are using the buffer.
    """

    cdef void * buffer

    def __cinit__(self, size):
        self.buffer = <DTYPE_t *> malloc(size)
        if not buffer:
            raise MemoryError()

    cdef void * get_buffer(self):
        return self.buffer

    def __dealloc__(self):
        free(self.buffer)

class buffer_array(np.ndarray):
    _memory_handler = None

def large_empty(shape):
    
    # Get all the dimensions in the right format.
    cdef int nd = len(shape)
    # Very important to use this type.
    cdef np.npy_intp * dims = <np.npy_intp *> malloc(nd * sizeof(np.npy_intp))
    cdef long size = 1
    for ii, dim in enumerate(shape):
        dims[ii] = dim
        size *= dim
    # Allocate all the memory we need in a buffer.
    cdef memory_container mem
    mem = memory_container(size * sizeof(DTYPE_t))
    # Get the buffer and convert it to a numpy array.
    cdef void * buf = mem.get_buffer()
    arr = np.PyArray_SimpleNewFromData(nd, dims, DTYPE_num, buf)
    # Store the only living reference (once this function returns) to the
    # memory container on the array.
    arr = arr.view(buffer_array)
    arr._memory_handler = mem
    
    free(dims)
    return arr

@cython.boundscheck(False)
@cython.wraparound(False)
def call_cholesky(np.ndarray[DTYPE_t, ndim=2, mode='c'] matrix not None):
    """Specialized cholesky that uses gpu acceleration and minimal memory.
    
    The cholesky is performed inplace without ever accessing or writing to the
    lower triangle.
    """
    
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Input must be square")
    cdef int n=matrix.shape[0]
    cholesky_(&n, <DTYPE_t **> matrix.data, &n)

@cython.boundscheck(False)
@cython.wraparound(False)
def up_tri_copy(np.ndarray[DTYPE_t, ndim=2, mode='c'] origional not None):
    """Makes a copy of only the upper triagular part of a noise matrix.

    The lower triangle remains uninitialized. However, that memory should
    remain unallowcated on a virtual memory system.
    """
    
    if origional.shape[0] != origional.shape[1]:
        raise ValueError("Inputs must be square")
    # Counters.
    cdef int ii, jj
    # Allocate the memory for the output.
    cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] out
    out = large_empty((origional.shape[0], origional.shape[1]))
    for ii in xrange(origional.shape[0]):
        for jj in xrange(ii, origional.shape[1]):
            out[ii,jj] = origional[ii,jj]
    return out

@cython.boundscheck(False)
@cython.wraparound(False)
def inv_diag_from_chol(np.ndarray[DTYPE_t, ndim=2, mode='c'] chol not None,
                       np.ndarray[DTYPE_t, ndim=1, mode='c'] out not None):
    """From an upper triangular cholesky factor, find the diagonal of the
    inverse of the factored matrix.
    
    Uses back substitution but destroys the input cholesky factor. 
    """
    
    cdef int ii, jj, kk, n=chol.shape[0]
    cdef DTYPE_t tmp
    # Replace the cholesky factor by its own inverse.
    tinv_(&n, <DTYPE_t **> chol.data, &n)
    # Square the appropriate parts of the factor's inverse to get the inverse
    # diagonal of the factored matrix.
    for ii in range(n):
        out[ii] = 0
    for ii in range(n):
        for jj in range(ii, n):
            out[ii] += chol[ii,jj]**2

@cython.boundscheck(False)
@cython.wraparound(False)
def cho_solve(np.ndarray[DTYPE_t, ndim=2, mode='c'] chol not None,
               np.ndarray[DTYPE_t, ndim=1, mode='c'] b not None):
    """Use back and forward substitution to solve a system of equations.

    `chol` is the upper triagular cholesky factor for the system we want to
    solve.  This is to be used for map making and the physical limits defined
    in `map.constants` are used.
    """
    
    cdef long n = chol.shape[0]
    if chol.shape[1] != n or b.shape[0] != n:
        raise ValueError("Incompatible array dimensions")
    cdef long ii, jj
    cdef DTYPE_t tmp
    # Limit of small information, about 0.1 K**-1. Modes with less information
    # are ignored.
    cdef DTYPE_t small = 1./sp.sqrt(T_infinity)
    cdef np.ndarray[DTYPE_t, ndim=1, mode='c'] x
    x = b.copy()
    # Forward substitution.
    for ii in range(n):
        # Protect against very low information singular modes.
        if chol[ii,ii] < small:
            x[ii] = 0
            continue
        tmp = 0
        for jj in range(0, ii):
            tmp += x[jj] * chol[jj,ii]
        x[ii] = (x[ii] - tmp) / chol[ii,ii]
    # Back substitution.
    for ii in range(n - 1, -1, -1):
        # Protect against very low information singular modes.
        if chol[ii,ii] < small:
            x[ii] = 0
            continue
        tmp = 0
        for jj in range(ii + 1, n):
            tmp += x[jj] * chol[ii,jj]
        x[ii] = (x[ii] - tmp) / chol[ii,ii]
    return x
