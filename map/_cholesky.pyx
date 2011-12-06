"""GPU aceclerated Cholesky for clean map maker."""

import numpy as np
import scipy as sp

# Cython import
cimport numpy as np
cimport cython

# We will do everything in double precision.
DTYPE = np.float
ctypedef np.float_t DTYPE_t

cdef extern void cholesky_(int *n, double **A, int *lda)
cdef extern void tinv_(int *n, double **A, int *lda)

def test_cholesky(n):

    cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] matrix
    cdef int n1 = n

    matrix = sp.eye(n)
    matrix *= (sp.arange(n) + 1)**2
    matrix[4, 3] = 10
    print matrix
    import time
    st = time.clock()
    cholesky_(&n1, <DTYPE_t **> matrix.data, &n1)
    print matrix
    print "Time was: ", time.clock() - st

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
def up_tri_copy(np.ndarray[DTYPE_t, ndim=2, mode='c'] origional not None,
             np.ndarray[DTYPE_t, ndim=2, mode='c'] out not None):
    """Makes a copy of only the upper triagular part of a noise matrix.

    The lower triangle remains uninitialized. However, that memory should
    remain unallowcated on a virtual memory system.
    """

    cdef int ii, jj
    if (origional.shape[0] != out.shape[0] 
        or origional.shape[1] != out.shape[1]):
        raise ValueError("Inputs must have same dimensions.")
    if origional.shape[0] != origional.shape[1]:
        raise ValueError("Inputs must be square")
    for ii in xrange(origional.shape[0]):
        for jj in xrange(ii, origional.shape[1]):
            out[ii,jj] = origional[ii,jj]

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
