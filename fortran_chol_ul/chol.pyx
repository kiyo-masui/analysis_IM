from numpy cimport ndarray
import numpy as np

from libc.stdint cimport int64_t

cdef extern:
    void c_cholesky(int64_t *n, double *A, int64_t *lda)

def p_cholesky(ndarray[double, mode='c', ndim=2] A):
    cdef int64_t n, lda
    n = A.shape[1]
    lda = A.shape[0]
    c_cholesky(&n, &A[0,0], &lda)


def test(n):
    A = np.eye(n)
    A *= np.arange(1, n+1)**2
    p_cholesky(A)
    print A.flat[::n+1]

