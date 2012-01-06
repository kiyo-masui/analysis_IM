"""Module holds fast optimized code for the map maker."""

import numpy as np
import scipy as sp

# Cython import
cimport numpy as np
cimport cython

# We will do everything in double precision.
DTYPE = np.float
ctypedef np.float_t DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
def update_map_noise_chan_ra_row(
    np.ndarray[DTYPE_t, ndim=2, mode='c'] diagonal_inv not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] freq_modes not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] time_modes not None, 
    np.ndarray[DTYPE_t, ndim=4, mode='c'] freq_mode_update not None, 
    np.ndarray[DTYPE_t, ndim=4, mode='c'] time_mode_update not None, 
    np.ndarray[DTYPE_t, ndim=4, mode='c'] cross_update not None, 
    np.ndarray[np.int_t, ndim=3, mode='c'] pointing_inds not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] pointing_weights not None, 
    f_ind_in,
    ra_ind_in,
    np.ndarray[DTYPE_t, ndim=4, mode='c'] map_noise_inv not None):
    """Convert noise to map space.

    We fill build the matrix one row at a time for performance reasons.
    This is implemented without using the algebra interface for the arrays,
    again for performance.
    """
    
    # Shapes.
    cdef int n_chan = diagonal_inv.shape[0]
    cdef int n_time = diagonal_inv.shape[1]
    cdef int n_pix_per_pointing = pointing_inds.shape[2]
    cdef int m = freq_modes.shape[0]
    cdef int q = time_modes.shape[0]
    # Indecies.
    cdef int time_ind
    cdef int f_ind = f_ind_in
    cdef int ra_ind = ra_ind_in
    cdef int dec_ind
    cdef int this_ra, this_dec
    # Counters.
    cdef int ii, jj, kk, pp, rr
    # working floats.
    cdef DTYPE_t tmp1, tmp2, weight
    # Allocate workspace.
    cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] update_term
    update_term = sp.empty((n_chan, n_time), dtype=DTYPE)
    # Need a reordered version of cross update and pointing_ind 
    # for efficient indexing.
    cdef np.ndarray[DTYPE_t, ndim=4, mode='c'] cross_transpose
    cross_transpose = sp.ascontiguousarray(
            cross_update.transpose((2, 3, 0, 1)))
    cdef np.ndarray[np.int_t, ndim=3, mode='c'] pointing_inds_transpose
    pointing_inds_transpose = sp.ascontiguousarray(
            pointing_inds.transpose((0, 2, 1)))
    # Now loop over the pointings.
    with nogil:
        for time_ind in xrange(n_time):
            # Check if this time touches the ra in we are updating.
            for jj in xrange(n_pix_per_pointing):
                if pointing_inds[time_ind,0,jj] == ra_ind:
                    break
            else:
                continue
            # Reset the time noise row to zero.
            for jj in xrange(n_chan):
                for kk in xrange(n_time):
                    update_term[jj,kk] = 0
            # Transform the four blocks of the update matrix to the right
            # space and add them.
            # The frequency modes
            for jj in xrange(m):
                tmp1 = freq_modes[jj,f_ind]
                for kk in xrange(m):
                    for pp in xrange(n_chan):
                        tmp2 = freq_modes[kk,pp] * tmp1
                        for rr in xrange(n_time):
                            update_term[pp,rr] -= \
                                    freq_mode_update[jj,time_ind,kk,rr] * tmp2
            # The time modes.
            for jj in xrange(q):
                tmp1 = time_modes[jj,time_ind]
                for kk in xrange(q):
                    for pp in xrange(n_chan):
                        tmp2 = time_mode_update[jj,f_ind,kk,pp] * tmp1
                        for rr in xrange(n_time):
                            update_term[pp,rr] -= time_modes[kk,rr] * tmp2
            # The cross blocks.
            for jj in xrange(q):
                tmp1 = time_modes[jj,time_ind]
                for kk in xrange(m):
                    for pp in xrange(n_chan):
                        tmp2 = freq_modes[kk,pp] * tmp1
                        for rr in xrange(n_time):
                            update_term[pp,rr] -= \
                                    cross_transpose[jj,f_ind,kk,rr] * tmp2
            for jj in xrange(m):
                tmp1 = freq_modes[jj,f_ind]
                for kk in xrange(q):
                    for pp in xrange(n_chan):
                        tmp2 = cross_update[jj,time_ind,kk,pp] * tmp1
                        for rr in xrange(n_time):
                            update_term[pp,rr] -= time_modes[kk,rr] * tmp2
            # Now weight by the diagonal.
            tmp1 = diagonal_inv[f_ind, time_ind]
            for jj in xrange(n_chan):
                for kk in xrange(n_time):
                    update_term[jj,kk] *= tmp1 * diagonal_inv[jj, kk]
            update_term[f_ind, time_ind] += diagonal_inv[f_ind, time_ind]
            # Change to the map domain and distribute to the output.
            for ii in xrange(n_pix_per_pointing):
                if pointing_inds[time_ind,0,ii] != ra_ind:
                    continue
                dec_ind = pointing_inds[time_ind,1,ii]
                weight = pointing_weights[time_ind,ii]
                # Loop over the time axes to convert to pixel and accumulate in
                # the output matrix.
                for jj in xrange(n_chan):
                    for kk in xrange(n_time):
                        tmp1 = update_term[jj,kk] * weight
                        for pp in xrange(n_pix_per_pointing):
                            this_ra = pointing_inds_transpose[kk,pp,0]
                            this_dec = pointing_inds_transpose[kk,pp,1]
                            tmp2 = tmp1 * pointing_weights[kk,pp]
                            map_noise_inv[dec_ind,jj,this_ra,this_dec] += tmp2

@cython.boundscheck(False)
@cython.wraparound(False)
def update_map_noise_independant_chan(
    np.ndarray[DTYPE_t, ndim=2, mode='c'] diagonal_inv not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] time_modes not None, 
    np.ndarray[DTYPE_t, ndim=3, mode='c'] time_mode_update not None, 
    np.ndarray[np.int_t, ndim=3, mode='c'] pointing_inds not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] pointing_weights not None, 
    f_ind_in,
    np.ndarray[DTYPE_t, ndim=4, mode='c'] map_noise_inv not None):
    """Convert noise to map space.
    
    This function is only used when ignoring frequency correlations.  The noise
    matrix is update one frequency slice at a time for performance.
    """
    
    # Shapes.
    cdef int n_chan = diagonal_inv.shape[0]
    cdef int n_time = diagonal_inv.shape[1]
    cdef int n_pix_per_pointing = pointing_inds.shape[2]
    cdef int q = time_modes.shape[0]
    # Indecies.
    cdef int time_ind
    cdef int f_ind = f_ind_in
    cdef int ra_ind
    cdef int dec_ind
    cdef int this_ra, this_dec
    # Counters.
    cdef int ii, jj, kk, pp, rr
    # working floats.
    cdef DTYPE_t tmp1, tmp2, weight
    # Allocate workspace.
    cdef np.ndarray[DTYPE_t, ndim=1, mode='c'] update_term
    update_term = sp.empty(n_time, dtype=DTYPE)
    # Need a reordered version of pointing_inds 
    # for efficient indexing.
    cdef np.ndarray[np.int_t, ndim=3, mode='c'] pointing_inds_transpose
    pointing_inds_transpose = sp.ascontiguousarray(
            pointing_inds.transpose((0, 2, 1)))
    # Now loop over the pointings.
    with nogil:
        for time_ind in xrange(n_time):
            # Reset the time noise row to zero.
            for jj in xrange(n_time):
                update_term[jj] = 0
            # The time modes.
            for jj in xrange(q):
                tmp1 = time_modes[jj,time_ind]
                for kk in xrange(q):
                    tmp2 = time_mode_update[f_ind,jj,kk] * tmp1
                    for rr in xrange(n_time):
                        update_term[rr] -= time_modes[kk,rr] * tmp2
            # Now weight by the diagonal.
            tmp1 = diagonal_inv[f_ind, time_ind]
            for jj in xrange(n_time):
                update_term[jj] *= tmp1 * diagonal_inv[f_ind,jj]
            update_term[time_ind] += diagonal_inv[f_ind, time_ind]
            # Change to the map domain and distribute to the output.
            for ii in xrange(n_pix_per_pointing):
                ra_ind = pointing_inds[time_ind,0,ii]
                dec_ind = pointing_inds[time_ind,1,ii]
                weight = pointing_weights[time_ind,ii]
                # Loop over the time axes to convert to pixel and accumulate in
                # the output matrix.
                for kk in xrange(n_time):
                    tmp1 = update_term[kk] * weight
                    for pp in xrange(n_pix_per_pointing):
                        this_ra = pointing_inds_transpose[kk,pp,0]
                        this_dec = pointing_inds_transpose[kk,pp,1]
                        tmp2 = tmp1 * pointing_weights[kk,pp]
                        map_noise_inv[ra_ind,dec_ind,this_ra,this_dec] += tmp2


