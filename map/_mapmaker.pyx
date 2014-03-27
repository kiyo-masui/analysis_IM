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

def get_noise_inv_diag(
    np.ndarray[DTYPE_t, ndim=2, mode='c'] diagonal_inv not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] freq_modes not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] time_modes not None, 
    np.ndarray[DTYPE_t, ndim=4, mode='c'] freq_mode_update not None, 
    np.ndarray[DTYPE_t, ndim=4, mode='c'] time_mode_update not None, 
    np.ndarray[DTYPE_t, ndim=4, mode='c'] cross_update not None):

    # Shapes.
    cdef int n_chan = diagonal_inv.shape[0]
    cdef int n_time = diagonal_inv.shape[1]
    cdef int m = freq_modes.shape[0]
    cdef int q = time_modes.shape[0]
    # Indicies.
    cdef int t_ind, f_ind
    # Counters.
    cdef int ii, jj
    # Working floats.
    cdef DTYPE_t tmp
    # Output.
    cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] noise_inv_diag
    noise_inv_diag = sp.zeros((n_chan, n_time), dtype=DTYPE)
    # First we construct the update term.
    # Frequency modes.
    for ii in range(m):
        for jj in range(m):
            for f_ind in range(n_chan):
                tmp = freq_modes[jj,f_ind] * freq_modes[ii,f_ind]
                for t_ind in range(n_time):
                    noise_inv_diag[f_ind,t_ind] += (tmp 
                            * freq_mode_update[ii,t_ind,jj,t_ind])
    # Time modes.
    for ii in range(q):
        for jj in range(q):
            for t_ind in range(n_time):
                tmp = time_modes[ii,t_ind] * time_modes[jj,t_ind]
                for f_ind in range(n_chan):
                    noise_inv_diag[f_ind,t_ind] += (tmp 
                            * time_mode_update[ii,f_ind,jj,f_ind])
    # Cross terms.
    for ii in range(m):
        for jj in range(q):
            for t_ind in range(n_time):
                for f_ind in range(n_chan):
                    tmp = freq_modes[ii,f_ind] * time_modes[jj,t_ind]
                    noise_inv_diag[f_ind,t_ind] += 2. * (tmp 
                            * cross_update[ii,t_ind,jj,f_ind])
    # Multiply the update by the initial on both sides.
    noise_inv_diag *= diagonal_inv**2
    # Add the initial.
    noise_inv_diag = diagonal_inv - noise_inv_diag
    return noise_inv_diag


@cython.boundscheck(False)
@cython.wraparound(False)
def update_map_noise_independant_chan(
    np.ndarray[DTYPE_t, ndim=2, mode='c'] diagonal_inv not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] time_modes not None, 
    np.ndarray[DTYPE_t, ndim=3, mode='c'] time_mode_update not None, 
    np.ndarray[np.int_t, ndim=3, mode='c'] pointing_inds not None, 
    np.ndarray[DTYPE_t, ndim=2, mode='c'] pointing_weights not None, 
    f_ind_in,
    np.ndarray[DTYPE_t, ndim=4, mode='c'] map_noise_inv not None,
    ra0_ind_range=None,
    dec0_ind_range=None,
    ):
    """Convert noise to map space.
    
    This function is only used when ignoring frequency correlations.  The noise
    matrix is update one frequency slice at a time for performance.
    
    The total noise covariance matrix is 5D (freq, ra0, dec0, ra1, dec1).  This 
    function caculates at a since frequency, a range of r0, a range of dec1,
    all ra1 and all dec1.

    Note that ``map_noise_inv.shape[0]`` should equal 
    ``ra0_ind_range[1] - ra0_ind_range[0]`` and likewise for 
    ``map_noise_inv.shape[1]`` and the dec range.

    """

    # TODO: Change the order of the arguments to a more sane scheme, and make 
    # corresponding change in dirty_map.py.
    # XXX: Currently the arguments/function is backward compatible with before
    # the ra0, dec0 range subdivision.

    # Shapes.
    cdef int n_chan = diagonal_inv.shape[0]
    cdef int n_time = diagonal_inv.shape[1]
    cdef int n_pix_per_pointing = pointing_inds.shape[2]
    cdef int q = time_modes.shape[0]
    cdef int ra_ind_start, ra_ind_end
    cdef int dec_ind_start, dec_ind_end
    if not ra0_ind_range:
        ra_ind_start = 0
        ra_ind_end = map_noise_inv.shape[0]
    else:
        ra_ind_start = ra0_ind_range[0]
        ra_ind_end = ra0_ind_range[1]
    if not dec0_ind_range:
        dec_ind_start = 0
        dec_ind_end = map_noise_inv.shape[1]
    else:
        dec_ind_start = dec0_ind_range[0]
        dec_ind_end = dec0_ind_range[1]
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
            # Check if this time touches the ra and dec range we are
            # updating.
            for jj in xrange(n_pix_per_pointing):
                if (pointing_inds_transpose[time_ind,jj,0] >= ra_ind_start
                    and pointing_inds_transpose[time_ind,jj,0] < ra_ind_end
                    and pointing_inds_transpose[time_ind,jj,1] >= dec_ind_start
                    and pointing_inds_transpose[time_ind,jj,1] < dec_ind_end
                    ):
                    break
            else:
                continue
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
                if not (ra_ind >= ra_ind_start
                        and ra_ind < ra_ind_end
                        and dec_ind >= dec_ind_start
                        and dec_ind < dec_ind_end
                        ):
                    continue
                ra_ind = ra_ind - ra_ind_start
                dec_ind = dec_ind - dec_ind_start
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


