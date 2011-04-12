"""Module contains classes representing the beam operator."""

import scipy as sp
from scipy import interpolate
from scipy.ndimage.filters import convolve

from core import algebra
import kiyopy.custom_exceptions as ce


class Beam(object) :
    """Object representing the beam operator.
    
    This is appropriate for all beams that are round and stationary (but may
    be frequency dependant).  Specific beam forms should inherit from this one
    whethar they be imperical or functional.
    """
    
    def apply(self, alg_object, wrap=False, right_apply=False) :
        """Apply the beam, as a linear operator, to vector or matrix.

        If the beam is viewed as a matrix, this operation is equivalent to
        matrix multiplication with the beam on the left.
        """

        if ((not 'freq' in alg_object.axes) 
            or (not 'ra' in alg_object.axes)
            or (not 'dec' in alg_object.axes)) :
                raise ce.DataError("Beam operation only works in frequncy, "
                                   "ra, dec, coords.")
        
        # allowcate memory for the output.
        out = algebra.zeros_like(alg_object)
        # Figure out the pixel sizes (in real degrees).
        dfreq = alg_object.info['freq_delta']
        dra = alg_object.info['ra_delta']
        dra /= sp.cos(alg_object.info['dec_centre']*sp.pi/180)
        ddec = alg_object.info['dec_delta']
        # Figure out the convolution mode.
        if wrap :
            mode = 'wrap'
        else :
            mode = 'constant'
        # Loop over frequencies and do convolution one frequency at a time.
        freq_array = alg_object.get_axis('freq')
        for ii, freq in enumerate(freq_array) :
            # How wide the kernal has to be.
            width = self.kernal_size(freq)
            # Make sure the dimensions are an odd number of pixels.
            nkx = width//dra
            if nkx%2 == 0 :
                nkx += 1
            nky = width//ddec
            if nky%2 == 0 :
                nky += 1
            # Calculate kernal lags.
            lagsx = (sp.arange(nkx, dtype=float) - (nkx - 1)//2)*dra
            lagsy = (sp.arange(nky, dtype=float) - (nky - 1)//2)*ddec
            lags_sq = lagsx[:, None]**2 + lagsy[None, :]**2
            # Make gaussian beam profile.
            kernal = dra*ddec*self.beam_function(lags_sq, freq, 
                                                 squared_delta=True)
            if isinstance(alg_object, algebra.vect) :
                if alg_object.axes != ('freq', 'ra', 'dec') :
                    raise ce.DataError("Vector axis names must be exactly "
                                       "('freq', 'ra', 'dec')")
                # Do the convolution.
                convolve(alg_object[ii, ...], kernal, out[ii], mode=mode,
                         cval=0)
            elif isinstance(alg_object, algebra.mat) :
                # If applying from the left, loop over columns and convolve
                # over rows.  If applying from the right, do the oposite.
                if right_apply :
                    if alg_object.col_names() != ('freq', 'ra', 'dec') :
                        raise ce.DataError("Matrix column axis names must be "
                                           "exactly ('freq', 'ra', 'dec')")
                    iterator = alg_object.iter_row_index()
                else :
                    if alg_object.row_names() != ('freq', 'ra', 'dec') :
                        raise ce.DataError("Matrix row axis names must be "
                                           "exactly ('freq', 'ra', 'dec')")
                    iterator = alg_object.iter_col_index()
                for index in iterator :
                    sub_mat = alg_object[index] # A view.
                    # Pick out this frequency.
                    sub_mat = sub_mat[ii, ...]
                    # make a view of the ouput array.
                    sub_out = out[index]
                    sub_out = sub_out[ii, ...]
                    # Do the convolution.
                    convolve(sub_mat, kernal, sub_out, mode=mode,
                             cval=0)
        return out

class GaussianBeam(Beam) :
    """Object representing the beam operator for the special case that it is an
    exact frequency dependant Gaussian.
    """

    def __init__(self, width, freq=None) :
        
        # Calculate the standard deviation.
        sig = width/(2*sp.sqrt(2*sp.log(2)))

        # Make function that returns the beam sigma.
        if freq is None :
            # Frequency independant case.
            def sigma(freq) :
                return float(sig)
            self.sigma = sigma
        else :
            self.sigma = interpolate.interp1d(freq, sig)

    def beam_function(self, delta_r, frequency, squared_delta=False) :
        """Returns the beam weight as a funciton of angular lag.
        
        This has units 1/strad.
        """
        
        # Get the width.
        sig = self.sigma(frequency)
        # Square the lags.
        if squared_delta :
            dr_sqr = delta_r
        else :
            dr_sqr = delta_r**2
        # Calculate the profile.
        profile = sp.exp(-dr_sqr/(2*sig**2))
        profile *= 1/(2*sp.pi*sig**2)

        return profile

    def kernal_size(self, frequency) :
        """Gets the minimum convolution kernal size in degrees."""

        # For gaussian, 5 sigma taper is good enough.
        return 10.0*self.sigma(frequency)

