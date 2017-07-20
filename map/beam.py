"""Module contains classes representing the beam operator."""
import scipy as sp
import numpy as np
from scipy import interpolate
from scipy.ndimage.filters import convolve

from core import algebra
import kiyopy.custom_exceptions as ce


class Beam(object):
    """Object representing the beam operator.

    This is appropriate for all beams that are round and stationary (but may
    be frequency dependent).  Specific beam forms should inherit from this one
    whether they be empirical or functional.
    """

    def apply(self, alg_ob, mode="constant", cval=0, right_apply=False):
        """Apply the beam, as a linear operator, to vector or matrix.

        This operation is equivalent to matrix multiplication by the beam
        matrix.  The matrix multiplication can be performed with the beam
        matrix on either the left or the right.  The beam matrix is a symmetric
        matrix.

        Parameters
        ----------
        alg_ob: `vect` or `mat` subclass.
            Object to which the beam will be applied.  This object must live in
            map space.  If it is a vect, alg_ob.axes must be ('freq', 'ra',
            'dec').  If it is a mat, and `right_apply` is False, then
            alg_ob.row_names() must return ('freq', 'ra', 'dec').  If
            `right_apply` is True, then alg_ob.col_names() should return
            this tuple.  Also, the meta data for these three axis must be set
            (see `algebra.alg_object.set_axis_info`).
        wrap: bool
            If tests True, use periodic boundary conditions for the
            convolution. Otherwise zero pad (default).
        right_apply: bool
            Whether to apply the beam operator with the from the left (False,
            default) or from the right (True).  If `alg_ob` is a vect subclass,
            this has no effect (because the beam matrix is symmetric).

        Returns
        -------
        out: `vect` or `mat` subclass same shape as `alg_ob`.
            Convolved map or matrix.

        Notes
        -----
        The missing feature here is the ability to preallocate memory or to be
        able to overwrite the input alg_ob in place.  This will be important
        for large matrices that we can only hold in memory one at a time.
        Also this would be pretty easy to implement.
        """
        if ((not 'freq' in alg_ob.axes)
            or (not 'ra' in alg_ob.axes)
            or (not 'dec' in alg_ob.axes)):
            raise ce.DataError("Beam operation only works in frequency, "
                               "ra, dec, coords.")

        out = algebra.zeros_like(alg_ob)

        # Figure out the pixel sizes (in real degrees).
        dra = abs(alg_ob.info['ra_delta'])
        dra /= sp.cos(alg_ob.info['dec_centre'] * sp.pi / 180.)
        ddec = abs(alg_ob.info['dec_delta'])

        # Loop over frequencies and do convolution one frequency at a time.
        freq_array = alg_ob.get_axis('freq')
        for ii, freq in enumerate(freq_array):
            width = self.kernel_size(freq)

            # Make sure the dimensions are an odd number of pixels.
            nkx = width // abs(dra)
            if nkx % 2 == 0:
                nkx += 1

            nky = width // abs(ddec)
            if nky % 2 == 0:
                nky += 1

            # Calculate kernel lags.
            lagsx = (sp.arange(nkx, dtype=float) - (nkx - 1) // 2) * dra
            lagsy = (sp.arange(nky, dtype=float) - (nky - 1) // 2) * ddec
            lags_sq = lagsx[:, None] ** 2. + lagsy[None, :] ** 2.

            kernel = dra * ddec * self.beam_function(lags_sq, freq,
                                                     squared_delta=True)
            kernel /= np.sum(kernel)

            if isinstance(alg_ob, algebra.vect):
                if alg_ob.axes != ('freq', 'ra', 'dec'):
                    raise ce.DataError("Vector axis names must be exactly "
                                       "('freq', 'ra', 'dec')")

                convolve(alg_ob[ii, ...], kernel, output=out[ii], mode=mode,
                         cval=cval)

            elif isinstance(alg_ob, algebra.mat):
                # If applying from the left, loop over columns and convolve
                # over rows.  If applying from the right, do the opposite.
                if right_apply:
                    if alg_ob.col_names() != ('freq', 'ra', 'dec'):
                        raise ce.DataError("Matrix column axis names must be "
                                           "exactly ('freq', 'ra', 'dec')")

                    iterator = alg_ob.iter_row_index()

                else:
                    if alg_ob.row_names() != ('freq', 'ra', 'dec'):
                        raise ce.DataError("Matrix row axis names must be "
                                           "exactly ('freq', 'ra', 'dec')")

                    iterator = alg_ob.iter_col_index()

                for index in iterator:
                    sub_mat = alg_ob[index]
                    # Pick out this frequency.
                    sub_mat = sub_mat[ii, ...]
                    # make a view of the output array.
                    sub_out = out[index]
                    sub_out = sub_out[ii, ...]

                    convolve(sub_mat, kernel, sub_out, mode=mode, cval=cval)
        return out

    def angular_transform(self, frequency):
        """Return angular beam Fourier transform function."""

        raise NotImplementedError("General function not written yet.")

    def radial_transform(self, width):
        """Return the radial beam Fourier transform function.

        In the radial direction the beam is just top hat function, so the
        Fourier transform is a sinc function.

        Parameters
        ----------
        width: float
            The frequency width of the beam function.

        Returns
        -------
        transform: function
            Call signature: transform(k_rad). Vectorized radial beam transform
            as a function of radial wave number.  Accepts an array of wave
            numbers (with units the reciprocal of those of `width`) and returns
            an array of the same shape.
        """

        factor = width / 2. / sp.pi
        return lambda k_rad: sp.sinc(k_rad * factor)

    def radial_real_space_window(self, width1, width2, return_limits=False):
        """Return the radial beam window function in real space.

        Gives the convolution of the radial parts of the beams for two pixels
        along the frequency axis. Since the beam is just a top hat in the
        radial direction, the window is just a trapezoidal function.

        Parameters
        ----------
        width1, width2: floats
            The frequency widths of the beam functions for the two pixels.
            Note that this information is not stored in this class (which only
            hold angular information about the beam).
        return_limits: bool
            Also return the integration limits that will capture most of the
            of the window. Default is False.

        Returns
        -------
        window: function
            Vectorized window function, accepting and returning numpy arrays
            of floats.  The function is centered at 0, not delta_f.
        limits: tuple
            Only returned if `return_limits` is True.  Limits of integration
            that will capture most of the window.
        """

        if width1 <= 0 or width2 <= 0:
            raise ValueError("Widths must be floats greater than 0.")

        if width2 < width1:
            tmp = width1
            width1 = width2
            width2 = tmp

        def window(freq):
            # Outside the window case.
            out = sp.empty(freq.shape, dtype=float)
            out[abs(freq) > (width1 + width2) / 2.] = 0.

            # Flat part of the function.
            plateau = 1. / width2
            out[abs(freq) < (width2 - width1) / 2.] = plateau

            # Climbing region.
            mask = ((freq <= (width1 - width2) / 2.) &
                    (freq >= -(width1 + width2) / 2.))
            out[mask] = (plateau * \
                        ((width2 + width1) / 2. + freq[mask]) / width1)

            # Falling region.
            mask = ((freq >= (width2 - width1) / 2.) &
                    (freq <= (width1 + width2) / 2.))

            out[mask] = (plateau * \
                        ((width2 + width1) / 2. - freq[mask]) / width1)

            return out

        if return_limits:
            return window, (-(width1 + width2) * 0.6, (width1 + width2) * 0.6)
        else:
            return window


class GaussianBeam(Beam):
    """Beam operator for Gaussian angular shape.

    This object inherits from the `Beam` object and gets much of its
    functionality from there.  The object contains all the information about
    the beam of an instrument (which may be frequency dependent) and provides
    operations like the convolution of a map with the beam and matrix
    operations with the beam.

    This function fixes only the angular parts of the beam.  Along the
    frequency axis, the beam is assumed to be a top hat.  The frequency width
    of the beam is not fixed.  The width of the objects that the beam works
    with is assumed.  This reflects the fact that 21 cm experiments typically
    have essentially perfect frequency resolution and therefore the beam along
    the frequency axis is just a matter of binning.

    Parameters
    ----------
    width: float or array
        The full width half max of the Gaussian beam.  If an array is provided,
        the beam is frequency dependent, with the corresponding frequencies
        given in `freq`. Linear interpolation is used to fill out the function.
        Generally any units should work as long as you are consistent.
    freq: None or array
        The frequencies that correspond to the the provided widths.  A
        frequency independent beam is given by freq=None and only a single
        number for `width`.  Any units will work as long as you are consistent.
    extrapolate: bool
        whether to allow extrapolation of the beam function beyond the provided
        frequencies. Default is False.
    """

    def __init__(self, width, freq=None, extrapolate=False):
        # Calculate the standard deviation.
        sig = width / (2. * sp.sqrt(2. * sp.log(2.)))

        # Make function that returns the beam sigma.
        # first: frequency independent case.
        if freq is None:
            sig = float(sig)

            def sigma(f):
                return sig

            self._sigma = sigma
        else:
            if not extrapolate:
                self._sigma = interpolate.interp1d(freq, sig)
            else:
                self._sigma = interpolate.InterpolatedUnivariateSpline(
                    freq, sig, k=1)

    def beam_function(self, delta_r, frequency, squared_delta=False):
        """Return the beam weight as a function of angular lag.

        This gives the angular beam as a function of angular lag.  This
        function is fully vectorized. The result has units 1/degrees**2.

        Parameters
        ----------
        delta_r: array
            Angular lags at which to evaluate the beam function.  Units of
            degrees.
        frequency: float or array
            Frequencies at which to evaluate the beam.  Either a single
            frequency can be provided or one for each angular lag.
        squared_delta: bool
            If True, assume that `delta_r` is the squared angular lag instead
            of just the angular lag.  Default is False.

        Returns
        -------
        profile: array
            The beam profile value at the given lags.
        """
        sig = self._sigma(frequency)

        if squared_delta:
            dr_sqr = delta_r
        else:
            dr_sqr = delta_r ** 2.

        # Calculate the profile.
        profile = sp.exp(-dr_sqr / (2. * sig ** 2.))
        profile *= 1. / (2. * sp.pi * sig ** 2.)

        return profile

    def angular_transform(self, frequency):
        """Return the angular beam Fourier transform function.

        Parameters
        ----------
        frequency: float
            For a frequency dependent beam, the frequency at which the Fourier
            transform should be calculated.

        Returns
        -------
        transform: function
            Call signature: transform(k_trans). Vectorized angular beam
            transform function. Accepts an array of wave numbers (generally
            with unites inverse degrees) and returns an array of the same
            shape.
        """
        frequency = float(frequency)
        sig = self._sigma(frequency)
        factor = sig ** 2. / 2.

        return lambda k_trans: sp.exp(-factor * k_trans ** 2.)

    def angular_real_space_window(self, f1, f2, return_limits=False):
        """Return a vectorized angular real space window function.

        Gives the convolution of the angular parts of the beams for two pixels
        along the frequency axis. Since the beam is just a Gaussian in the
        radial direction, the window is also a Gaussian.

        Parameters
        ----------
        f1, f2: floats
            Center frequencies of the two pixels for which to construct the
            window.
        return_limits: bool
            Also return the integration limits that will capture most of the
            of the window. Default is False.

        Returns
        -------
        Window: function
            Vectorized window function.  Maps angular lag (usually in degrees)
            onto the window function (usually inverse degrees).  Both inputs
            and outputs should be arrays.  Window function is centered at zero,
            not the angular separation of the pixels.
        limits: tuple
            Only returned if `return_limits` is True.  Radial limits of
            integration that will capture most of the window.
        """
        v1 = self._sigma(f1) ** 2.
        v2 = self._sigma(f2) ** 2.
        v = v1 + v2
        norm = 1. / 2. / sp.pi / v
        factor = 1. / 2. / v

        if return_limits:
            return lambda lag: norm * sp.exp(-lag ** 2. * factor), \
                               (0, 5. * sp.sqrt(v))

        else:
            return lambda lag: norm * sp.exp(-lag ** 2. * factor)

    def kernel_size(self, frequency, taper=4.):
        """Return the minimum convolution kernel size in degrees.

        Gives the width of the box that will adequately represent the whole
        beam.

        Parameters
        ----------
        frequency: array
            Frequencies at which to get the kernel size.

        Returns
        -------
        widths: array
            The width of the required box at each frequency.  In degrees if
            `width` was given in degrees on initialization of this class.
        """
        return 2. * taper * self._sigma(frequency)
