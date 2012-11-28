"""Polarized beam model."""

import math

import numpy as np
from scipy import special

class SimpleBeam(object):
    """Simple real space representation of the GBT beam.
    
    Represents the beam using a 2D hermite series in real space.  Only
    represents the I -> {XX, U, V, YY} leakage.
    """
    
    # Array index (reverse) priority: (freq, mode, mode, pol, x, y)

    def __init__(self, freq=1):
        freq = np.atleast_1d(freq)
        if freq.ndim != 1:
            raise ValueError("Frequency should be a 1D array.")
        self.freq = freq.copy()

    def set_width(self, width=1):
        """Set the overall gaussian width of the beam as a function of
        frequency."""
        width, tmp_f = np.broadcast_arrays(width, self.freq)
        if self.freq.shape != width.shape:
            raise ValueError("Width not broadcastable to same shape as"
                             " frequency.")
        self.width = width.copy()

    def set_coefficients(self, coefficients):
        """Set the coefficients of the all components of beam."""
        
        # Dimension checks.
        if coefficients.ndim != 4:
            msg = "Expected 4 dimensional data (freq, n, n, pol)"
            raise ValueError(msg)
        if coefficients.shape[0] != len(self.freq):
            msg = "Amplitude shapes incompatible number of channels."
            raise ValueError(msg)
        if coefficients.shape[3] != 4:
            msg = "Expected 4 polarizations."
            raise ValueError(msg)
        if coefficients.shape[2] != coefficients.shape[1]:
            msg = "Expected the same number of x and y modes."
            raise ValueError(msg)
        # Set the parameter.
        self.coefficients = coefficients.copy()

    def get_skewer(self, x, y):
        m = self.coefficients.shape[1]
        n_chan = len(self.freq)
        # Convert, x and y to 1D arrays broadcastable to the length of the
        # frequency axis.
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        if (x.shape != (1,) and x.shape != (n_chan,)
            or y.shape != (1,) and y.shape != (n_chan,)):
            msg = "Coordinate arrays must be scalars of 1D arrays the length"
            msg += " of the frequency axis."
            raise ValueError(msg)
        x.shape = x.shape + (1,)
        y.shape = y.shape + (1,)
        out = self.get_slice(x, y)
        out.shape = (n_chan, 4)
        return out
        
    def get_slice(self, x, y):
        m = self.coefficients.shape[1]
        n_chan = len(self.freq)
        # Convert, x and y to 2D arrays broadcastable to the length of the
        # frequency axis.
        x, y = np.broadcast_arrays(x, y)
        if x.ndim == 0:
            shape = (1, 1)
        elif x.ndim == 1:
            shape = (1, x.shape[0])
        elif x.ndim == 2:
            if x.shape[0] != n_chan and x.shape[0] != 1:
                msg = "Coordinate shapes not compatible with number of"
                msg += " channels."
                raise ValueError(msg)
            shape = (x.shape[0], x.shape[1])
        else:
            msg = "Coordinates must be broadcastable to (n_chan, n_coord)."
            raise ValueError(msg)
        x.shape = shape
        y.shape = shape
        n_points = x.shape[1]
        # Get the Gaussian factor.
        sigma = self.width / (2 * np.sqrt(2 * np.log(2)))
        sigma.shape = (n_chan, 1)
        gauss = np.exp(-(x**2 + y**2) / (2 * sigma**2))
        # Following line normalizes integral, but we want normalized amplitude.
        # gauss /= sigma * np.sqrt(2 * np.pi)
        # In the Hermite polynomials, the coordinates have to be scaled my the
        # width of the Gaussian **squared** (the weight function).
        window_sigma = sigma # Physicist Hermite polynomials.
        # Scale coordinates.
        x = x / window_sigma
        y = y / window_sigma
        # Allocate memory for output.
        out = np.zeros((n_chan, 4, n_points), dtype=np.float64)
        # Loop over the Hermite polynomials and add up the contributions.
        for ii in range(m):
            poly_x = special.eval_hermite(ii, x)[:,None,:]
            norm_x = 1. / math.sqrt(2**ii * math.factorial(ii))
            for jj in range(m - ii):
                poly_y = special.eval_hermite(jj, y)[:,None,:]
                norm = norm_x / math.sqrt(2**jj * math.factorial(jj))
                factor = norm * gauss[:,None,:]
                out += (self.coefficients[:,ii,jj,:,None]
                        * poly_x * poly_y * factor)
        return out


    def get_full(self, n_side, size):
        """Get full representation of the beam."""
        out = np.empty((len(self.freq), 4, n_side, n_side), dtype=np.float64)
        grid = (np.arange(n_side) - n_side/2.0 + 0.5) / n_side * size
        x, y = np.meshgrid(grid, grid)
        out = self.get_slice(x.flat, y.flat)
        out.shape = out.shape[:-1] + (n_side, n_side)
        return out


class Beam(object):
    """Class is a model for an approximately round single dish beam.
    
    Beam is modelled in the aperture plane.
    """
    
    def __init__(self, freq):
        # Set the observing frequencies.
        self.frequency = freq

    def set_aperture(self, D, R_p=1., R_c=1.):
        pass
    

    # Can be f dependant.
    def set_illumination(self, W, R_p=1., R_c=1.):
        pass


