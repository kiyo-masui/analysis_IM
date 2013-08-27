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
            shape = x.shape
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


class HermiteBasis(object):
    """Hermite basis fucntions for beam maps.
    
    Parameters
    ----------
    freq : 1D array of length n
        Frequency corresponding to each channel.
    center : n x 2 array
        Beam center offsets.
    width : 1D array of length n
        Beam width.
    """

    def __init__(self, freq, center=None, width=None):
        freq = np.atleast_1d(freq)
        if freq.ndim != 1:
            msg = "Expected `freq` to be a 1D array."
            raise ValueError(msg)
        n_chan = len(freq)
        if center is None:
            center = [0, 0]
        if width is None:
            width = [1]
        # Broadcast the center and width up to the right shape.
        width, tmp_f = np.broadcast_arrays(width, freq)
        if width.shape != (n_chan,):
            raise ValueError("Width not broadcastable to same shape as"
                             " frequency.")
        center, tmp_f = np.broadcast_arrays(center, freq[:,None])
        if center.shape != (n_chan, 2):
            raise ValueError("Center offsets not broadcastable to same "
                             "length as frequency.")
        self.freq = np.array(freq, copy=True, dtype=np.float64)
        self.width = np.array(width, copy=True, dtype=np.float64)
        self.center = np.array(center, copy=True, dtype=np.float64)

    def eval_basis(self, mode, x, y):
        """Evaluate a basis function at all frequencies.
        
        Parameters
        ----------
        mode : (int, int)
            What Hermite polynomial to evalutate.
        x, y : arrays
            Coordinate locations to evalutate the polynomial. If 1D array,
            evaluate the polynomial at the same location at for each channel.
            If 2D array, evaluate at a different location for each channel.
            First dimension much be compatible with number of channels. 
        """

        n_chan = len(self.freq)
        # Convert, x and y to 2D arrays broadcastable to the length of the
        # frequency axis.
        x, y = np.broadcast_arrays(x, y)
        return_1D = False
        if x.ndim == 0:
            shape = (1, 1)
            return_1D = True
        elif x.ndim == 1:
            shape = (1, x.shape[0])
        elif x.ndim == 2:
            if x.shape[0] != n_chan and x.shape[0] != 1:
                msg = "Coordinate shapes not compatible with number of"
                msg += " channels."
                raise ValueError(msg)
            shape = x.shape
        else:
            msg = "Coordinates must be broadcastable to (n_chan, n_coord)."
            raise ValueError(msg)
        x.shape = shape
        y.shape = shape
        n_points = x.shape[1]
        # Allocate memory for output.
        out = np.zeros((n_chan, n_points), dtype=np.float64)
        # Loop over channels and evaluate basis function in each.
        for ii in range(n_chan):
            if x.shape[0] == n_chan:
                this_x = x[ii,:]
                this_y = y[ii,:]
            else:
                this_x = x[0,:]
                this_y = y[0,:]
            out[ii,:] = self.eval_basis_chan(mode, this_x, this_y, ii)
        if return_1D:
            out.shape = (n_chan,)
        return out

    def eval_basis_chan(self, mode, x, y, chan):
        """Evaluate a basis function at a specific frequency channel.
        
        Parameters
        ----------
        mode : (int, int)
            What Hermite polynomial to evalutate.
        x, y : arrays
            Coordinate locations to evalutate the polynomial. If 1D array,
            evaluate the polynomial at the same location at for each channel.
        chan : int
            Index of the frequency channel.
        """

        x, y = np.broadcast_arrays(x, y)
        # Include the center offsets.
        x = x - self.center[chan,0]
        y = y - self.center[chan,1]
        # Get Gaussian factor.
        sigma = self.width[chan] / (2 * np.sqrt(2 * np.log(2)))
        gauss = np.exp(-(x**2 + y**2) / (2 * sigma**2))
        # In the Hermite polynomials, the coordinates have to be scaled by the
        # width of the Gaussian **squared** (the weight function).
        window_sigma = sigma # Physicist Hermite polynomials.
        x /= window_sigma
        y /= window_sigma
        # Construct the basis function.
        poly_x = special.eval_hermite(mode[0], x)
        norm = 1. / math.sqrt(2**mode[0] * math.factorial(mode[0]))
        poly_y = special.eval_hermite(mode[1], y)
        norm = norm / math.sqrt(2**mode[1] * math.factorial(mode[1]))
        out = poly_x * poly_y * norm * gauss
        return out
    
    def get_basis_grid(self, mode, n_side, size):
        """Get a full grid representation of the basis function."""
        
        # Grid points.
        grid = (np.arange(n_side) - n_side/2.0 + 0.5) / n_side * size
        x, y = np.meshgrid(grid, grid)
        # Flatten.
        x.shape = y.shape = (n_side**2)
        out = self.eval_basis(mode, x, y)
        # Reshape into a grid.
        out.shape = (len(self.freq), n_side, n_side)
        return out

class LinearBeam(object):
    """Beam object built from a set of basis functions and thier coefficients.
    """
    
    def __init__(self, Basis, coefficients):
        
        self.n_chan = coefficients.shape[0]
        self.n_pol = coefficients.shape[1]
        # Some checks to make sure that the basis and the coefficients are
        # compatible.
        # First check if we have one basis or a list of basis objects (one per
        # pol).
        #try:
        #    if len(Basis) != self.n_pol:
        #        raise RuntimeError()
        #except ValueError:
        #    Basis = [Basis] * self.n_pol
        #except RuntimeError:
        #    msg = ("Must either pass one Basis Object or a list of them "
        #           "with length n_pol.")
        #    raise ValueError(msg)
        # Make sure all Basis objects have the right length frequency axis.
        #for B in Basis:
        #    if len(B.freq) != self.n_chan:
        #        msg = ("Basis objects mush have frequency axis length matching"
        #               " coefficients.")
        #        raise ValueError(msg)
        if len(Basis.freq) != self.n_chan:
            raise ValueError("Basis object mush have frequency axis length"
                             " matching coefficients.")
        # Now check the coefficients.
        self.n_modes = coefficients.shape[2]
        if coefficients.shape != (self.n_chan, self.n_pol,
                                  self.n_modes, self.n_modes):
            raise ValueError("Beam_coefficients not the right shape.")
        
        # Store parameters.
        self.Basis = Basis
        self.coefficients = coefficients

    def get_skewer(self, x, y):
        
        out = np.zeros((self.n_chan, self.n_pol), dtype=np.float64)
        for ii in range(self.n_modes):
            for jj in range(self.n_modes):
                if not np.all(self.coefficients[:,:,ii,jj] == 0):
                    basis_function = self.Basis.eval_basis((ii, jj), x, y)
                    out += (basis_function[:,None]
                            * self.coefficients[:,:,ii,jj])
        return out

    def get_full_beam(self, n_side, size):

        out = np.zeros((self.n_chan, self.n_pol, n_side, n_side),
                       dtype=np.float64)
        for ii in range(self.n_modes):
            for jj in range(self.n_modes):
                if not np.all(self.coefficients[:,:,ii,jj] == 0):
                    basis_function = self.Basis.get_basis_grid((ii, jj), 
                                                               n_side, size)
                    out += (basis_function[:,None,:,:]
                            * self.coefficients[:,:,ii,jj,None,None])
        return out
                    

def plot_beam_map(beam_map, color_map=1, side=None, normalize=None,
                  rotate=None):
    """
    Plotting function for precomputed beam maps.
    
    Parameters
    ----------
    beam_map : array with shape [n_pol, n_side, n_side]
        Map of the beam at a single frequency but all polarizations.
    color_map : scalar
        power law index for the color scale. 1 for linear color, .5 for sqrt,
        etc.
    side : float
        size of the side of the map, e.g. 1 degree
    normalize : list or string
        How to normalize the polarizations.  Either pass a list of scalars (one
        for each polarization) or a key such as 'max03' which normalized to the
        maximum of the 0th and 3rd polarization (appropriate for XX,XY,YX,YY
        polarizations).
    rotate : string
        How to rotate the polarizations (after normalizing).  Currently only
        'XXYYtoIQ' is implemented.
    """

    import matplotlib.pyplot as plt

    n_pol = beam_map.shape[0]
    n_side = beam_map.shape[1]
    # Assume square.
    if beam_map.shape[2] != n_side:
        raise ValueError("maps should be square.")
    if side is None:
        side = n_side
    side = float(side)
    beam_map = beam_map.copy()

    if normalize:
        if normalize == 'max03':
            if n_pol != 4:
                raise ValueError("`normalize='max03'` requires 4 polarization"
                                 " beam.")
            XX_amp = np.amax(beam_map[0])
            YY_amp = np.amax(beam_map[3])
            XY_amp = np.sqrt(XX_amp * YY_amp)
            normalize = [XX_amp, XY_amp, XY_amp, YY_amp]
        elif len(normalize) == n_pol:
            pass
        elif len(normalize) == 1:
            normalize = normalize * n_pol
        else:
            raise ValueError("`normalize` must be a list of scalar factors"
                             " with length `n_pol` or 1; or the option"
                             " 'max03'.")
        for ii in range(n_pol):
            beam_map[ii] /= normalize[ii]
    if rotate:
        if rotate == "XXYYtoIQ":
            if n_pol != 4:
                raise ValueError("`rotate='XXYYtoIQ'` requires 4 polarization"
                                 " beam.")
            rotate = np.zeros((4, 4))
            rotate[0,0] = 0.5
            rotate[0,3] = 0.5
            rotate[1,0] = 0.5
            rotate[1,3] = -0.5
            rotate[2,1] = 1
            rotate[3,2] = 1
        else:
            msg = ("`rotate` must be 'XXYYtoIQ' or a 2D array with second"
                   " dimension of length `n_pol`")
            try:
                if rotate.ndim != 2:
                    raise ValueError()
                if rotate.shape[1] != n_pol:
                    raise ValueError()
            except ValueError:
                raise ValueError(msg)
            except AttributeError:
                raise TypeError(msg)
        broadcast = ((slice(None), slice(None))
                     + (None,) * (beam_map.ndim - 1))
        beam_map = np.sum(rotate[broadcast] * beam_map, axis=1)

    # Put though the colormap.
    neg = beam_map < 0
    beam_map = abs(beam_map) ** color_map
    beam_map[neg] *= -1

    beam_map.shape = (n_pol * n_side, n_side)
    plt.figure()
    plt.imshow(beam_map.T, extent=[0, n_pol * side, 0, side])
    #plt.colorbar()

def compare_beam_maps(map1, map2, color_map=1):

    import matplotlib.pyplot as plt

    if map1.shape != map2.shape:
        raise ValueError("Map shapes must agree.")
    n_pol = map1.shape[0]
    n_side = map1.shape[1]
    # Assume square.
    if map1.shape[2] != n_side:
        raise ValueError("maps should be square.")
    # Normalize each map.  Assume polarizations are XX, XY, YX, YY.
    # Map 1.
    norm_x = map1[0,n_side//2,n_side//2]
    norm_y = map1[3,n_side//2,n_side//2]
    norm_cross = np.sqrt(norm_x * norm_y)
    map1 = map1.copy()
    map1[0] /= norm_x
    map1[3] /= norm_y
    map1[[1,2]] /= norm_cross
    # Map 2.
    norm_x = map2[0,n_side//2,n_side//2]
    norm_y = map2[3,n_side//2,n_side//2]
    norm_cross = np.sqrt(norm_x * norm_y)
    map2 = map2.copy()
    map2[0] /= norm_x
    map2[3] /= norm_y
    map2[[1,2]] /= norm_cross
    # We want to plot the difference.
    beam_map = map1 - map2
    # Put though the colormap.
    neg = beam_map < 0
    beam_map = abs(beam_map) ** color_map
    beam_map[neg] *= -1
    # Plot.
    beam_map.shape = (n_pol * n_side, n_side)
    plt.figure()
    plt.imshow(beam_map.T)
    plt.colorbar()


