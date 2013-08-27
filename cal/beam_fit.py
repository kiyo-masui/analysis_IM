"""Module for fitting a beam model to point source scans."""

import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
from scipy import optimize, linalg

from utils import misc
import source



class FormattedData(object):
    """All data that goes into a fit.

    This class holds GBT SDfits data in a format that makes it easy for
    fitting.
    
    Parameters
    ----------
    Blocks : list of DataBlock objects.
        Data to be used for fit.
    """

    def __init__(self, Blocks):
        self._data = self.preprocess_data(Blocks)
        # Figure out all dimensions.
        self.n_chan = self._data[0]['data'].shape[0]
        self.n_pol = self._data[0]['data'].shape[1]
        n_time = 0
        for data in self._data:
            s = data['data'].shape
            n_time += s[2]
            if s[0] != self.n_chan or s[1] != self.n_pol:
                msg = "Data blocks have conflicting dimensions."
                raise ValueError(msg)
        self.n_time = n_time
        self.n_scan = len(self._data)

    def get_data_weight_chan(self, chan_ind):
        """Get all data and weights for a single channel."""
        
        data_arr = np.empty((self.n_pol, self.n_time), dtype=np.float64)
        weights_arr = np.empty((self.n_pol, self.n_time), dtype=np.float64)
        for data, start, end in self.iterate_scans():
            this_chan_data = data['data'][chan_ind,:,:]
            this_chan_weights = data['weight'][chan_ind,:,:]
            data_arr[:,start:end] = this_chan_data
            weights_arr[:,start:end] = this_chan_weights
        return data_arr, weights_arr

    def get_all_data_weight(self):
        data_arr = np.empty((self.n_chan, self.n_pol, self.n_time),
                            dtype=np.float64)
        weights_arr = np.empty((self.n_chan, self.n_pol, self.n_time),
                               dtype=np.float64)
        for data, start, end in self.iterate_scans():
            this_chan_data = data['data']
            this_chan_weights = data['weight']
            data_arr[:,:,start:end] = this_chan_data
            weights_arr[:,:,start:end] = this_chan_weights
        return data_arr, weights_arr

    def get_all_field(self, field):
        """Get all of a specific data field."""

        format = self._data[0][field].dtype
        out = np.empty(self.n_time, dtype=format)
        for data, start, end in self.iterate_scans():
            this_field = data[field]
            out[start:end] = this_field
        return out

    def generate_scan_basis_polys(self, order):
        """Generate basis time polynomials for each scan up to order.
        
        Generate basis polynomials as a function of time for each scan 
        up to the passed maximum order.

        Parameters
        ----------
        order : integer
            Maximum order of the polynomials.

        Returns
        -------
        polys : ndarray of shape (n_scan, `order`, n_time)
        """

        polys = np.zeros((self.n_scan, order, self.n_time),
                         dtype=np.float64)
        scan_ind = 0
        for data, start, end in self.iterate_scans():
            time = data['time']
            this_poly = misc.ortho_poly(time, order)
            polys[scan_ind,:,start:end] = this_poly
            scan_ind +=1
        return polys

    def iterate_scans(self):
        """Returns an iterator over the data scans, with start and end time
        indeces."""

        class iterator(object):
            
            def __init__(self_itr):
                self_itr.scan = 0
                self_itr.time_ind = 0
                self_itr._data = self._data

            def __iter__(self_itr):
                return self_itr

            def next(self_itr):
                if len(self_itr._data) <= self_itr.scan:
                    raise StopIteration()
                start = self_itr.time_ind
                data = self._data[self_itr.scan]
                self_itr.scan += 1
                end = start + data['data'].shape[2]
                self_itr.time_ind = end
                return data, start, end
        
        return iterator()

    def preprocess_data(self, Blocks):
        """Puts data into convienient format and calculates weights."""
        
        # Where we will store formated data.
        data_list = []
        for Data in Blocks:
            # First do some basic checks on data format.
            if not tuple(Data.field['CRVAL4']) == (-5, -7, -8, -6):
                msg = "Expected data polarizations to be XX, XY, YX, YY."
                raise NotImplementedError(msg)
            if not Data.data.shape[2] == 1:
                msg = "Multiple cal states, ignoring all but first."
                raise Warning(msg)

            # Rearrange data.
            data = np.swapaxes(Data.data[:,:,0,:], 0, 2)
            mask = ma.getmaskarray(data).copy()
            data = data.filled(0)
            # XX and YY data should all be positive.
            mask[data[:,[0],:] < 0] = True
            mask[data[:,[3],:] < 0] = True
            # If any polarizations are masked, all polarizations should be.
            mask[np.any(mask, 1)[:,None,:]] = True
            data[mask] = 0
            
            # Make the noise weights.
            weight = data.copy()
            cross_norm = np.sqrt(weight[:,0,:] * weight[:,3,:])
            weight[:,[1,2],:] = cross_norm[:,None,:]
            weight[mask] = 1
            weight = 1 / weight
            weight[mask] = 0

            # Get the pointing.
            Data.calc_pointing()
            ra = Data.ra
            dec = Data.dec
            # Telescope frame.
            az = Data.field['CRVAL2']
            el = Data.field['CRVAL3']
            UT = Data.field['DATE-OBS']
            # Get the time array.
            Data.calc_time()
            time = Data.time
            
            # Store it all and add to list.
            this_data = {}
            this_data['data'] = data
            this_data['weight'] = weight
            this_data['ra'] = ra
            this_data['dec'] = dec
            this_data['az'] = az
            this_data['el'] = el
            this_data['time'] = time
            this_data['UT'] = UT
            data_list.append(this_data)
        return data_list

def fit_simple_gaussian(BeamData, Source):
    """Perform nonlinear fit of a gaussian to the Beam.
    
    Fits a simple, non-linear gaussian model to the data, using only the XX and
    YY polarizations.  Model includes Gaussian width, amplitude, centroid
    offset, and baseline system temperature.
    """
    
    # Get source relative coordinates for the data.
    UT = BeamData.get_all_field('UT')
    az_s, el_s = Source.azelGBT(UT)
    az_factor = np.cos(np.mean(el_s) * np.pi / 180)
    az = (BeamData.get_all_field('az') - az_s) * az_factor
    el = BeamData.get_all_field('el') - el_s
    
    # Construct the model function.
    def model(pars):
        az_off = pars[0]
        el_off = pars[1]
        width = pars[2]
        XX_amp = pars[3]
        YY_amp = pars[4]
        XX_Tsys = pars[5]
        YY_Tsys = pars[6]
        sigma = width / (2 * np.sqrt(2 * np.log(2)))
        gauss = np.exp(-((az - az_off)**2 + (el - el_off)**2) / (2 * sigma**2))
        model = np.array([XX_amp, YY_amp])[:,None] * gauss
        model += np.array([XX_Tsys, YY_Tsys])[:,None]
        return model
    
    # Allocate memory for outputs.
    all_pars = np.empty((BeamData.n_chan, 7), dtype=np.float64)
    # Initialization parameters.
    init_pars = [0, 0, 0.25, 10, 10, 10, 10]
    # Loop over channels and fit each one independantly.
    for ii in xrange(BeamData.n_chan):
        # Get the data for this channel.
        data, weights = BeamData.get_data_weight_chan(ii)
        # Use only XX and YY polarizations.
        data = data[[0,3],:]
        weights = weights[[0,3],:]
        # Construct the residual funciton.
        get_residuals = lambda p: ((data - model(p)) * weights).flat[:]
        fit_pars, ier = optimize.leastsq(get_residuals, init_pars, ftol=1e-5, 
                                         xtol=1e-5)
        all_pars[ii] = fit_pars
    
    # Unpack the parameters and return.
    source = all_pars[:,0:2].copy()
    width = all_pars[:,2].copy()
    amps = all_pars[:,3:5].copy()
    Tsys = all_pars[:,5:7].copy()
    return source, width, amps, Tsys


def linear_fit(BeamData, Basis, Source, beam_modes=2, time_modes=2):
    """Perform a linear fit to a set of basis function."""
    
    n_time = BeamData.n_time
    n_chan = BeamData.n_chan
    n_pol = BeamData.n_pol
    n_scan = BeamData.n_scan
    # Get source relative coordinates for the data.
    UT = BeamData.get_all_field('UT')
    az_s, el_s = Source.azelGBT(UT)
    az_factor = np.cos(np.mean(el_s) * np.pi / 180)
    az = (BeamData.get_all_field('az') - az_s) * az_factor
    el = BeamData.get_all_field('el') - el_s
    
    # Get time domain basis function for each scan. This is frequency
    # independant.
    time_basis = BeamData.generate_scan_basis_polys(time_modes)
    n_scan_parameters = n_scan * time_modes

    # Memory for holding basis functions.  The basis functions are frequency
    # channel dependant do to the dependance of the beam width and beam center.
    n_beam_basis = (beam_modes * (beam_modes + 1)) / 2
    beam_basis = np.empty((n_beam_basis, n_time), dtype=np.float64)
    # Number of parameters per polarization and per channel (i.e. per fit).
    n_parameters = n_beam_basis + n_scan_parameters
    # Memory for the coefficient matrix.
    coefficients = np.empty((n_time, n_parameters), dtype=np.float64)
    # Fill in the scan basis functions which are the same for each channel.
    scan_coeffs = np.reshape(time_basis, (n_scan_parameters, n_time)).T
    coefficients[:,n_beam_basis:] = scan_coeffs
    # Memory for output parameters.
    params = np.empty((n_chan, n_pol, n_parameters), dtype=np.float64)
    # Memory for storing the fit model data.
    model_data = np.empty((n_chan, n_pol, n_time), dtype=np.float64)

    # Loop over channels and fit each independantly.
    for ii in xrange(BeamData.n_chan):
        # Get the data for this channel.
        data, weights = BeamData.get_data_weight_chan(ii)
        # Evaluate the basis functions for this channel at the telescope
        # pointing.
        index = 0
        for jj in range(beam_modes):
            for kk in range(beam_modes - jj):
                ba_fun = Basis.eval_basis_chan((jj, kk), az, el, ii)
                beam_basis[index,:] = ba_fun
                index += 1
        # Construct the coefficient matrix, which is the same for each
        # polarization.
        coefficients[:,:n_beam_basis] = beam_basis.T
        # Perform fit for each polarization.
        for jj in range(n_pol):
            # Weight coefficients and data.
            this_data = data[jj,:] * weights[jj,:]
            this_coefficients = coefficients * weights[jj,:,None]
            # Perform the fit.
            this_params, residues, rank, s = linalg.lstsq(this_coefficients,
                                                          this_data)
            # TODO: can probably use the rank to figure out if there is any
            # data on that slice, do something extra with the gaps.
            params[ii,jj,:] = this_params
            model_data[ii,jj,:] = np.dot(coefficients, this_params)
    # Now unpack the parameters.
    scan_params = params[:,:,n_beam_basis:]
    scan_params.shape = (n_chan, n_pol, n_scan, time_modes)
    scan_params = np.rollaxis(scan_params, 2, 0)
    beam_params = np.zeros((n_chan, n_pol, beam_modes, beam_modes),
                           dtype=np.float64)
    index = 0
    for jj in range(beam_modes):
        for kk in range(beam_modes - jj):
            beam_params[:,:,jj,kk] = params[:,:,index]
            index += 1
    return beam_params, scan_params, model_data


