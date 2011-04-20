"""Module for calculating the map covarance matrix."""

import scipy as sp
from scipy import interpolate

from core import algebra
import beam
from kiyopy import parse_ini

# Put input parameters here along with thier default values.
params_init = {
    # Unit System ('transvers-radial') can be "deg-freq", "Mpc/h-Mpc/h", 
    # "deg-log_freq" or "Mpc/h-mu".
    "unit_system" : "deg-freq",
    # File to read the power spectrum from.
    "power_file_name" : "",
    # The format of the power spectrum to be read. Lost of these will be
    # incompatible with certain unit systems.
    "power_format" : "P(k)",
    # Map file from which to get the axis information.  
    # Perhaps should also have
    # the option to specify that covariance shape and such explicitly.
    "map_file" : "",
    # Output file name.
    "out_file_name" : ""
    }

prefix = 'cv_'

class Covariance(object) :
    """Driver program for covariance matrix calculation.
    
    Parameters
    ----------
    parameter_file_or_dict : str or dict
        The input parameters for the program which will be parsed by
        `kiyopy.parse_ini`.  See that module for a descrition of the input
        format.
    feedback : int
        Feedback level for the program.  Default is 2.
    """
    
    def __init__(self, parameter_file_or_dict=None, feedback=2) :
        # Read the parameter file, store in dictionary named parameters.
        self.params = parse_ini.parse(parameter_file_or_dict, params_init, 
                                      prefix=prefix, feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1) :
        """Function that does all the work.
        
        Parameters
        ----------
        nprocesses : int
            Number of threads to use.  So far this program is not threaded.
            this argument is only included for compatibility with other
            modules.
        """
        params = self.params

        #### Front end.  Read in the power spectrum, convert it to 2d correlation
        # function, etc.
        # corr_function = f(params['unit_system'], params['power_file_name'],
        #                   params['power_format'])

        # Temporary face correlation to make the rest of the code run.
        corr_function = lambda rho, f : (1.0e6/(rho**2 + 0.001)
                                         /(abs(f) + 1.0e3))

        #### Units.  Figure out the axes in the proper units etc.
        # Start by getting the map from which we will be getting our axes.
        map = algebra.open_memmap(params["map_file"], mode='r')
        map = algebra.make_vect(map)
        if map.axes != ('freq', 'ra', 'dec') :
            raise ce.DataError("Expeceted map to be in freq-ra-dec space.")
        # Next figure out all the axes.
        nfreq = map.shape[0]
        nra = map.shape[1]
        ndec = map.shape[2]
        freq = map.get_axis('freq')
        ra = map.get_axis('ra')*sp.cos(sp.pi*map.info['dec_centre']/180.0)
        dec = map.get_axis('dec')
        # Coordinate dependant conversion factors.
        system = params['unit_system']
        if system == "deg-freq" :
            z = freq
            z_width = sp.ones(nfreq)*map.info['freq_delta']
            Drho_Ddeg = sp.ones(nfreq)
        elif system == "Mpc/h-Mpc/h" :
            # Do a bunch of cosmology dependant stuff.
            # Not written yet, so bail here.
            return
        else :
            raise ValueError('Unit system must be one of "deg-freq", '
                             '"Mpc/h-Mpc/h", "deg-log_freq" or "Mpc/h-mu".')
        # Get the beam object.
        # The following is temporary.  Eventually need to read beam data from
        # file.
        beam_freq = sp.arange(600.0e6, 1000.0e6, 50.0e6)
        beam_width = 0.3*600.0e6/beam_freq
        Beam = beam.GaussianBeam(beam_width, beam_freq)
        # Figure out the range of angular lags we need to consider in our
        # coordinate system.
        # First the max lag.
        # This is inefficient if the range of Drho_Ddeg is very wide.
        max_lag = (max(ra) - min(ra))**2 + (max(dec) - min(dec))**2
        max_lag = sp.sqrt(max_lag)*max(Drho_Ddeg)
        # Figure out the minimum lag step.
        lag_step = min([abs(map.info['ra_delta']) *
                        sp.cos(sp.pi*map.info['dec_centre']/180.0), 
                        abs(map.info['dec_delta'])])
        if max_lag/lag_step > 10*(ndec + nra) :
            raise RunTimeError("Dynamic range is very large.  There will be "
                               "too many integrals to do.")
        # There is probably a more efficient lag binning than this... peicewise
        # linear then log?
        divisions = 10.0
        lag = sp.arange(0.0, max_lag + lag_step/divisions, lag_step/divisions)
        sq_lag = lag**2
        nlag = len(lag)

        #### Integral.  Loop over all possible lags and do the integral.
        # Allowcate memory for all f,f',lag combinations.
        integrals = sp.empty((nfreq, nfreq, nlag), dtype=float)
        if self.feedback >= 2 :
            print "Starting integrals."
        for find1 in xrange(nfreq) :
            for find2 in xrange(nfreq) :
                for lind in xrange(nlag) :
                    # Get separation in radial direction.
                    delta_z = abs(z[find1] - z[find2])
                    # Get the window functions.
                    W, rho_limits = Beam.angular_real_space_window(freq[find1],
                            freq[find2], return_limits=True)
                    Q, z_limits = Beam.radial_real_space_window(z_width[find1],
                            z_width[find2], return_limits=True)
                    # Construct the integrand.
                    int = integrand(corr_function, W, Q, lag[lind], delta_z)
                    # Integrate.
                    # Reiman sum is the most efficient integration algorithm
                    # known to mankind.
                    z_vals = sp.arange(z_limits[0], z_limits[1], 
                                       (z_limits[1] - z_limits[0])/20)
                    z_vals = z_vals[None, :]
                    rho_vals = sp.arange(rho_limits[0], rho_limits[1], 
                                       (rho_limits[1] - rho_limits[0])/20)
                    rho_vals = rho_vals[:, None]
                    result = (sp.sum(int(rho_vals, z_vals)) 
                              * (z_limits[1] - z_limits[0])/20
                              * (rho_limits[1] - rho_limits[0])/20)
                    # Store the result.
                    integrals[find1, find2, lind] = result
        if self.feedback >= 2 :
            print "Integrals done."

        #### Assignment.  Allocate memory, loop over elements and assign.
        # Allowcate memory and build final coraviance matrix object.
        covariance = algebra.open_memmap(params["out_file_name"], 
                        mode='w+', dtype=float, shape=(nfreq, nra, ndec, nfreq,
                        nra, ndec))
        covariance = algebra.make_mat(covariance, axis_names=("freq", "ra", 
                        "dec","freq", "ra", "dec"), 
                        row_axes=(0, 1, 2), col_axes=(3, 4, 5))
        covariance.copy_axis_info(map)
        # Make a matrix of angular pairwise lags.
        sq_angles = (dec[None, :, None, None] - dec[None, None, None, :])**2
        sq_angles = sq_angles + \
                (ra[:, None, None, None] - ra[None, None, :, None])**2
        if self.feedback >= 2 :
            print "Filling covariance matrix by interpolation."
        # Now fill in the elements by interpolating integrals.
        for find1 in xrange(nfreq) :
            for find2 in xrange(nfreq) :
                # The pairwise angular lags in the unit system.
                this_sq_lag = (sq_angles*(Drho_Ddeg[find1] + 
                                                Drho_Ddeg[find2])**2/4)
                # The interpolation function.  Perhaps there is a better
                # algorithm than cubic?
                interpolator = interpolate.interp1d(sq_lag,
                                integrals[find1,find2,:], bounds_error=True,
                                kind='cubic')
                covariance[find1,:,:,find2,:,:] = interpolator(this_sq_lag)
        del covariance, map
        if self.feedback >= 2 :
            print "Done"


def integrand(correlation_function, W, Q, delta_rho, delta_f) :
    """Returns integrand function for the covariance element.
    
    Constructs the integrand to be integrated to calculate an
    element of a pixel-pixel covariance matrix.  The result is a function of
    the transverse and radial wave numbers.  Many of the parameters described
    below are functions that factors in the returned integrand.  These
    are functions of either the transverse or radial lags.  These functions
    must be vectorized.

    All inputs must be brought to the same units before being passed.  The
    covaraince calculation can be done in either Mpc units or degree-hertz
    units, but all inputs must be consistant.

    Parameters
    ----------
    correlation_function : vectorized function (r, f)
        Correlation function to be converted to a covariance.  This function 
        will normally be constructed by interpolating a data set.
    W : vectorized function (r)
        Beam angular window function.
    Q: vectorized function (f)
        Beam radial window function.
    delta_rho : float
        Angular separation of the two pixels.  Units of Mpc/h or degrees.
    delta_f : float
        Radial separation of the two pixels.  Units of Mpc/h or Hertz.
        
    Returns
    -------
    f : vectorized function (r, f)
        Integrand for the covariance matrix.  The integral of this function
        over all space is a covariance matrix element.  Note that the aximuthal
        integral has already been performed and only the radial and transverse
        integrals remain.
    """

    return lambda r, f : ((2*sp.pi)*correlation_function(r+delta_rho,f+delta_f)
                          * W(r) * Q(f) * r)

# If this file is run from the command line, execute the main function.
if __name__ == "__main__":
    import sys
    Covariance(str(sys.argv[1])).execute()

