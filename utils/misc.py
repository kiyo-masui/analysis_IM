"""Some general utilities that are useful."""

import datetime
import time
import sys

import copy
import scipy as sp
import numpy as np
from scipy.interpolate import interp1d
from scipy import linalg, special
import ephem
from scipy.stats import chisqprob

def elaz2radec_lst(el, az, lst, lat = 38.43312) :
    """DO NOT USE THIS ROUTINE FOR ANTHING THAT NEEDS TO BE RIGHT.  IT DOES NOT
    CORRECT FOR PRECESSION.

    Calculates the Ra and Dec from elavation, aximuth, LST and Latitude.

    This function is vectorized with numpy so should be fast.  Standart numpy
    broadcasting should also work.

    All angles in degrees, lst in seconds. Latitude defaults to GBT.
    """

    # Convert everything to radians.
    el = sp.radians(el)
    az = sp.radians(az)
    lst = sp.array(lst, dtype = float)*2*sp.pi/86400
    lat = sp.radians(lat)
    # Calculate dec.
    dec = sp.arcsin(sp.sin(el)*sp.sin(lat) +
                    sp.cos(el)*sp.cos(lat)*sp.cos(az))
    # Calculate the hour angle
    ha = sp.arccos((sp.sin(el) - sp.sin(lat)*sp.sin(dec)) /
                   (sp.cos(lat)*sp.cos(dec)))
    ra = sp.degrees(lst - ha) % 360

    return ra, sp.degrees(dec)

def get_ephem_GBT(UT=None):
    GBT = ephem.Observer()
    GBT.long = '-79:50:23.4'
    GBT.lat = '38:25:59.23'
    # According to the document below, corrections for refraction should
    # already be present.
    # http://www.gb.nrao.edu/GBT/MC/doc/dataproc/gbtAntFits/gbtAntFits.pdf
    GBT.pressure = 0 # Turn off refraction.
    GBT.epoch = ephem.J2000
    if UT:
        GBT.date = UT2ephem_date(UT)
    return GBT

def UT2ephem_date(UT):
    UT_wholesec, partial_sec = UT.split('.', 1)
    time_obj = time.strptime(UT_wholesec, "%Y-%m-%dT%H:%M:%S")
    UT_reformated = time.strftime("%Y/%m/%d %H:%M:%S", time_obj)
    return UT_reformated + "." + partial_sec

def azel2radecGBT(az, el, UT):
    """Calculates the Ra and Dec from the azimuth, elavation and UT for an
    observer at GBT.

    All input should be formated to correspond to the data in a GBT fits file.
    El, and Az in degrees and UT a string like in the GBT DATE-OBS field.

    Largely copied from Kevin's code.
    """

    GBT = get_ephem_GBT(UT)

    el_r = el*sp.pi/180.0
    az_r = az*sp.pi/180.0
    ra, dec = GBT.radec_of(az_r,el_r)

    return ra*180.0/sp.pi, dec*180.0/sp.pi

def elaz2radecGBT(el, az, UT):
    return azel2radecGBT(az, el, UT)

def radec2azelGBT(ra, dec, UT):
    """Calculates the Ra and Dec from the elevation, azimuth and UT for an
    observer at GBT.

    All input should be formated to correspond to the data in a GBT fits file.
    El, and Az in degrees and UT a string like in the GBT DATE-OBS field.

    Largely copied from Kevin's code.
    """

    GBT = get_ephem_GBT(UT)

    source = ephem.FixedBody()
    source._ra = ra * sp.pi / 180
    source._dec = dec * sp.pi / 180
    source._epoch = ephem.J2000
    
    source.compute(GBT)
    return source.az * 180 / sp.pi, source.alt * 180 / sp.pi

def LSTatGBT(UT) :
    """Calculates the LST from the UT of an observer at GBT.

    All input should be formated to correspond to the data in a GBT fits file.
    UT a string like in the GBT DATE-OBS field.

    Largely copied from Kevin's code.
    """

    GBT = get_ephem_GBT(UT)
    LST = GBT.sidereal_time() #IN format xx:xx:xx.xx ?
    return LST*180.0/sp.pi

def azel2pGBT(az, el, UT):
    """Converts the azimuth and elevation to a parallactic angle for GBT.
    
    This function is includes precession, but probably isn't
    accurate to the nutation level.
    """
    
    # Step size for finite difference, degrees.
    delta = 0.1
    # Get ra and dec at the centre, right and left of the point.
    ra_c, dec_c = azel2radecGBT(az, el, UT)
    ra_r, dec_r = azel2radecGBT(az - delta, el, UT)
    ra_l, dec_l = azel2radecGBT(az + delta, el, UT)
    # Get the deltas for the 2 coordinates, in real degrees on sky.
    delta_ra = (ra_l - ra_r) * np.cos(dec_c * np.pi / 180)
    delta_dec = dec_l - dec_r
    # Calculate the Paralactic angle.
    p_rad = np.arctan2(delta_dec, -delta_ra)
    return p_rad * 180 / np.pi

def radec_to_sexagesimal(ra, dec):
    """
    Accepts float values for ra and dec
    returns ra, dec in the form (hh,mm,ss), abs(deg,mm,ss), sign
    where the absolute value of the declination is returned and the
    sign is the sign of the declination ('+' or '-')
    """
    if ra < 0:
        ra += 360.

    ra_hr = int(ra / 15.)
    minsec = np.mod(ra, 15.) / 15. * 60
    ra_min = int(minsec)
    ra_sec = np.mod(minsec, 1.) * 60

    if dec >= 0.:
        sign = '+'
    else:
        sign = '-'

    dec_deg = int(abs(dec))
    minsec = np.mod(abs(dec), 1.) * 60
    dec_min = int(minsec)
    dec_sec = np.mod(minsec, 1.) * 60
    return (ra_hr, ra_min, ra_sec), (dec_deg, dec_min, dec_sec), sign

def sexagesimal_to_radec(ra, dec):
    """
    Assumes ra, dec arguments are of the form "hh:mm:ss", "deg:mm:ss"
    Returns float values for ra and dec
    """
    ra_hr, ra_min, ra_sec= ra.split(':')

    sign = ra_hr[0]
    if sign == '-':
        sign_ra = -1.
        ra_hr = ra_hr[1:]
    else:
        sign_ra = 1.

    ra = sign_ra * 15. * ( float(ra_hr) + float(ra_min)/60. + \
                           float(ra_sec)/3600. )

    dec_deg, dec_min, dec_sec= dec.split(':')

    sign = dec_deg[0]
    if sign == '-':
        sign_dec = -1.
        dec_deg  = dec_deg[1:]
    else:
        sign_dec = 1.

    dec =  sign_dec * (float(dec_deg) + float(dec_min)/60. + \
                       float(dec_sec)/3600.)

    return ra, dec

def time2float(UT) :
    """Calculates float seconds from a time string.

    Convert a time string in format %Y-%m-%dT%H:%M:%S.partial to a float number
    of seconds ignoring many corrections (such as leap seconds).  However it
    should be internally consistant, and it should be safe to use this method
    to convert to a float, do some operations, then use `float2time` to convert
    back.  Works for array inputs."""

    if not isinstance(UT, sp.ndarray) :
        UT = sp.array([UT])
        single_time = True
    else:
        single_time = False

    time_array = sp.empty(UT.shape, dtype=float)
    
    for ii in xrange(UT.size):
        split_UT = UT.flat[ii].split('.', 1)
        if len(split_UT) == 2:
            UT_wholesec, partial_sec = split_UT
        else:
            UT_wholesec = split_UT[0]
            partial_sec = '0'
        to = datetime.datetime(*time.strptime(UT_wholesec,
                                              "%Y-%m-%dT%H:%M:%S")[:6])
        epoch_start = datetime.datetime(2000, 1, 1)
        td = to - epoch_start
        td_seconds = ((td.microseconds + (td.seconds + td.days * 24 * 3600) *
                       10**6) / 10**6)
        time_array.flat[ii] = td_seconds + float("0." + partial_sec)
    if single_time :
        return time_array[0]
    else :
        return time_array


def float2time(t) :
    """Does the reverse operation of `time2float`."""

    if not isinstance(t, sp.ndarray) :
        t = sp.array([t])
        single_time = True
    else:
        single_time = False
    
    time_str_array = sp.empty(t.shape, dtype='S22')

    for ii in range(t.size) :
        td = datetime.timedelta(seconds=t.flat[ii])
        epoch_start = datetime.datetime(2000, 1, 1)
        time_obj = epoch_start + td
        
        time_str = time_obj.strftime("%Y-%m-%dT%H:%M:%S.%f")
        # Truncate the fraction seconds to hundredths.
        time_str = time_str[:22]
        time_str_array.flat[ii] = time_str
    if single_time :
        return time_str_array[0]
    else :
        return time_str_array


def mk_map_grid(centre, shape, spacing) :
    """Make a grid of coordinates in Ra and Dec.

    This function accepts a field centre (tuple 2 floats), map-shape (tuple 2
    ints) and pixel spacing (float) for a map.  It returns two arrays with
    the specified map shape, and meshed like sp.meshgrid.  However, the spacing
    in Ra is converted to real degrees (devided by cos(dec) at field centre).

    All units in degrees, values are pixel centres.
    """

    dec = centre[1] + spacing*sp.arange(-(shape[1]-1.)/2., shape[1]/2.)
    ra = centre[0] + (spacing/sp.cos(centre[1]*sp.pi/180.) *
                      sp.arange(-(shape[0]-1.)/2., shape[0]/2.))

    grid_ra, grid_dec = sp.meshgrid(ra, dec)

    return grid_ra, grid_dec


def get_beam(freq) :
    """Get the GBT beam width at a frequency (or an array of frequencies).

    This is currently pretty rough and only uses on scans worth of data.
    """

    # This data is pretty rough.  Just cut and pasted from one scan, not
    # averaged.
    beam_data = [0.316148488246, 0.306805630985, 0.293729620792,
                 0.281176247549, 0.270856788455, 0.26745856078,
                 0.258910010848, 0.249188429031]
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905], dtype=float)
    freq_data *= 1.0e6
    f = interp1d(freq_data, beam_data, bounds_error=False, fill_value = -1)
    b = f(freq)
    b[b<0] = 0.316148488246
    return b


def polint2str(pol_int) :
    """Convert an interger representing a polarization to a representing the
    polarization.

    This is based on the SDfits convention that I pulled from: 
        https://safe.nrao.edu/wiki/bin/view/Main/SdfitsDetails

    Here are the return values based on the passed integer.

    RR  -1
    LL 	-2
    RL 	-3
    LR 	-4
    XX 	-5
    YY 	-6
    XY 	-7
    YX 	-8
    I 	1
    Q 	2
    U 	3
    V 	4
    Otherwise raises a ValueError.
    """

    if pol_int == -1 :
        return 'RR'
    elif pol_int == -2 :
        return 'LL'
    elif pol_int == -3 :
        return 'RL'
    elif pol_int == -4 :
        return 'LR'
    elif pol_int == -5 :
        return 'XX'
    elif pol_int == -6 :
        return 'YY'
    elif pol_int == -7 :
        return 'XY'
    elif pol_int == -8 :
        return 'YX'
    elif pol_int == 1 :
        return 'I'
    elif pol_int == 2 :
        return 'Q'
    elif pol_int == 3 :
        return 'U'
    elif pol_int == 4 :
        return 'V'
    else :
        raise ValueError("Polarization integer must be in range(-8, 5) and "
                         "nonzero")

#<<<<<<< HEAD
#=======

#>>>>>>> master
def ampfit(data, covariance, theory, rank_thresh=1e-12, diag_only=False):
    """Fits the amplitude of the theory curve to the data.

    Finds `amp` such that `amp`*`theory` is the best fit to `data`.

    Returns
    -------
    amp : float
        Fitted amplitude.
    error : float
        Error on fitted amplitude.
    """

    data = sp.asarray(data)
    covariance = sp.asarray(copy.deepcopy(covariance))
    theory = sp.asarray(theory)

    if len(data.shape) != 1:
        raise ValueError("`data` must be a 1D vector.")

    n = len(data)

    if data.shape != theory.shape:
        raise ValueError("`theory` must be the same shape as `data`.")

    if covariance.shape != (n,n):
        msg = "`covariance` must be a square matrix compatible with data."
        raise ValueError(msg)

    if diag_only:
        covariance = sp.diag(sp.diag(covariance))
        print data
        print sp.diag(sp.sqrt(covariance))

    covariance_inverse = linalg.inv(covariance)
    weighted_data = sp.dot(covariance_inverse, data)
    amp = sp.dot(theory, weighted_data)
    normalization = sp.dot(covariance_inverse, theory)
    normalization = sp.dot(theory, normalization)
    amp /= normalization
    error = sp.sqrt(1/normalization)

    u, s, v = np.linalg.svd(covariance)
    dof = np.sum(s > rank_thresh)
    resid = data - amp * theory
    chi2 = sp.dot(covariance_inverse, resid)
    chi2 = sp.dot(resid, chi2)
    pte = sp.stats.chisqprob(chi2, dof - 1)

    return {"amp":amp, "error":error, \
            "chi2":chi2, "dof":dof - 1, \
            "pte":pte}

def rebin_1D(array, reduce=4, axis=-1):
    shape = array.shape
    if axis < 0:
        axis = array.ndim + axis
    array.shape = shape[:axis] + (shape[axis]//reduce, reduce) + shape[axis+1:]
    out = sp.zeros(shape[:axis] + (shape[axis]//reduce,) + shape[axis+1:],
                   dtype=array.dtype)
    for ii in range(reduce):
        index = [slice(None),]*array.ndim
        index[axis+1] = ii
        index = tuple(index)
        out += array[index]
    out /= reduce
    array.shape = shape
    return out


def ortho_poly(x, n, window=1., axis=-1):
    """Generate orthonormal basis polynomials.

    Generate the first `n` orthonormal basis polynomials over the given domain
    and for the given window using the Gram-Schmidt process.
    
    Parameters
    ----------
    x : 1D array length m
        Functional domain.
    n : integer
        number of polynomials to generate. `n` - 1 is the maximum order of the
        polynomials.
    window : 1D array length m
        Window (weight) function for which the polynomials are orthogonal.

    Returns
    -------
    polys : n by m array
        The n polynomial basis functions. Normalization is such that
        np.sum(polys[i,:] * window * polys[j,:]) = delta_{ij}
    """
    
    if np.any(window < 0):
        raise ValueError("Window function must never be negative.")
    # Check scipy versions. If there is a stable polynomial package, use it.
    s_ver = sp.__version__.split('.')
    major = int(s_ver[0])
    minor = int(s_ver[1])
    if major <= 0 and minor < 8:
        new_sp = False
        if n > 20:
            raise NotImplementedError("High order polynomials unstable.")
    else:
        new_sp = True
    # Get the broadcasted shape of `x` and `window`.
    # The following is the only way I know how to get the broadcast shape of
    # x and window.
    # Turns out I could use np.broadcast here.  Fix this later.
    shape = np.broadcast(x, window).shape
    m = shape[axis]
    # Construct a slice tuple for up broadcasting arrays.
    upbroad = [slice(sys.maxsize)] * len(shape)
    upbroad[axis] = None
    upbroad = tuple(upbroad)
    # Allocate memory for output.
    polys = np.empty((n,) + shape, dtype=float)
    # For stability, rescale the domain.
    x_range = np.amax(x, axis) - np.amin(x, axis)
    x_mid = (np.amax(x, axis) + np.amin(x, axis)) / 2.
    x = (x - x_mid[upbroad]) / x_range[upbroad] * 2
    # Reshape x to be the final shape.
    x = np.zeros(shape, dtype=float) + x
    # Now loop through the polynomials and construct them.
    # This array will be the starting polynomial, before orthogonalization
    # (only used for earlier versions of scipy).
    if not new_sp:
        basic_poly = np.ones(shape, dtype=float) / np.sqrt(m)
    for ii in range(n):
        # Start with the basic polynomial.
        # If we have an up-to-date scipy, start with nearly orthogonal
        # functions.  Otherwise, just start with the next polynomial.
        if not new_sp:
            new_poly = basic_poly.copy()
        else:
            new_poly = special.eval_legendre(ii, x)
        # Orthogonalize against all lower order polynomials.
        for jj in range(ii):
            new_poly -= (np.sum(new_poly * window * polys[jj,:], axis)[upbroad]
                         * polys[jj,:])
        # Normalize, accounting for possibility that all data is masked. 
        norm = np.array(np.sqrt(np.sum(new_poly**2 * window, axis)))
        if norm.shape == ():
            if norm == 0:
                bad_inds = np.array(True)
                norm = np.array(1.)
            else:
                bad_inds = np.array(False)
        else:
            bad_inds = norm == 0
            norm[bad_inds] = 1.
        new_poly /= norm[upbroad]
        new_poly[bad_inds[upbroad]] = 0
        # Copy into output.
        polys[ii,:] = new_poly
        # Increment the base polynomial with another power of the domain for
        # the next iteration.
        if not new_sp:
            basic_poly *= x
    return polys

def ortho_poly_2D(x, y, n, window=1):
    """Generate a 2D orthonormal polynomial basis up to order `n - 1`.

    Generate the first `n(n + 1)/2` orthonormal basis polynomials over 
    the given domain and for the given window using the Gram-Schmidt 
    process.
    
    Parameters
    ----------
    x : 2D array
        Functional domain, x coordinate.
    y : 2D array
        Functional domain, y coordinate. Shape must be the same as `x` or
        broadcastable to the same shape.
    n : integer
        number of polynomials to generate. `n - 1` is the maximum order of the
        polynomials.
    window : 2D array
        Window (weight) function for which the polynomials are orthogonal.
        Shape must be the same as `x` or broadcastable to the same shape.

    Returns
    -------
    polys : array of shape (n, n) + x.shape
        The `n(n + 1)/2` polynomial basis functions. Normalization is
        such that np.sum(polys[i, j,:] * window * polys[k, l,:]) = delta_{ik}
        delta_{jl}. Note that half the array is not used and is left empty. 
    """
    
    # Check input shapes
    b = np.broadcast(x, y, window)
    if b.nd != 2:
        raise ValueError("Inputs not broadcastable to a 2D array.")

    # This doesn't actually work and it might not be possible to make it work.
    msg = "This function is not working and the theory behind it is dubious."
    raise NotImplementedError(msg)

    out = np.zeros((n, n) + b.shape, dtype=float)
    # Check scipy versions. If there is a stable polynomial package, use it.
    s_ver = sp.__version__.split('.')
    major = int(s_ver[0])
    minor = int(s_ver[1])
    if major <= 0 and minor < 8:
        new_sp = False
        if n > 20:
            raise NotImplementedError("High order polynomials unstable.")
    else:
        new_sp = True
    # For stability, rescale the domain.
    x_range = np.amax(x) - np.amin(x)
    x_mid = (np.amax(x) + np.amin(x)) / 2.
    x = (x - x_mid) / x_range * 2
    y_range = np.amax(y) - np.amin(y)
    y_mid = (np.amax(y) + np.amin(y)) / 2.
    y = (y - y_mid) / y_range * 2
    # Initialize basic polynomials for x and y domains if using old scipy.
    if not new_sp:
        basic_x = np.ones_like(x)
    # Loop over the polynomial indices and generate them.
    for ii in range(n):
        if not new_sp:
            basic_y = np.ones_like(y)
        # If using the new scipy, start with stably evaluated Legendre.
        # polynomial.
        if new_sp:
            basic_x = special.eval_legendre(ii, x)
        for jj in range(n - ii):
            if new_sp:
                basic_y = special.eval_legendre(jj, y)
            # This polynomial begins as the product of the basic polynomials.
            new_poly = basic_x * basic_y
            # Orthogonalize against all lower order polynomials.
            #for kk in range(ii):
            #    for mm in range(jj):
            #        new_poly -= out[kk,mm] * np.sum(out[kk,mm] * window *
            #                                        new_poly)
            # Normalize.
            new_poly /= np.sqrt(np.sum(new_poly**2))
            # Copy to output array.
            out[ii,jj,...] = new_poly
            # If using old scipy, update the basic polynomial to the next
            # order.
            if not new_sp:
                basic_y *= y
        if not new_sp:
            basic_x *= x
    return out



#<<<<<<< HEAD
#=======
if __name__ == "__main__":
    import doctest

    # run some tests
    OPTIONFLAGS = (doctest.ELLIPSIS |
                   doctest.NORMALIZE_WHITESPACE)
    doctest.testmod(optionflags=OPTIONFLAGS)

#>>>>>>> master
