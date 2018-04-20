import numpy as np
import scipy as sp

#Takes Parkes beam fit data, which has some odd sharp features, and replaces them with interpolations on the remaining data.

#Currently replaces data points more than 30% away from 0.26 degree FWHM.

def interpolate_on_outliers(FWHMs, tol = 0.3, FWHM = 0.26):
    orig = np.array(FWHMs)
    freqs = np.arange(len(orig))
    orig[np.abs((orig/FWHM)-1)>tol] = 0
    good_points = np.nonzero(orig)
    interp = sp.interpolate.interp1d(freqs[good_points],orig[good_points], fill_value = 'extrapolate')
    #interp = sp.interpolate.UnivariateSpline(freqs[good_points],orig[good_points])
    new_FWHMs = interp(freqs)
    return new_FWHMs

def interpolate_on_steepness(FWHMs, slope_max = 0.002):
    orig = np.array(FWHMs)
    freqs = np.arange(len(orig))
    slopes = np.diff(orig)/np.diff(freqs)
    slope_keep = np.abs(slopes)<slope_max
    bool = np.append(slope_keep[0],slope_keep)
    interp = sp.interpolate.interp1d(freqs[bool],orig[bool], fill_value = 'extrapolate')
    new_FWHMs = interp(freqs)
    return new_FWHMs

def butter_lowpass_filter(FWHMs, butter_order=5, butter_cutoff=0.1):
    b, a = sp.signal.butter(butter_order, butter_cutoff, 'low')
    new_FWHMs = sp.signal.filtfilt(b, a, FWHMs)
    return new_FWHMs
