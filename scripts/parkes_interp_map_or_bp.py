import numpy as np
import scipy
import scipy.interpolate
import core.algebra as al
import copy

#Interpolates parkes maps or bandpasses that have been binned to lower freq. res.

def interpolate_1d(data, factor = 16.):
    #Interpolate data, presumably a function of frequency.
    #Increase frequency sampling by factor=factor.
    #For now using linear interp., since extrapolation doesn't work with spline.
    f_dim = len(data)
    freq_orig = np.arange(f_dim)
    freq_new = np.arange(0,f_dim,1./factor)
    interp = scipy.interpolate.interp1d(freq_orig, data, fill_value = 'extrapolate')
    #print freq_orig
    #print freq_new
    #print data
    data_new = interp(freq_new)
    return data_new

def interpolate_bp(bp_path, factor = 16.):
    bp_orig = np.load(bp_path)
    bp_new = np.zeros((bp_orig.shape[0],bp_orig.shape[1]*factor))
    for ind in range(bp_orig.shape[0]):
        bp_new[ind,:] = interpolate_1d(bp_orig[ind,:],factor)
        #Renormalize (keep bp mean the same)
        #bp_new[ind,:] *= np.mean(bp_orig[ind,:])/np.mean(bp_new[ind,:])
    return bp_new

def interpolate_map(map_path, factor = 16.):
    #Loop for now.  There may be a faster way.
    map = al.make_vect(al.load(map_path))
    new_map = np.zeros((map.shape[0]*factor, map.shape[1], map.shape[2]))
    new_map = al.info_array(new_map)
    new_map.info = copy.deepcopy(map.info)
    orig_freq_centre = map.info['freq_centre']
    new_map.info['freq_delta'] /= factor
    print map.info['freq_delta']
    print new_map.info['freq_delta']
    #map_centre is the frequency centre of index n//2, so it must shift.
    new_map.info['freq_centre'] = orig_freq_centre - (map.info['freq_delta'] - new_map.info['freq_delta'])/2.
    for pair in ([x,y] for x in range(map.shape[1]) for y in range(map.shape[2])):
         new_map[:, pair[0], pair[1]] = interpolate_1d(map[:,pair[0],pair[1]],factor)
    return new_map
