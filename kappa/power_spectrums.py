import numpy as np
import core.algebra as al
from quadratic_products.pwrspec_estimator import cross_power_est


def power_specs(map_s, bins, m1=None, m2=None, verbose=True):
    r"""Bin different power spectrums for the result of the estimator,
    or the same power spectrums for maps m1, m2
    
    Input parameters :
        - map_s : KappaEstimator object, result of the estimator
        - bins : boundaries of bins in k space used for binning 1D power spectra
        - m1 : optional, replaces overdensity map delta
        - m2 : optional, replaces kappa
    Output : dictionary containing the following data :
        - p_delta : delta power spectrum
        - p_kappa : reconstructed kappa power spectrum
        - crossp : delta x kappa cross correlation spectrum
        - r : delta x kappa correlation function
    """

    if m1 is None: m1 = map_s.map
    m1 = al.make_vect(m1)
    m1.info = map_s.info
    if m2 is None: m2 = map_s.clean
    m2 = al.make_vect(m2)
    m2.info = map_s.info
    weight = np.ones_like(map_s.map)

    if verbose: print "Computing power spectrums..."
    pow_delta = cross_power_est(m1, m1, weight, weight, window=None, nonorm=True)
    if verbose: print "pow_delta : done"
    pow_kappa = cross_power_est(m2, m2, weight, weight, window=None, nonorm=True)
    if verbose: print "pow_kappa : done"
    crossp = cross_power_est(m1, m2, weight, weight, window=None, nonorm=True)
    if verbose: print "crossp : done"

    kz = map_s.get_k_axis('z')
    kx = map_s.get_k_axis('x')
    ky = map_s.get_k_axis('y')
    l, m, n = map_s.map.shape
    rad_3d = (kz*kz)[:, None, None] * np.ones((1, m, n))
    rad_3d += (kx*kx)[None, :, None] * np.ones((l, 1, n))
    rad_3d += (ky*ky)[None, None, :] * np.ones((l, m, 1))
    rad = np.sqrt(rad_3d)
    count = np.histogram(rad, bins=bins)[0]

    if verbose: print "Binning..."
    p_delta = np.histogram(rad, bins=bins, weights=pow_delta)[0]
    if verbose: print "pow_delta : done"
    p_kappa = np.histogram(rad, bins=bins, weights=pow_kappa)[0]
    if verbose: print "pow_kappa : done"
    c_p = np.histogram(rad, bins=bins, weights=crossp)[0]
    if verbose: print "crossp : done"

    res_dict = {'p_delta': p_delta/count,
                'p_kappa': p_kappa/count,
                'crossp': c_p/count}
    res_dict['r'] = res_dict['crossp']/np.sqrt(res_dict['p_delta']*res_dict['p_kappa'])
    return res_dict
