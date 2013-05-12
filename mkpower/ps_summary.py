#! /usr/bin/env python 

import numpy as np
from core import algebra
from simulations import corr21cm

def load_theory_ps(k_bins_centre):
    c21 = corr21cm.Corr21cm()
    #ps_theory  = c21.T_b(0.8)**2*c21.get_pwrspec(k_bins_centre)
    #ps_theory *= k_bins_centre**3/2./np.pi**2

    power_th = np.loadtxt('/Users/ycli/Code/analysis_IM/simulations/data/wigglez_halofit_z0.8.dat')
    power_th[:,0] *= 0.72
    ps_theory  = (c21.T_b(0.8)*1.e-3)**2*power_th[:,1]/0.72**3
    ps_theory *= power_th[:,0]**3/2./np.pi**2
    return ps_theory, power_th[:,0]

def get_2d_k_bin_centre(ps_2d):

    k_p = np.arange(ps_2d.shape[0]) - ps_2d.shape[0]//2
    k_p = ps_2d.info['k_p_centre'] * ps_2d.info['k_p_delta']**k_p

    k_v = np.arange(ps_2d.shape[1]) - ps_2d.shape[1]//2
    k_v = ps_2d.info['k_v_centre'] * ps_2d.info['k_v_delta']**k_v

    return k_p, k_v

def get_2d_k_bin_edges(ps_2d):

    k_p, k_v = get_2d_k_bin_centre(ps_2d)

    k_p_edges  = np.append(k_p, k_p[-1]*ps_2d.info['k_p_delta'])
    k_p_edges /= np.sqrt(ps_2d.info['k_p_delta'])

    k_v_edges  = np.append(k_v, k_v[-1]*ps_2d.info['k_v_delta'])
    k_v_edges /= np.sqrt(ps_2d.info['k_v_delta'])

    return k_p_edges, k_v_edges

def load_transfer_function(rf_root, tr_root, auto=True):
    
    ps_rf = algebra.make_vect(algebra.load(rf_root))
    ps_tr = algebra.make_vect(algebra.load(tr_root))

    ps_tr[ps_tr==0] = np.inf

    transfer_function = ps_rf/ps_tr

    if auto:
        transfer_function **=2

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_rf)

    return transfer_function, k_p_edges, k_v_edges

def load_weight(ns_root, transfer_function=None):
    
    ps_ns = algebra.make_vect(algebra.load(ns_root))
    kn_ns = algebra.make_vect(algebra.load(ns_root.replace('pow', 'kmn')))

    kn_ns[kn_ns==0] = np.inf

    gauss_noise = 2. * ps_ns / np.sqrt(kn_ns)
    if transfer_function != None:
        gauss_noise *= transfer_function

    gauss_noise[gauss_noise==0] = np.inf
    weight = (1./gauss_noise)**2

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_ns)

    return weight, k_p_edges, k_v_edges

def load_power_spectrum_err(ps_root):
    
    ps_2derr = algebra.make_vect(algebra.load(ps_root.replace('pow', 'err')))

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2derr)

    return ps_2derr, k_p_edges, k_v_edges

def load_power_spectrum(ps_root):
    
    ps_2d = algebra.make_vect(algebra.load(ps_root))
    kn_2d = algebra.make_vect(algebra.load(ps_root.replace('pow', 'kmn')))

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2d)

    return ps_2d, kn_2d, k_p_edges, k_v_edges

def convert_2dps_to_1dps_sim(rf_root, ns_root=None, truncate_range=None ):

    power_spectrum = algebra.make_vect(algebra.load(rf_root))
    power_spectrum_err = algebra.make_vect(algebra.load(rf_root.replace('pow', 'err')))
    k_mode_number = algebra.make_vect(algebra.load(rf_root.replace('pow', 'kmn')))
    if ns_root != None:
        weight, k_p_edges, k_v_edges = load_weight(ns_root )
        #k_mode_number *= weight
        k_mode_number = weight

    k_p_edges, k_v_edges = get_2d_k_bin_edges(power_spectrum)
    k_p_centre, k_v_centre = get_2d_k_bin_centre(power_spectrum)
    k_centre = np.sqrt(k_p_centre[:,None]**2 + k_v_centre[None, :]**2)
    k_edges_1d = k_p_edges

    power_spectrum *= k_mode_number
    power_spectrum_err *= k_mode_number

    power_spectrum_err **= 2

    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        power_spectrum_err, kpe, kve = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, power_spectrum_err)
        power_spectrum, kpe, kve     = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, power_spectrum)
        k_centre, kpe, kve           = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, k_centre)
        k_mode_number, kpe, kve      = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, k_mode_number)

    k_centre = k_centre.flatten()

    #k_mode,k_e = np.histogram(k_centre, k_edges_1d, weights=np.ones_like(k_centre))
    k_mode,k_e = np.histogram(k_centre, k_edges_1d, weights=k_mode_number.flatten())
    k_mode_2,k_e = np.histogram(k_centre, k_edges_1d, weights=(k_mode_number**2).flatten())
    power, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum.flatten())
    err,   k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum_err.flatten())
    k_mode[k_mode==0] = np.inf
    k_mode_2[k_mode_2==0] = np.inf
    power_1d = power/k_mode
    power_1d_err = np.sqrt(err/k_mode_2)
    return power_1d, power_1d_err, k_p_centre

def truncate_2dps(k_p_range, k_v_range, k_p_edges, k_v_edges, power_spectrum):
    
    k_p_range_idx = np.digitize(k_p_range, k_p_edges) - 1
    k_v_range_idx = np.digitize(k_v_range, k_v_edges) - 1

    k_p_range_idx[k_p_range_idx<0] = 0
    k_v_range_idx[k_v_range_idx<0] = 0
    #print k_p_range_idx
    #print k_v_range_idx

    return power_spectrum[k_p_range_idx[0]:k_p_range_idx[1],
                          k_v_range_idx[0]:k_v_range_idx[1]],\
           k_p_edges[k_p_range_idx[0]:k_p_range_idx[1]+1],\
           k_v_edges[k_v_range_idx[0]:k_v_range_idx[1]+1]

def convert_2dps_to_1dps(ps_root, ns_root, tr_root, rf_root, truncate_range=None):

    transfer_function, k_p_edges, k_v_edges = load_transfer_function(rf_root, tr_root)

    weight, k_p_edges, k_v_edges = load_weight(ns_root, transfer_function)

    power_spectrum, k_mode_number, k_p_edges, k_v_edges = load_power_spectrum(ps_root)

    power_spectrum_err, k_p_edges, k_v_edges = load_power_spectrum_err(ps_root)

    #k_mode_number *= weight

    power_spectrum *= transfer_function
    power_spectrum *= weight
    #power_spectrum *= k_mode_number

    power_spectrum_err *= transfer_function
    power_spectrum_err *= weight
    #power_spectrum_err *= k_mode_number
    power_spectrum_err **= 2

    k_p_centre, k_v_centre = get_2d_k_bin_centre(power_spectrum)
    k_centre = np.sqrt(k_p_centre[:,None]**2 + k_v_centre[None, :]**2)

    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        power_spectrum_err, kpe, kve = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, power_spectrum_err)
        power_spectrum, kpe, kve     = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, power_spectrum)
        k_centre, kpe, kve           = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, k_centre)
        weight, kpe, kve             = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, weight)
        k_mode_number, kpe, kve      = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                                      k_v_edges, k_mode_number)

    #print power_spectrum.shape
    power_spectrum_err = power_spectrum_err.flatten()
    power_spectrum     = power_spectrum.flatten()
    k_centre           = k_centre.flatten()
    weight             = weight.flatten()
    k_mode_number      = k_mode_number.flatten()

    k_edges_1d = k_p_edges

    #k_mode,k_e = np.histogram(k_centre, k_edges_1d, weights=np.ones_like(k_centre))
    #k_mode, k_e = np.histogram(k_centre, k_edges_1d, weights=k_mode_number)
    #k_mode_2, k_e = np.histogram(k_centre, k_edges_1d, weights=(k_mode_number**2))
    normal, k_e = np.histogram(k_centre, k_edges_1d, weights=weight)
    normal_2, k_e = np.histogram(k_centre, k_edges_1d, weights=(weight**2))
    power, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum)
    err, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum_err)

    #k_mode[k_mode==0] = np.inf
    #k_mode_2[k_mode_2==0] = np.inf
    normal[normal==0] = np.inf
    normal_2[normal_2==0] = np.inf
    power_1d = power/normal
    power_1d_err = np.sqrt(err/normal_2)
    #power_1d = power/k_mode
    #power_1d_err = np.sqrt(err/k_mode_2)

    #normal[normal==0] = np.inf
    #normal_2[normal_2==0] = np.inf
    #power_1d = power/normal/k_mode
    #power_1d_err = np.sqrt(err/normal_2)/k_mode

    return power_1d, power_1d_err, k_p_centre
