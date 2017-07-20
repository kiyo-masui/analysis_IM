#! /usr/bin/env python 

import numpy as np
from core import algebra
from simulations import corr21cm
import os

def load_2df_ps_from_paper(name='Percival'):
    data = np.loadtxt('/Users/ycli/Code/analysis_IM/simulations/data/ps_2df_%s.dat'%name)
    data[:,1] = data[:,1]*data[:,0]**3/2./np.pi**2.
    data[:,2] = data[:,2]*data[:,0]**3/2./np.pi**2.

    return data[:,0], data[:,1], data[:,2]

def load_camb():
    data = np.loadtxt('/Users/ycli/Code/analysis_IM/mkpower/cambini/_matterpower.dat')
    data[:,0] *= 0.72
    data[:,1] /= 0.72**3

    return data
    #return data[:,0], data[:,1]

def load_theory_ps(k_bins_centre, redshift=0.08, cross=False, in_Jy=False):
    c21 = corr21cm.Corr21cm(redshift=redshift, nu_upper=1400., nu_lower=1200.)
    T_b = c21.T_b(redshift) * 1.e-3
    if in_Jy:
        jy = 2.*1.38e-23/1.e-26/((1.+redshift)*0.21)**2
    else:
        jy = 1.

    if cross:
        ps_theory  = jy*c21.get_pwrspec(k_bins_centre)/T_b
    else:
        ps_theory  = jy**2.*c21.get_pwrspec(k_bins_centre)


    ps_theory *= k_bins_centre**3./2./(np.pi**2.)

    #print ps_theory

    return ps_theory

    #power_th = np.loadtxt('/Users/ycli/Code/analysis_IM/simulations/data/wigglez_halofit_z0.8.dat')
    #ps_theory *= power_th[:,0]**3/2./np.pi**2
    #return ps_theory, power_th[:,0]

def get_1d_k_bin_centre(ps_1d):

    k = np.arange(ps_1d.shape[-1]) - ps_1d.shape[-1]//2
    k = ps_1d.info['k_centre'] * ps_1d.info['k_delta']**k

    return k

def get_2d_k_bin_centre(ps_2d):

    k_p = np.arange(ps_2d.shape[0]) - ps_2d.shape[0]//2
    k_p = ps_2d.info['k_p_centre'] * ps_2d.info['k_p_delta']**k_p

    k_v = np.arange(ps_2d.shape[1]) - ps_2d.shape[1]//2
    k_v = ps_2d.info['k_v_centre'] * ps_2d.info['k_v_delta']**k_v

    return k_p, k_v

def get_1d_k_bin_edges(ps_1d):
    
    k = get_1d_k_bin_centre(ps_1d)

    k_edges  = np.append(k, k[-1]*ps_1d.info['k_delta'])
    k_edges /= np.sqrt(ps_1d.info['k_delta'])

    return k_edges

def get_2d_k_bin_edges(ps_2d):

    k_p, k_v = get_2d_k_bin_centre(ps_2d)

    k_p_edges  = np.append(k_p, k_p[-1]*ps_2d.info['k_p_delta'])
    k_p_edges /= np.sqrt(ps_2d.info['k_p_delta'])

    k_v_edges  = np.append(k_v, k_v[-1]*ps_2d.info['k_v_delta'])
    k_v_edges /= np.sqrt(ps_2d.info['k_v_delta'])

    return k_p_edges, k_v_edges

def load_transfer_function(rf_root, tr_root, h5obj=None, auto=True):
    
    if h5obj==None:
        ps_rf = algebra.make_vect(algebra.load(rf_root))
        kn_rf = algebra.make_vect(algebra.load(rf_root.replace('pow', 'kmn')))
        ps_tr = algebra.make_vect(algebra.load(tr_root))
        kn_tr = algebra.make_vect(algebra.load(tr_root.replace('pow', 'kmn')))
    else:
        ps_rf = algebra.make_vect(algebra.load_h5(h5obj, rf_root))
        kn_rf = algebra.make_vect(algebra.load_h5(h5obj, rf_root.replace('pow', 'kmn')))
        ps_tr = algebra.make_vect(algebra.load_h5(h5obj, tr_root))
        kn_tr = algebra.make_vect(algebra.load_h5(h5obj, tr_root.replace('pow', 'kmn')))

    kn = kn_rf * kn_tr

    ps_rf[ps_rf==0] = np.inf
    transfer_function = ps_tr/ps_rf

    transfer_function[kn == 0] = 0.
    transfer_function[transfer_function<0.] = 0.
    #transfer_function[transfer_function>1.] = 1.

    transfer_function[transfer_function == 0.] = np.inf
    #transfer_function[transfer_function > 1.] = 1.

    transfer_function = 1./transfer_function

    if auto:
        transfer_function *= transfer_function

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_rf)

    return transfer_function, k_p_edges, k_v_edges

def load_weight_sec(ns_root, transfer_function=None):
    sec_pair_list = ['_AB', '_BA', '_CA', '_DA']

    noise = []
    for sec_pair in sec_pair_list:
        ps_ns = algebra.make_vect(algebra.load(ns_root + sec_pair))
        noise.append(ps_ns)
    noise = np.array(noise)
    noise = np.mean(noise, axis=0)

    noise[noise==0] = np.inf

    weight = []
    #for sec_pair in sec_pair_list:
    for sec_pair in ['_AB',]:
        kn_ns = algebra.make_vect(algebra.load(ns_root.replace('pow', 'kmn') + sec_pair))
        weight_sec = np.abs(kn_ns)/noise**2
        if transfer_function!=None:
            weight_sec /= transfer_function**2
        weight_sec[np.isnan(weight_sec)] = 0.
        weight_sec[np.isinf(weight_sec)] = 0.
        weight_sec[kn_ns == 0] = 0.
        weight.append(weight_sec)

    weight = np.array(weight)
    weight = np.mean(weight, axis=0)

    return weight

def load_weight(ns_root, transfer_function=None):
    
    ps_ns = algebra.make_vect(algebra.load(ns_root))
    kn_ns = algebra.make_vect(algebra.load(ns_root.replace('pow', 'kmn')))


    weight = kn_ns/ps_ns**2
    if transfer_function != None:
        weight /= transfer_function**2

    weight[ np.isnan(weight) ] = 0.
    weight[ np.isinf(weight) ] = 0.
    weight[ kn_ns == 0. ] = 0.

    #kn_ns[kn_ns==0] = np.inf
    #gauss_noise = 2. * ps_ns / np.sqrt(kn_ns)
    #if transfer_function != None:
    #    gauss_noise *= transfer_function

    #gauss_noise[gauss_noise==0] = np.inf
    #weight = (1./gauss_noise)**2

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_ns)

    return weight, k_p_edges, k_v_edges

def load_noise_err(ne_root):

    noise_error = algebra.make_vect(algebra.load(ne_root))

    k_p_edges, k_v_edges = get_2d_k_bin_edges(noise_error)

    return noise_error, k_p_edges, k_v_edges


def load_power_spectrum_1d(ps_root, raw=False):
    
    ps_1d = algebra.make_vect(algebra.load(ps_root))

    if not raw:
        ps_1derr = algebra.make_vect(algebra.load(sn_root.replace('pow', 'err')))
        ps_1dkmn = algebra.make_vect(algebra.load(ps_root.replace('pow', 'kmn')))
        k_centre = get_1d_k_bin_centre(ps_1d)
        return ps_1d, ps_1derr, ps_1dkmn, k_centre
    else:
        k_centre = get_1d_k_bin_centre(ps_1d)
        return ps_1d, k_centre
        

def load_power_spectrum_opt_1d(ps_root, sn_root):

    ps_1d = algebra.make_vect(algebra.load(ps_root))
    sn_1d = algebra.make_vect(algebra.load(sn_root))
    ps_1derr = algebra.make_vect(algebra.load(sn_root.replace('pow', 'err')))
    ps_1dkmn = algebra.make_vect(algebra.load(ps_root.replace('pow', 'kmn')))

    k_edges = get_1d_k_bin_edges(ps_1d)

    ps_1d -= sn_1d

    return ps_1d, ps_1derr, ps_1dkmn, k_edges

def load_power_spectrum_opt(ps_root, sn_root):

    ps_2d = algebra.make_vect(algebra.load(ps_root))
    sn_2d = algebra.make_vect(algebra.load(sn_root))
    ps_2derr = algebra.make_vect(algebra.load(sn_root.replace('pow', 'err')))
    ps_2dkmn = algebra.make_vect(algebra.load(ps_root.replace('pow', 'kmn')))

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2d)

    # subtract short noise
    ps_2d -= sn_2d

    return ps_2d, ps_2derr, ps_2dkmn, k_p_edges, k_v_edges, sn_2d

def load_power_spectrum(ps_root, h5obj=None):

    if h5obj == None:
        ps_2d = algebra.make_vect(algebra.load(ps_root))
        kn_2d = algebra.make_vect(algebra.load(ps_root.replace('pow', 'kmn')))
    else:
        ps_2d = algebra.make_vect(algebra.load_h5(h5obj, ps_root))
        kn_2d = algebra.make_vect(algebra.load_h5(h5obj, ps_root.replace('pow', 'kmn')))

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2d)

    return ps_2d, kn_2d, k_p_edges, k_v_edges

def load_power_spectrum_err(ps_root):
    
    ps_2derr = algebra.make_vect(algebra.load(ps_root.replace('pow', 'err')))

    k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2derr)

    return ps_2derr, k_p_edges, k_v_edges

def load_power_spectrum_sec(ps_root, ne_root, sec=['A', 'B', 'C', 'D'], 
                            transfer_function=None):
    
    power_list = []
    k_num_list = []
    weight_list = []

    k_p_edges = []
    k_v_edges = []

    for i in range(len(sec)):
        for j in range(i+1, len(sec)):
            tab_l = sec[i] + sec[j]
            weight_l = load_weight(ne_root%tab_l, transfer_function)[0]
            tab_r = sec[j] + sec[i]
            weight_r = load_weight(ne_root%tab_r, transfer_function)[0]

            weight = weight_l * weight_r
            weight_list.append(weight)

            tab = tab_l + 'x' + tab_r
            ps_2d = algebra.make_vect(algebra.load(ps_root%tab))
            kn_2d = algebra.make_vect(algebra.load((ps_root%tab).replace('pow','kmn')))

            if k_p_edges == []:
                k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2d)

            ps_2d *= weight
            kn_2d *= weight

            power_list.append(ps_2d)
            k_num_list.append(kn_2d)

    power = np.sum(np.array(power_list), axis=0)
    k_num = np.sum(np.array(k_num_list), axis=0)
    power_err = np.sum((power[None, ...] - np.array(power_list))**2., axis=0)
    weight = np.sum(np.array(weight_list), axis=0)
    weight_square = np.sum(np.array(weight_list)**2, axis=0)
    weight_mean = np.mean(np.array(weight_list), axis=0)

    weight[weight==0] = np.inf
    power /= weight
    k_num /= weight
    power_err /= weight_square
    power_err = np.sqrt(power_err)

    power = algebra.make_vect(power, axis_names=power_list[0].info['axes'])
    power.info = power_list[0].info
    k_num = algebra.make_vect(k_num, axis_names=k_num_list[0].info['axes'])
    power.info = k_num_list[0].info

    return power, power_err, k_num, k_p_edges, k_v_edges, weight_mean

def convert_2dps_to_1dps_sim(rf_root, ns_root=None, truncate_range=None ):

    power_spectrum = algebra.make_vect(algebra.load(rf_root))
    power_spectrum_err = algebra.make_vect(algebra.load(rf_root.replace('pow', 'err')))
    k_mode_number = algebra.make_vect(algebra.load(rf_root.replace('pow', 'kmn')))
    if ns_root != None:
        weight, k_p_edges, k_v_edges = load_weight(ns_root)
        #weight = np.ones(power_spectrum.shape)
        #k_mode_number *= weight
        #k_mode_number = weight
    else:
        weight = np.ones(power_spectrum.shape)

    k_p_edges, k_v_edges = get_2d_k_bin_edges(power_spectrum)
    k_p_centre, k_v_centre = get_2d_k_bin_centre(power_spectrum)
    k_centre = np.sqrt(k_p_centre[:,None]**2 + k_v_centre[None, :]**2)
    k_edges_1d = k_p_edges

    power_spectrum *= weight
    power_spectrum_err *= weight

    power_spectrum_err **= 2

    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        power_spectrum_err, kpe, kve=truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                          k_v_edges, power_spectrum_err)
        power_spectrum, kpe, kve    =truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                          k_v_edges, power_spectrum)
        weight, kpe, kve            =truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                          k_v_edges, weight)
        k_centre, kpe, kve          =truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                          k_v_edges, k_centre)
        k_mode_number, kpe, kve     =truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                          k_v_edges, k_mode_number)

    k_mode_number = k_mode_number.flatten()
    power_spectrum_err = power_spectrum_err.flatten()#[k_mode_number!=0.]
    power_spectrum = power_spectrum.flatten()#[k_mode_number!=0.]
    weight = weight.flatten()#[k_mode_number!=0]
    k_centre = k_centre.flatten()#[k_mode_number!=0]

    normal,k_e = np.histogram(k_centre, k_edges_1d, weights=weight)
    normal_2,k_e = np.histogram(k_centre, k_edges_1d, weights=weight**2)
    power, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum)
    err,   k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum_err)
    normal[normal==0] = np.inf
    normal_2[normal_2==0] = np.inf
    power_1d = power/normal
    power_1d_err = np.sqrt(err/normal_2)
    return power_1d, power_1d_err, k_p_centre

def truncate_2dps(k_p_range, k_v_range, k_p_edges, k_v_edges, power_spectrum):
    
    k_p_range_idx = np.digitize(k_p_range, k_p_edges) - 1
    k_v_range_idx = np.digitize(k_v_range, k_v_edges) - 1

    k_p_range_idx[k_p_range_idx<0] = 0
    k_v_range_idx[k_v_range_idx<0] = 0
    #print k_p_range_idx
    #print k_v_range_idx

    #return power_spectrum[k_p_range_idx[0]:k_p_range_idx[1],
    #                      k_v_range_idx[0]:k_v_range_idx[1]],\
    #       k_p_edges[k_p_range_idx[0]:k_p_range_idx[1]+1],\
    #       k_v_edges[k_v_range_idx[0]:k_v_range_idx[1]+1]
    return power_spectrum[k_v_range_idx[0]:k_v_range_idx[1],
                          k_p_range_idx[0]:k_p_range_idx[1]],\
           k_p_edges[k_p_range_idx[0]:k_p_range_idx[1]+1],\
           k_v_edges[k_v_range_idx[0]:k_v_range_idx[1]+1]

def convert_2dps_to_1dps_opt(ps_root, sn_root, truncate_range=None):

    power_spectrum, power_spectrum_err, power_spectrum_kmn,\
    k_p_edges, k_v_edges, shortnoise = load_power_spectrum_opt(ps_root, sn_root)

    power_spectrum_err **= 2

    k_p_centre, k_v_centre = get_2d_k_bin_centre(power_spectrum)
    k_centre = np.sqrt(k_v_centre[:,None]**2 + k_p_centre[None, :]**2)

    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        power_spectrum_kmn, kpe, kve = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, power_spectrum_kmn)
        power_spectrum_err, kpe, kve = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, power_spectrum_err)
        power_spectrum, kpe, kve     = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, power_spectrum)
        k_centre, kpe, kve           = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, k_centre)
        shortnoise, kpe, kve         = truncate_2dps( k_p_range, k_v_range, k_p_edges, 
                                            k_v_edges, shortnoise)
    power_spectrum_kmn = power_spectrum_kmn.flatten()
    power_spectrum_err = power_spectrum_err.flatten()[power_spectrum_kmn!=0]
    power_spectrum     = power_spectrum.flatten()[power_spectrum_kmn!=0]
    k_centre           = k_centre.flatten()[power_spectrum_kmn!=0]
    shortnoise         = shortnoise.flatten()[power_spectrum_kmn!=0]
    power_spectrum_kmn = power_spectrum_kmn[power_spectrum_kmn!=0]

    #power_spectrum *= power_spectrum_kmn
    #power_spectrum_err *= power_spectrum_kmn**2
    #shortnoise *= power_spectrum_kmn

    k_edges_1d = k_p_edges

    #normal, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum_kmn)
    #normal_2, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum_kmn**2)
    normal, k_e = np.histogram(k_centre, k_edges_1d)
    power, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum)
    err, k_e = np.histogram(k_centre, k_edges_1d, weights=power_spectrum_err)
    short, k_e = np.histogram(k_centre, k_edges_1d, weights=shortnoise)

    normal = normal.astype(float)
    normal[normal==0] = np.inf
    #normal_2[normal_2==0] = np.inf
    power_1d = power/normal
    #power_1d_err = np.sqrt(err/normal_2)
    power_1d_err = np.sqrt(err)/normal
    shortnoise_1d = short/normal

    return power_1d, power_1d_err, k_p_centre, shortnoise_1d

def convert_2dps_to_1dps_each(ps_list, kn_list, weight, truncate_range=None, order=0):

    ps_1d_list = []
    k_centre = None

    if not isinstance(ps_list, list):
        ps_list = [ps_list, ]
        kn_list = [kn_list, ]

    for i in range(len(ps_list)):
        ps_2d = ps_list[i]
        kn_2d = kn_list[i]
        if weight != None:
            wt_2d = weight
        else:
            wt_2d = np.ones(ps_2d.shape)
            wt_2d[ps_2d==0] = 0.

        k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2d)
        k_p_centr, k_v_centr = get_2d_k_bin_centre(ps_2d)

        k_edges = k_p_edges
        k_centre = k_p_centr

        kc_2d = np.sqrt(k_v_centr[:,None]**2 + k_p_centr[None, :]**2)
        mu_2d = k_p_centr[None, :] / kc_2d

        # for dipole
        coef = [0, 0, 0]
        coef[order] = 1
        P = np.polynomial.legendre.Legendre(coef)
        #ps_2d *= (mu_2d ** order)
        ps_2d = ps_2d * P(mu_2d)
        #wt_2d *= P(mu_2d)

        if truncate_range != None:
            k_p_range = [truncate_range[0], truncate_range[1]]
            k_v_range = [truncate_range[2], truncate_range[3]]

            restrict_v = np.where(np.logical_or(k_v_centr<k_v_range[0], 
                                                k_v_centr>k_v_range[1]))
            restrict_p = np.where(np.logical_or(k_p_centr<k_p_range[0], 
                                                k_p_centr>k_p_range[1]))
            wt_2d[restrict_v, :] = 0.
            wt_2d[:, restrict_p] = 0.


        #Commented only for min test
        ps_2d = ps_2d * wt_2d
        kn_2d = kn_2d * wt_2d

        ps_2d[kn_2d==0] = 0.
        wt_2d[kn_2d==0] = 0.

        kc_2d = kc_2d.flatten()
        ps_2d = ps_2d.flatten()
        wt_2d = wt_2d.flatten()
        kn_2d = kn_2d.flatten()

        wt_1d, k_e = np.histogram(kc_2d, k_edges, weights=wt_2d)
        ps_1d, k_e = np.histogram(kc_2d, k_edges, weights=ps_2d)
        kn_1d, k_e = np.histogram(kc_2d, k_edges, weights=kn_2d)
        kc_2d = np.array(kc_2d)        
        #try:
        #    ps_1d[0] = np.median(ps_2d[np.where((kc_2d <= k_edges[0]) & (kc_2d > 0))])
        #    print ps_1d[0]
        #except ValueError:
        #    pass
        #for i in range(len(k_edges)-2):
        #    try:
        #        ps_1d[i+1] = np.median(ps_2d[np.where((kc_2d <= k_edges[i+1]) & (kc_2d >  k_edges[i]))])
        #        print ps_1d[i+1]
        #    except ValueError:
        #        pass

        wt_1d[wt_1d==0] = np.inf
        ps_1d /= wt_1d
       
        ps_1d_list.append(ps_1d)

    ps_1d_list = np.ma.array(ps_1d_list)
    ps_1d_list[np.isnan(ps_1d_list)] = np.ma.masked
    ps_1d_list[np.isinf(ps_1d_list)] = np.ma.masked

    ps_1d_std = np.ma.std(ps_1d_list, axis=0)
    ps_1d_mean = np.mean(ps_1d_list, axis=0)
    #ps_1d_std = 1 / np.sqrt(wt_1d)
    return ps_1d_mean, ps_1d_std, k_centre


def convert_2dps_to_1dps_sec(ps_root, ns_root, tr_root, rf_root, 
                             sec=['A', 'B', 'C', 'D'], truncate_range=None):

    ps_1d_list = []
    k_centre = None

    if os.path.exists(rf_root) and os.path.exists(tr_root):
        transfer_function = load_transfer_function(rf_root, tr_root)[0]
    else:
        print "Note: data for transfer function estimation not exists."
        transfer_function = None

    #noise  = algebra.make_vect(algebra.load(ns_root.replace('_%s', '')))
    #noise  = load_weight(ns_root.replace('_%s', ''), transfer_function)[0]
    noise  = load_weight_sec(ns_root.replace('_%s', ''), transfer_function)
    #signal = algebra.make_vect(algebra.load(ps_root.replace('_%s', '')))
    #weight = 1./(noise - signal)**2
    weight = noise

    for i in range(len(sec)):
        for j in range(i+1, len(sec)):
            tab_l = sec[i] + sec[j]
            weight_l = load_weight(ns_root%tab_l, transfer_function)[0]
            tab_r = sec[j] + sec[i]
            weight_r = load_weight(ns_root%tab_r, transfer_function)[0]

            #wt_2d = weight_l * weight_r
            #wt_2d = np.ones(wt_2d.shape)
            wt_2d = weight

            tab = tab_l + 'x' + tab_r
            ps_2d = algebra.make_vect(algebra.load(ps_root%tab))
            kn_2d = algebra.make_vect(algebra.load((ps_root%tab).replace('pow','kmn')))

            if transfer_function != None:
                ps_2d *= transfer_function
                kn_2d /= transfer_function**2

            k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2d)
            k_p_centr, k_v_centr = get_2d_k_bin_centre(ps_2d)

            k_edges = k_p_edges
            k_centre = k_p_centr

            kc_2d = np.sqrt(k_v_centr[:,None]**2 + k_p_centr[None, :]**2)

            if truncate_range != None:
                k_p_range = [truncate_range[0], truncate_range[1]]
                k_v_range = [truncate_range[2], truncate_range[3]]
                #ps_2d = truncate_2dps(k_p_range,k_v_range,k_p_edges,k_v_edges,ps_2d)[0]
                #kc_2d = truncate_2dps(k_p_range,k_v_range,k_p_edges,k_v_edges,kc_2d)[0]
                #wt_2d = truncate_2dps(k_p_range,k_v_range,k_p_edges,k_v_edges,wt_2d)[0]
                #kn_2d = truncate_2dps(k_p_range,k_v_range,k_p_edges,k_v_edges,kn_2d)[0]
                restrict_v = np.where(np.logical_or(k_v_centr<k_v_range[0], 
                                                    k_v_centr>k_v_range[1]))
                restrict_p = np.where(np.logical_or(k_p_centr<k_p_range[0], 
                                                    k_p_centr>k_p_range[1]))
                wt_2d[restrict_v, :] = 0.
                wt_2d[:, restrict_p] = 0.

            ps_2d *= wt_2d
            kn_2d *= wt_2d

            ps_2d[kn_2d==0] = 0.
            wt_2d[kn_2d==0] = 0.

            kc_2d = kc_2d.flatten()
            ps_2d = ps_2d.flatten()
            wt_2d = wt_2d.flatten()
            kn_2d = kn_2d.flatten()

            wt_1d, k_e = np.histogram(kc_2d, k_edges, weights=wt_2d)
            ps_1d, k_e = np.histogram(kc_2d, k_edges, weights=ps_2d)
            kn_1d, k_e = np.histogram(kc_2d, k_edges, weights=kn_2d)

            ps_1d /= wt_1d

            ps_1d_list.append(ps_1d)

    ps_1d_list = np.ma.array(ps_1d_list)
    ps_1d_list[np.isnan(ps_1d_list)] = np.ma.masked
    ps_1d_list[np.isinf(ps_1d_list)] = np.ma.masked

    ps_1d_mean = np.ma.mean(ps_1d_list, axis=0)
    ps_1d_std = np.ma.std(ps_1d_list, axis=0)

    return ps_1d_mean, ps_1d_std, k_centre

def convert_2dps_to_1dps(ps_root, ns_root, tr_root, rf_root, ne_root=None, 
                         sec=['A', 'B', 'C', 'D'], truncate_range=None):

    if os.path.exists(rf_root):
        transfer_function = load_transfer_function(rf_root, tr_root)[0]
    else:
        print "Reference file or Transfer file are not esits, ignore compensation"
        transfer_function = None

    if ne_root != None:
        power_spectrum, power_spectrum_err, k_mode_number, k_p_edges, k_v_edges,\
            weight = load_power_spectrum_sec(ps_root, ne_root, sec,transfer_function)
    else: 
        weight, k_p_edges, k_v_edges = load_weight(ns_root, transfer_function)
        power_spectrum, k_mode_number, k_p_edges, k_v_edges =\
            load_power_spectrum(ps_root)
        power_spectrum_err, k_p_edges, k_v_edges = load_power_spectrum_err(ps_root)

    if transfer_function != None:
        power_spectrum *= transfer_function
        power_spectrum_err *= transfer_function

    power_spectrum *= weight
    power_spectrum_err *= weight

    power_spectrum_err **= 2.

    k_p_centre, k_v_centre = get_2d_k_bin_centre(power_spectrum)
    #k_centre = np.sqrt(k_p_centre[:,None]**2 + k_v_centre[None, :]**2)
    k_centre = np.sqrt(k_v_centre[:,None]**2 + k_p_centre[None, :]**2)

    if truncate_range != None:
        k_p_range = [truncate_range[0], truncate_range[1]]
        k_v_range = [truncate_range[2], truncate_range[3]]
        power_spectrum_err, kpe, kve = truncate_2dps(k_p_range,k_v_range,k_p_edges, 
                                            k_v_edges, power_spectrum_err)
        power_spectrum, kpe, kve     = truncate_2dps(k_p_range,k_v_range,k_p_edges, 
                                            k_v_edges, power_spectrum)
        k_centre, kpe, kve           = truncate_2dps(k_p_range,k_v_range,k_p_edges, 
                                            k_v_edges, k_centre)
        weight, kpe, kve             = truncate_2dps(k_p_range,k_v_range,k_p_edges, 
                                            k_v_edges, weight)
        k_mode_number, kpe, kve      = truncate_2dps(k_p_range,k_v_range,k_p_edges, 
                                                      k_v_edges, k_mode_number)

    #print power_spectrum.shape
    #k_mode_number      = k_mode_number.flatten()
    #power_spectrum_err = power_spectrum_err.flatten()[k_mode_number!=0]
    #power_spectrum     = power_spectrum.flatten()[k_mode_number!=0]
    #weight             = weight.flatten()[k_mode_number!=0]
    #k_centre           = k_centre.flatten()[k_mode_number!=0]

    #power_spectrum_err = power_spectrum_err[weight!=0]
    #power_spectrum     = power_spectrum[weight!=0]
    #k_centre           = k_centre[weight!=0]
    #weight             = weight[weight!=0]

    power_spectrum_err[ k_mode_number == 0 ] = 0.
    power_spectrum[ k_mode_number == 0 ] = 0.
    weight[ k_mode_number == 0 ] = 0.

    power_spectrum_err = power_spectrum_err.flatten()
    power_spectrum = power_spectrum.flatten()
    weight = weight.flatten()
    k_centre = k_centre.flatten()

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

def covar_1d (ps_list, kn_list, weight, truncate_range=None, order=0):

    ps_1d_list = []
    k_centre = None

    if not isinstance(ps_list, list):
        ps_list = [ps_list, ]
        kn_list = [kn_list, ]

    for i in range(len(ps_list)):
        ps_2d = ps_list[i]
        kn_2d = kn_list[i]
        if weight != None:
            wt_2d = weight
        else:
            wt_2d = np.ones(ps_2d.shape)
            wt_2d[ps_2d==0] = 0.

        k_p_edges, k_v_edges = get_2d_k_bin_edges(ps_2d)
        k_p_centr, k_v_centr = get_2d_k_bin_centre(ps_2d)

        k_edges = k_p_edges
        k_centre = k_p_centr

        kc_2d = np.sqrt(k_v_centr[:,None]**2 + k_p_centr[None, :]**2)
        mu_2d = k_p_centr[None, :] / kc_2d

        # for dipole
        coef = [0, 0, 0]
        coef[order] = 1
        P = np.polynomial.legendre.Legendre(coef)
        #ps_2d *= (mu_2d ** order)
        ps_2d = ps_2d * P(mu_2d)
        #wt_2d *= P(mu_2d)

        if truncate_range != None:
            k_p_range = [truncate_range[0], truncate_range[1]]
            k_v_range = [truncate_range[2], truncate_range[3]]

            restrict_v = np.where(np.logical_or(k_v_centr<k_v_range[0],
                                                k_v_centr>k_v_range[1]))
            restrict_p = np.where(np.logical_or(k_p_centr<k_p_range[0],
                                                k_p_centr>k_p_range[1]))
            wt_2d[restrict_v, :] = 0.
            wt_2d[:, restrict_p] = 0.

        ps_2d = ps_2d * wt_2d
        kn_2d = kn_2d * wt_2d

        ps_2d[kn_2d==0] = 0.
        wt_2d[kn_2d==0] = 0.

        kc_2d = kc_2d.flatten()
        ps_2d = ps_2d.flatten()
        wt_2d = wt_2d.flatten()
        kn_2d = kn_2d.flatten()

        wt_1d, k_e = np.histogram(kc_2d, k_edges, weights=wt_2d)
        ps_1d, k_e = np.histogram(kc_2d, k_edges, weights=ps_2d)
        kn_1d, k_e = np.histogram(kc_2d, k_edges, weights=kn_2d)

        wt_1d[wt_1d==0] = np.inf
        ps_1d /= wt_1d

        ps_1d_list.append(ps_1d)

    ps_1d_list = np.ma.array(ps_1d_list)
    ps_1d_list[np.isnan(ps_1d_list)] = np.ma.masked
    ps_1d_list[np.isinf(ps_1d_list)] = np.ma.masked

    ps_1d_corr = np.corrcoef(ps_1d_list.T)
    return ps_1d_corr, k_centre, k_edges

