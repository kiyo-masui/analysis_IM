#! /usr/bin/env python 

import numpy as np
import ps_summary
import matplotlib.pyplot as plt

def chisquare_estimation(ps_power_list, ps_error_list, ps_kbins_list,
                         redshift=0.8, cross=False):

    noise_vector = np.array(ps_error_list).flatten()[:, None]
    power_vector = np.array(ps_power_list).flatten()[:, None]
    kbins_vector = np.array(ps_kbins_list[0])[:,None]

    #power_vector = np.fabs(power_vector)
    nonzero = power_vector != 0
    noninf = np.logical_and(np.isfinite(power_vector), np.isfinite(noise_vector))
    power_vector = power_vector[np.logical_and(nonzero, noninf)][:,None]
    noise_vector = noise_vector[np.logical_and(nonzero, noninf)][:,None]
    kbins_vector = kbins_vector[np.logical_and(nonzero, noninf)]

    ps_th = ps_summary.load_theory_ps(kbins_vector, 
                                      redshift=redshift, 
                                      cross=cross)
    ps_th = ps_th[None, :]


    N = np.dot(noise_vector, noise_vector.T) * np.eye(noise_vector.shape[0])
    W = np.linalg.inv(N)
    A = np.repeat(ps_th, len(ps_power_list), axis=0).reshape(-1,1)

    P = np.linalg.inv(np.dot(np.dot(A.T, W),A))
    Q = np.dot(A.T, W)
    bias = np.dot(np.dot(P, Q), power_vector)

    bias_error = np.sqrt(np.dot(np.dot(P, np.dot(np.dot(A.T, W), A)), P.T))

    if not cross:
        bias = np.sqrt(bias)
        bias_error = np.sqrt(bias_error)

    return bias[0,0], bias_error[0,0]

def chisquare_fitting(ps_power_list, ps_error_list, ps_kbins_list, 
                      bias=[0.1, 1.5, 100], redshift=0.8, cross=False, plot=False):
    '''
        bias : [min, max, number]
    '''
    # get the theoretial power spectrum  (unitless)

    b = np.linspace( bias[0], bias[1], bias[2] )

    chi_sq_list = []
    data_number = 0.

    for i in range(len(ps_power_list)):

        ps_power = ps_power_list[i]
        ps_error = ps_error_list[i]
        ps_kbins = ps_kbins_list[i]

        ps_th = ps_summary.load_theory_ps(ps_kbins, redshift=redshift, cross=cross)

        ps_error[ps_error==0] = np.inf
        ps_power = np.fabs(ps_power)

        data_number += len(ps_error[ps_error!=0])

        chi_sq  = (b[:, None]*ps_th[None, :] - ps_power[None, :])**2.
        chi_sq /= ps_error[None, :]**2.
        chi_sq  = np.sum(chi_sq, axis=1)
        chi_sq_list.append(chi_sq)

    chi_sq = np.sum(np.array(chi_sq_list), axis=0)

    L = np.exp(-0.5 * chi_sq)

    L_argmax = L.argmax()
    L_step = L[1] - L[0]
    # searching upper limit
    p_upper_arg = L_argmax
    p_upper_total = L[L_argmax:-1].sum() * L_step - L[L_argmax] * 0.5 * L_step
    p_upper = L[L_argmax] * 0.5 * L_step
    for i in range(L_argmax + 1, len(L)):
        p_upper += L[i] * L_step
        if p_upper/p_upper_total > 0.68:
            p_upper_arg = i
            break
    # searching lower limit
    p_lower_arg = L_argmax
    p_lower_total = L[0:L_argmax+1].sum() * L_step - L[L_argmax] * 0.5 * L_step
    p_lower = L[L_argmax] * 0.5 * L_step
    for i in range(L_argmax - 1, 0, -1):
        p_lower += L[i] * L_step
        if p_lower/p_lower_total > 0.68:
            p_lower_arg = i
            break

    if plot:
        plt.figure(figsize=(8,8))
        plt.axvspan(xmin=b[p_lower_arg], xmax=b[p_upper_arg], 
                    ymin=0, alpha=0.2)
        plt.step(b, L, where='mid')
        plt.ylim(ymin=0, ymax=L.max())
        plt.savefig('./png/bias_fitting.png')
        #plt.show()


    return b[chi_sq.argmin()],\
           b[p_lower_arg]/np.sqrt(data_number),\
           b[p_upper_arg]/np.sqrt(data_number)
