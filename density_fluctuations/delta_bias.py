import growth_factor as gf
import camb_caller as cc
import numpy as np
import numdifftools as nd
from utils import batch_handler

def bias(z):
    """
    Calculate bias
    """
    return 2.

@batch_handler.memoized
def pwrspec_f_nl(vect, z = 0.5):
    """
    Read cosmological parameters at vect and calculate powerspectrum,
    accounting for the f_NL parameter of the primordial perturbations.
    """
    print vect
    k           = vect[0]
    f_NL        = vect[1]
    hubble_0    = 70
    omega_m     = 0.27
    b           = bias(z)
    delta_c     = 1.686
    c           = 300000.
    params_dict = {'hubble': hubble_0, 'transfer_redshift(1)': z}
    P,T         = cc.camb_pipe(params_dict)
    # put T_k back
    #T_k         = 10**T(-4)/10**T(k)
    delta_b     = (1+z)*3*f_NL*(b - 1)*omega_m*hubble_0**2*delta_c/ \
                  (c**2*(10**k)**2*gf.g(1/(1+z)))
    return 10**P(k)*(b + delta_b)**2

def pwrspec_deriv(vect):
    """
    Calculate gradient of the powerspectrum function.
    """
    deriv_pwr = nd.Gradient(pwrspec_f_nl, step_max=0.005)
    return deriv_pwr(vect)

def pwrspec_grad(vect, step_factor = 0.001):
    num      = len(vect)
    gradient = np.zeros(num)
    for ind in range(num):
        plus_vect       = np.array(vect)
        plus_vect[ind]  = np.array(vect)[ind]*(1.+step_factor)
        minus_vect      = np.array(vect)
        minus_vect[ind] = np.array(vect)[ind]*(1.-step_factor)
        gradient[ind]   = (pwrspec_f_nl(plus_vect) - pwrspec_f_nl(minus_vect))/(2*np.array(vect)[ind]*step_factor)
    return gradient
