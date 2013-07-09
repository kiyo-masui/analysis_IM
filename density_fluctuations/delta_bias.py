import growth_factor as gf
import camb_caller as cc
import numpy as np
import numdifftools as nd

def bias(z):
    """
    Calculate bias
    """
    return 2.

def pwrspec_f_nl(vect, z = 0.5):
    """
    Read cosmological parameters at vect and calculate powerspectrum,
    accounting for the f_NL parameter of the primordial perturbations.
    """
    k           = vect[0]
    f_NL        = vect[1]
    hubble_0    = 70
    omega_m     = 0.27
    b           = bias(z)
    delta_c     = 1.686
    c           = 300000.
    params_dict = {'hubble': hubble_0, 'transfer_redshift(1)': z}
    P,T = cc.camb_pipe(params_dict)
    delta_b     = 10**T((-4))*(1+z)*3*f_NL*(b - 1)*omega_m*hubble_0**2*delta_c/ \
                  (c**2*10**T(k)*(10**k)**2*gf.g(1/(1+z)))
    print "\n paramters are:", [k,f_NL], "\n"
    return 10**P(k)*(b + delta_b)**2

def pwrspec_deriv(vect):
    """
    Calculate gradient of the powerspectrum function.
    """
    deriv_pwr = nd.Gradient(pwrspec_f_nl, stepMax=0.1)
    return deriv_pwr(vect)

