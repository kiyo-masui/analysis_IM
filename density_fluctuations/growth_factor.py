import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Set cosmological parameters
Omega_m = 0.27
Omega_lambda = 0.73
Omega_r = 0.
#Omega_k = 1. - (Omega_m + Omega_lambda + Omega_r)
Omega_k = 0.
Hubble_0 = 70

def hubble(a):
    """
    Calculate hubble constant for the scalefactor a.
    """
    return Hubble_0*np.sqrt(Omega_r*a**(-4) + Omega_m*a**(-3) +
                         Omega_k*a**(-2) + Omega_lambda)

# Part for the linear growth factor

def integration_f(a):
    """
    Integrand
    """
    return (a*hubble(a)/Hubble_0)**(-3)

def g(a):
    """
    Linear growth factor
    """
    return (5*Omega_m*hubble(a)/(a*2*Hubble_0))*\
                  quad(integration_f,0,a)[0]

def dg_dlna(array, delta_a):
    """
    Derivative with respecto to natural log scale factor
    """
    return array*(g(array + delta_a) - g(array))/(delta_a)

# Analytic expression for check from Dodelson
def var(a):
    return (a*(1-Omega_m))/Omega_m
def D(array):
    b = np.zeros_like(array)
    for num in range(len(array)):
        x = array[num]
        b[num] = 5*Omega_m/(2*(1-Omega_m))*(3*((1+var(x))**(0.5))/(var(x)**(1.5))*np.log((1+var(x))**0.5-var(x)**0.5) + 1 + 3/var(x))
    return b

