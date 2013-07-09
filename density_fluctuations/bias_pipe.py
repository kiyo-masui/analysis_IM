import numpy as np
import camb_caller as cc
import delta_bias as db
from matplotlib import rc
import matplotlib.pyplot as plt
import time
"""
rc('text', usetex=True)
xlabel(r"$\mathrm{k}\left[\frac{Mpc}{h}\right]$",size=15)
"""
start = time.time()
# Set of parameters
k_min   = -3
k_max   = -1
number  = 20
output  = np.zeros(number)
k_range = np.linspace(k_min, k_max, number)
delta_k   = np.ones_like(k_range)*(10**(k_max)-10**(k_min))/number
for num in range(number):
    cosmo_params = [k_range[num], 10.]
    output[num]  = db.pwrspec_deriv(cosmo_params)[1]
P, T    = cc.camb_pipe({'hubble': 70})
P_k     = 10**P(k_range)
fisher  = np.sum(((10**k)**2)*(P_k+2.5*10**3)**(-2)*output*delta_k*output*8*10**9)/(2*9.8)
print "fisher matrix is :", fisher
print "\n\ntime elapsed:", start - time.time(), "sec."
