import numpy as np
import camb_caller as cc
import delta_bias as db
from matplotlib import rc
import matplotlib.pyplot as plt
import time

start = time.time()
# Set of parameters

k_min   = -3.
k_max   = -1.
number  = 50
output  = np.zeros(number)
k_range = np.linspace(k_min, k_max, number)
delta_k = np.ones_like(k_range)*(10**(k_max)-10**(k_min))/number
pwrspec = np.zeros(number)

for num in range(number):
    cosmo_params = [k_range[num], 0.000001]
    output[num]  = db.pwrspec_grad(cosmo_params)[1]
    pwrspec[num] = db.pwrspec_f_nl(cosmo_params)
P, T    = cc.camb_pipe({'hubble': 70})
P_k     = 10**P(k_range)
"""
cosmo_params = [-3, 0.001]
val1 = db.pwrspec_deriv(cosmo_params)[1]
val2 = db.pwrspec_grad(cosmo_params)[1]
"""

part1 = ((10**k_range)**2/(2*9.8))*delta_k*0.34*4.*10**10
part2 = (P_k+2.5*10**3)**(-2)
part3 = output**2
fisher  = np.sum(part1*part2*part3)
print "fisher matrix is :", fisher
print "sigma", 1/np.sqrt(fisher)
print "\n\ntime elapsed:", time.time() - start, "sec."
plt.plot(10**k_range, pwrspec)
plt.xscale('log')
plt.yscale('log')
plt.show()
