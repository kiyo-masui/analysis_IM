from scipy.integrate import *
import numpy as np
from abundance_matching import tinker_f
import math
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib import rc

# Data from the Tinker code
file = open("/cita/h/home-2/mufma/code/pilot/data/halo_mf.dndM", "r")
tdata = np.genfromtxt(file)
log_halo_mf = interp1d(np.log10(tdata[:,0]),np.log10(tdata[:,1]))
halo_mf = interp1d(tdata[:,0], tdata[:,1])

# Function to do the log integral
def log_log_f(x):
    return 10.**log_halo_mf(x)*10.**x*math.log(10)

def log_f(x):
    return halo_mf(10.**x)*10.**x*math.log(10)

# Printing the interpolated function
if sys.argv[1] == "inter":
    x = np.arange(11.,16.,0.01)
    y = log_log_f(x)
    plt.plot(x,y)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()


# Comparing different integration methods
if sys.argv[1] == "print":
    power = float(sys.argv[2])
    bound = 10**power
    print "quad : ", quad(halo_mf, bound, 1e16)
    print "fixed_quad : ", fixed_quad(halo_mf, bound, 1e16)
    print "quadrature : ", quadrature(halo_mf, bound, 1e16)
    print "romberg : ", romberg(halo_mf, bound, 1e16)
    print "log_f quad : ", quad(log_f, power, 16.)
    print "log_f fixed_quad : ", fixed_quad(log_f, power, 16.)
    print "log_f quadrature : ", quadrature(log_f, power, 16.)
    print "log_f romberg : ", romberg(log_f, power, 16.)

# Thorough test for the log_f methods
if sys.argv[1] == "plot_err":
    qual = 200
    start = 14.0
    quadrature_info = np.zeros(qual)
    quad_info = np.zeros(qual)
    quad_log_info = np.zeros(qual)
    quadrature_log_info = np.zeros(qual)
    step = 1./qual
    points = np.arange(0,1,step)

    for k in range(qual):
        num_1, num_2 = quadrature(log_log_f, start + k*step, 16.)
        quadrature_log_info[k] = num_1/num_2
        num_1, num_2 = quad(log_log_f, start + k*step, 16.)
        quad_log_info[k] = num_1/num_2
        num_1, num_2 = quad(log_f, start + k*step, 16.)
        quad_info[k] = num_1/num_2
        num_1, num_2 = quadrature(log_f, start + k*step, 16.)
        quadrature_info[k] = num_1/num_2

    plt.plot(points, quadrature_info, 'b-o')
    plt.plot(points, quadrature_log_info, 'g-o')
    plt.plot(points, quad_info, 'r-o')
    plt.plot(points, quad_log_info, 'y-o')
    plt.show()

if sys.argv[1] == "tinker_expect":
    qual = 200
    start = 14.0
    vol = 126.25**3
    quadrature_info = np.zeros(qual)
    quad_info = np.zeros(qual)
    quad_log_info = np.zeros(qual)
    quadrature_log_info = np.zeros(qual)
    step = 1./qual
    points = np.arange(0,1,step)

    for k in range(qual):
        num_1 = quadrature(log_log_f, start + k*step, 16.)[0]
        quadrature_log_info[k] = num_1
        num_1 = quad(log_log_f, start + k*step, 16.)[0]
        quad_log_info[k] = num_1
        num_1 = quad(log_f, start + k*step, 16.)[0]
        quad_info[k] = num_1
        num_1 = quadrature(log_f, start + k*step, 16.)[0]
        quadrature_info[k] = num_1
    rc('text', usetex=True)
    plt.figure()
    plt.subplot(1,2,1)
    plt.title(r'$\textrm{Cumulative Tinker}$')
    plt.plot(points, vol*quadrature_info, 'b-o')
    plt.plot(points, vol*quadrature_log_info, 'g-o')
    plt.plot(points, vol*quad_info, 'r-o')
    plt.plot(points, vol*quad_log_info, 'y-o')
    tick_locs = [0,0.2,0.4,0.6,0.8,1.0]
    tick_lbls = [14.0, 14.2, 14.4, 14.6, 14.8, 15.0]
    plt.xticks(tick_locs, tick_lbls)
    plt.xlabel(r'$\mathrm{log} \left[\frac{M}{h^{-1}M_{\bigodot}} \right ]$')
    plt.ylabel(r'$\left \langle N \right \rangle \textrm{number of clusters}$')
    plt.subplot(1,2,2)
    plt.title(r'$\textrm{Volume} \left[\mathrm{126h^{-1}Mpc}\right ]^{3}$')
    plt.plot(points, vol*quadrature_info, 'b-o')
    plt.plot(points, vol*quadrature_log_info, 'g-o')
    plt.plot(points, vol*quad_info, 'r-o')
    plt.plot(points, vol*quad_log_info, 'y-o')
    tick_locs = [0,0.2,0.4,0.6,0.8,1.0]
    tick_lbls = [14.0, 14.2, 14.4, 14.6, 14.8, 15.0]
    plt.xticks(tick_locs, tick_lbls)
    plt.xlabel(r'$\mathrm{log} \left[\frac{M}{h^{-1}M_{\bigodot}} \right]$')
    #plt.ylabel("volume * integral(dn/dM)")
    plt.show()
