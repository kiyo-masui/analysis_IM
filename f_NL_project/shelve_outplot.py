import shelve
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from scipy.interpolate import interp1d
rc('text', usetex=True)

def open_shelve(filename):
    """
    Open shelve file and return list of data
    """
    file = shelve.open(filename)
    data = file['sim:phys'][1][1]
    points = data['bin_center'][:]
    values = data['binavg'][:]
    file.close()
    return points, values

# Start a figure. Set scales and labels.
fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)
axs.set_title(r"$\mathrm{Output~of~program~pwrspec~combinations.~z\approx0.8}$",size=20)
axs.set_ylabel(r"$\frac{\mathit{P}(\mathrm{k})\mathrm{k}^{3}}{2\pi^{2}}$",size=22)
axs.set_xlabel(r"$\mathrm{k}\left[\frac{Mpc}{h}\right]$",size=15)
axs.set_xscale('log')
axs.set_yscale('log')

def plot_r_val(filename_1, filename_2, filename_3):
    """
    Plot r values and output it for analysis
    """
    map1 = open_shelve(filename_1)
    map2 = open_shelve(filename_2)
    map3 = open_shelve(filename_3)
    r_val = map1[1]/np.sqrt(map2[1]*map3[1])
    axs.set_ylabel(r"$\frac{P_{gal\times 21}}{\sqrt{P_{gal\times gal}P_{21\times 21}}}$",size=22)
    p, = axs.plot(map1[0], r_val, 'r-o')
    return p

def plot_pwrspec(filename):
    k_val, pwrspec = open_shelve(filename)
    p, = axs.plot(k_val, pwrspec, 'o')
    return p

def plot_camb():
    root_c = '/mnt/raid-project/gmrt/mufma/CAMB/'
    camb = np.loadtxt(root_c + 'test_matterpower0.889.dat')
    arr_k = camb[:,0]
    value_k = camb[:,1]
    value_k = value_k*arr_k**3/(2*np.pi**2)
    p, = axs.plot(arr_k,value_k)
    return p

def plot_legend(list_p, list_info):
    axs.legend(list_p, list_info, loc=2)

def plot_bias(filename1):
    k_val1, pwrspec1 = open_shelve(filename1)
    root_c = '/mnt/raid-project/gmrt/mufma/CAMB/'
    camb = np.loadtxt(root_c + 'test_matterpower0.889.dat')
    arr_k = camb[:,0]
    value_k = camb[:,1]
    pwrspec_camb = value_k*arr_k**3/(2*np.pi**2)
    P = interp1d(arr_k,pwrspec_camb)
    pwrspec2 = P(k_val1)
    p, = axs.plot(k_val1, pwrspec1/pwrspec2)
    return p
