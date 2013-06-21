from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import math
import matplotlib.pyplot as plt
import numpy as np
import h5py
from pylab import *

########  Part for the cumulative halo mass function  ########

def tinker_f(arr):
    """ (array) -> function
    This function takes an array and returns interpolated function.
    """
    # linear version
    #f = interp1d(arr[:,0], arr[:,1])
    # log-log interpolation
    f = interp1d(np.log10(arr[:,0]), np.log10(arr[:,1]))
    return f

def integr_f(x, arr):
    """ (number, array) -> float
    This function prepares tinker_f for the integration. Just transfer to the
    log scale.
    """
    g = tinker_f(arr)
    return 10**g(x)*math.log(10.)*10.**x


def cumul_halo_mf(x, arr):
    """ (array, number) -> float
    This function takes data for interpolation and integrates log of the function.
    """
    return quad(integr_f, x, 16., epsabs=1.49e-08, args=(arr))[0]

directory = "/cita/h/home-2/mufma/code/Tinker/test.dndM"
file = open(directory, "r")
halo_mf_data = np.genfromtxt(file)

########  Part for the cumulative HI mass function  ########

def schechter_f(x, norm, alpha, mass_star):
    """ (number, number, number, number) -> float
    This function calculates the differential distribution for the HI.
    We approximate it by the Schechter function.
    In our notation 10**x is a mass of the cluster. This notation
    helps to do integration, since function decrease "rapidly".
    """
    ratio = (10.**x)/mass_star
    return (norm/mass_star)*math.log(10)*\
                (10.**x)*(ratio)**(alpha)*math.exp(-ratio)

def cumul_HI_mf(x, schechter_par1,
                schechter_par2, schechter_par3):
    """ (number, number, number, number) -> float
    This function integrates the log of the schechter function and returns
    the value of the integral, with lower bound x.
    """
    return quad(schechter_f, x, 16.,
                args=(schechter_par1, schechter_par2, schechter_par3))[0]


"""
# debugging part
# integration methods don't reach the convergence
print cumul_HI_mf(10.4345,0.006,-1.37,10**9.8)
print cumul_halo_mf(14.0,halo_mf_data)
"""

########  Part for the abundance matching. Finding halo mass from the given HI mass. ########

def equality(x, arr, HI_mass,
             param_1, param_2, param_3):
    """ (number, array, number, number, number, number) -> number
    This function would be used to find zeros.
    """
    return cumul_HI_mf(HI_mass,param_1, param_2, param_3) -\
           cumul_halo_mf(x, arr)

# newton method is failing to converge after 50 steps.
def abundance_match(arr, axis_1):
    """ (array, array) -> array
    This function for each given HI_mass in axis_1 finds corresponding Halo_mass
    and put it in axis_2, which is returned.
    """
    axis_2 = np.zeros_like(axis_1)
    i = 0
    for HI_mass in axis_1:
        axis_2[i] = brentq(equality, 10., 16., args=(arr,
                           HI_mass, 0.006,-1.37, 10**9.8), xtol=10**(-8))
        i = i + 1
    return axis_2

########  Part to work with catalogs.  ########

def interpol_f():
    """ (None) -> function
    Interpolation of Tinker mass function
    which is a differential distribution for the halos
    we used Tinker halo mass function(halo_mf) for this part
    """
    directory = "/cita/h/home-2/mufma/code/Tinker/test.dndM"
    file = open(directory, "r")
    halo_mf_data = np.genfromtxt(file)
    z = np.arange(8.,12.,0.03)
    x = abundance_match(halo_mf_data, z)
    f = interp1d(x, z)
    file.close()
    return f

def M_halo_to_M_HI(filename, func, regim):
    """ (str, function) -> Nonetype
    This function takes catalog and add a dataset of HI_masses.
    """
    catalog = h5py.File(filename,'%s'%regim)
    catalog.create_group('HI_Masses')
    Halo_mass = catalog['Halo_Masses']['Halo_Masses'][:]
    HI_mass = 10**(func(np.log10(Halo_mass)))
    catalog['HI_Masses'].create_dataset('HI_Masses', data=HI_mass)
    catalog.close()

def M_HI_to_21(filename, redshift, regim):
    """ (str) -> Nonetype
    This function takes catalog and add dataset of 21_brightness.
    """
    catalog = h5py.File(filename,'%s'%regim)
    HI_mass = catalog['HI_Masses']['HI_Masses'][:]
    omega_m, omega_lambda, omega_HI = \
    0.27, 0.7, 4e-4
    HI_aver = np.mean(HI_mass)
    z = redshift
    coef = 0.29*1e3*((1.+z)*0.37/1.8)**0.5*omega_HI*(omega_m +
(1.+z)**(-3)*omega_lambda)**(-0.5)
    #in microKelvins
    HI_bright = coef*(HI_mass - HI_aver)/HI_aver
    catalog.create_group('T_b21')
    catalog['T_b21'].create_dataset('T_b21', data=HI_bright)
    catalog.close()


if __name__ == '__main__':
    # here is the pipeline for the joe's catalogs
    """
    path = "/mnt/raid-project/gmrt/mufma/h5py_catalogs/"
    redshift = [0.042,0.086,0.130,0.220,0.267,0.316,0.416,0.468,0.523,0.636,
                0.696,0.758,0.889,0.958,1.030,1.185,1.267,1.353,1.538,1.638,
                1.742,1.969,2.092,2.222,2.505,2.660,2.825]
    for list in range(1,66),range(67,77),range(78,101):
        for num in list:
            print "finished catalog number %d"%num
            for red in redshift:
                M_halo_to_M_HI(path + 'simulation%d/'%num +
                               '%.3fhalo_catalog.hdf5'%red, f)
                M_HI_to_21(path + 'simulation%d/'%num +
                           '%.3fhalo_catalog.hdf5'%red, red)
    """
