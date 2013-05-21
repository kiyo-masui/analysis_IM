from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import math
import matplotlib.pyplot as plt
import numpy as np
import h5py
from pylab import *

########  Part for the cumulative Halo mass function  ########

def tinker_f(arr):
    """(array) -> function
    This function takes an array and returns interpolated function.
    """
    # linear version
    #f = interp1d(arr[:,0], arr[:,1])
    # log-log interpolation
    f = interp1d(np.log10(arr[:,0]), np.log10(arr[:,1]))
    return f

def integr_f(x, arr):
    """(number, array) -> float
    This function prepares tink_func for the integration. Just transfer to the log scale.
    """
    g = tinker_f(arr)
    return g(x)*math.log(10.)*10.**x


def cumul_halo_mf(x, arr):
    """(array, number) -> float
    This function takes data for interpolation and integrates log of the function.
    """
    return quad(integr_f, x, 16., epsabs=1.49e-08, args=(arr))

directory = "/cita/h/home-2/mufma/code/Tinker/test.dndM"
file = open(directory, "r")
halo_mf_data = np.genfromtxt(file)
#print cumul_halo_mf(15., halo_mf_data)

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
    """(number, number, number, number) -> float
    This function integrates the log of the schechter function and returns
    the value of the integral, with lower bound x.
    """
    return quad(schechter_f, x, 16.,
                args=(schechter_par1, schechter_par2, schechter_par3))[0]


'''
# debugging part
# integration methods don't reach the convergence
#print cumul_HI_mf(schechter_function,8.,0.006,-1.37,10**9.8)
#print cumul_halo_mf(halo_mf_data,14.8)
'''
########  Part for the abundance matching. Finding Halo mass from the given HI mass. ########

def equality(x, arr, HI_mass,
             param_1, param_2, param_3):
    """(number, array, number, number, number, number) -> number
    This function would be used to find zeros.
    """
    return cumul_HI_mf(HI_mass,param_1, param_2, param_3) -\
                cumul_halo_mf(x, arr)

# newton method is failing to converge after 50 steps.
def abundance_match(arr, axis_1):
    """(array, array, number, number, number) -> array
    This function for each given HI_mass in axis_1 finds corresponding Halo_mass
    and put it in axis_2, which is returned.
    """
    axis_2 = np.zeros_like(axis_1)
    i = 0
    for HI_mass in axis_1:
        axis_2[i] = brentq(equality, 10., 16., args=(arr,
                           HI_mass, 0.006,-1.37, 10**9.8), xtol=10**(-4))
        i = i + 1
    return axis_2

########  Part to work with catalogs  ########

def M_halo_to_M_HI(filename):
    """(str) -> Nonetype
    This function takes catalog and add a dataset of HI_masses.
    """
    catalog = h5py.File(filename)
    catalog.create_group('HI_Masses')
    Halo_mass = catalog['Halo_Masses']['Halo_Masses'][:]
    HI_mass = np.zeros(len(Halo_mass))
    for index in range(len(Halo_mass)):
        HI_mass[index] = 10**(func(math.log10(Halo_mass[index])))
    catalog['HI_Masses'].create_dataset('HI_Masses', data=HI_mass)

def M_HI_to_21(filename):
    """(str) -> Nonetype
    This function takes catalog and add dataset of 21_brightness.
    """
    catalog = h5py.File(filename)
    HI_mass = catalog['HI_Masses']['HI_Masses'][:]
    omega_m, omega_lambda, omega_HI = \
    0.27, 0.7, 4e-4
    HI_aver = np.sum(HI_mass)/len(HI_mass)
    z = 0.696
    coef = 0.29*1e3*((1.+z)*0.37/1.8)**0.5*omega_HI*(omega_m +
(1.+z)**(-3)*omega_lambda)**(-0.5)
    HI_bright = np.zeros(len(HI_mass))
    for index in range(len(HI_mass)):
        #in microKelvins
        HI_bright[index] = coef*(HI_mass[index] -
HI_aver)/HI_aver
    catalog.create_group('21_brightness_probe')
    catalog['21_brightness_probe'].create_dataset('21_probe',
data=HI_bright)


if __name__ == '__main__':
    #M_HI_to_21('/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5')
    #M_halo_to_M_HI('/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5')

    # interpolation of Tinker mass function
    # differential distribution for the halos
    # we used Tinker halo mass function(halo_mf) for this part
    directory = "/cita/h/home-2/mufma/code/Tinker/test.dndM"
    file = open(directory, "r")
    halo_mf_data = np.genfromtxt(file)
    z = np.zeros(50)
    for i in range(len(z)):
        z[i] = 8. + i*0.06
    x = abundance_match(halo_mf_data, z)
    func = interp1d(x, z)
    x = arange(11.,15.,0.01)
    y = func(x)
    plot(x, y)
    show()
