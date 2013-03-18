from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import math
import pylab as pl
import numpy as np
import h5py

# differential distribution for the neutral hydrogen
# we approximate it by the Schechter function


def schechter_function(observable, normalization, alpha, observable_star):
    ratio = (10.**observable)/observable_star
    return (normalization/observable_star)*math.log(10)*
                (10.**observable)*(ratio)**(alpha)*math.exp(-ratio)

# differential distribution for the halos
# we used Tinker halo mass function(halo_mf) for this part
directory = "/cita/h/home-2/mufma/code/Tinker/test.dndM"
file = open(directory, "r")
halo_mf_data = np.genfromtxt(file)

# interpolation of Tinker mass function


def tink_func(halo_array):
    x = np.zeros(len(halo_array))
    y = np.zeros(len(halo_array))
    for i in range(len(halo_array)):
        x[i] = halo_array[i][0]
        y[i] = halo_array[i][1]
    f = interp1d(x, y)
    return f

# defining function for the cumulative halo mass function


def integr_func(x, halo_array):
    g = tink_func(halo_array)
    return g((10.**x)*math.log(10.)*10.**x)

# cumulative mass function for the neutral hydrogen


def cumul_HI_mf(differential_distribution, mass, schechter_par1,
                schechter_par2, schechter_par3):
    return quad(differential_distribution, mass, 16.,
                args=(schechter_par1, schechter_par2, schechter_par3))[0]

# cumulative mass function for the halos
# we do interpolation of Tinker mass function and then integrate it.


def cumul_halo_mf(halo_array, mass):
    return quad(integr_func, mass, 16., epsabs=1.49e-08, args=(halo_array))[0]

"""
# debugging part
# integration methods don't reach the convergence
print cumul_HI_mf(schechter_function,8.,0.006,-1.37,10**9.8)
print cumul_halo_mf(halo_mf_data,15.5)
"""

# auxillary part for root finding
# we will look for such x that equality(x)=0


def equality(halo_array, variable_mass_halo, mass_HI,
             param_1, param_2, param_3):
    return cumul_HI_mf(schechter_function, mass_HI,
                       param_1, param_2, param_3) -
                cumul_halo_mf(halo_array, variable_mass_halo)


# these part equate cumulative mass funcitons
# required work: how to find root of the equation equality(x)=0
# newton method is failing to converge after 50 iterations

def abundance_match(halo_array, axis_1, schechter_param1, schechter_param2,
                    schechter_param3):
    axis_2 = np.zeros_like(axis_1)
    i = 0
    for mass_1 in axis_1:
        axis_2[i] = brentq(equality, 10., 16., args=(halo_array, mass_1, 0.006,
                           -1.37, 10**9.8), xtol=10**(-4))
        i = i + 1
    return axis_2
        # find mass_2 such that cumulative_2(mass_2) = cumulative_1(mass_1)
# these part to determine the dependency of hydrogen mass on halos mass
#def interpolation(array1,array2):
#    function = interp1d(array1,array2)
    #return interpolated function

z = np.zeros(50)
for i in range(len(z)):
    z[i] = 8. + i*0.06
x = abundance_match(halo_mf_data, z, 0.006, -1.37, 10**9.8)
func = interp1d(x, z)


def M_halo_to_M_HI(filename):
    create_file = h5py.File(filename)
    create_file.create_group('HI_Masses')
    Halo_mass = create_file['Halo_Masses']['Halo_Masses'][:]
    print Halo_mass.max()
    print Halo_mass.min()
    HI_mass = np.zeros(len(Halo_mass))
    for index in range(len(Halo_mass)):
        HI_mass[index] = 10**(func(math.log10(Halo_mass[index])))
    create_file['HI_Masses'].create_dataset('HI_Masses', data=HI_mass)
M_halo_to_M_HI('/cita/h/home-2/mufma/code/analysis_IM/Tryout.hdf5')
