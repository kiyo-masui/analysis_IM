from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import interpld
from math import exp

# differential distribution for the hydrogen
# approximating it by the Schechter function
def schechter_function(observable,normalization,alpha,observable_star):
    return normalization*(observable/observable_star)**(alpha)*
           exp(-observable/observable_star)*(1/observable_star)

# differential distribution for the halos
# work required: determine which function which we use
def halos_function():
    return

# cumulative mass function for the neutral hydrogen
def cumulative_1(mass,differential_distribution,schechter_par1,
                 schechter_par2,schechter_par3):
    return quad(differential_distribution,mass,Inf,
                args=(schechter_par1,schechter_par2,schechter_par3))[0]

# cumulative mass function for the halos
def cumulative_2(mass,differential_distribution):
    return quad(differential_distribution,mass,Inf)[0]

# auxillary part to utilize the brentq
def equality(variable_mass,mass_1,param_1,param_2,param_3):
    return cumulative_1(mass_1,schechter_function,param_1,param_2,param_3) -
cumulative_2(variable_mass,halos_function)

# these part equate cumulative mass funcitons
# required work: determine what to put instead of 1000 for a good working of
# brentq
def abundance_match(axis_1,schechter_param1,schechter_param2,
                    schechter_param3):
    axis_2 = np.zeros_like(axis_1)
    i = 0
    for mass_1 in axis_1:
        axis_2[i] = brentq(equality,0,1000,args=(mass_1,schechter_param1,schechter_param2,schechter_param3))
        i = i + 1
    return axis_2
        # find mass_2 such that cumulative_2(mass_2) = cumulative_1(mass_1)
# these part to determine the dependency of hydrogen mass on halos mass
def interpolation(array1,array2):
    function = interpld(array1,array2)
    #return interpolated function

# Prototype:
# Mass_1_array = [1,2,3,4,5,6,7,8,9,10]
# print abundance_match(Mass_1_array,param1,param2,param3)   -  it returns a
# new array with matched masses for the hydrogen
# then we can interpolate these two arrays to some function
