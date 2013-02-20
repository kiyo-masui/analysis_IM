"""This script makes plots for putting limits on the HI powerspectrum from
both the cross power and the auto power."""

import numpy as np
import matplotlib.pyplot as plt

plt.ion()
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

data_file_root = '/cita/h/home-2/eswitzer/code/analysis_IM/pwrspec_plots/'

group_15hr = ('GBT_15hr_map_oldcalpolstack_x_GBT_15hr_map_oldcalpolstack'
                '_onesided_ixi_conv1p4_clipnoise_bfwindow_analysis/')

group_1hr = ('GBT_1hr_map_oldcalpolstack_x_GBT_1hr_map_oldcalpolstack'
               '_onesided_ixi_conv1p4_clipnoise_analysis/')

sim_group = ("GBT_1hr_map_oldcalpolstack_x_GBT_1hr_map_oldcalpolstack"
             "_onesided_iqu_conv1p4_clipnoise_puresignal/sim_nobeam/")
sim_case = "power_1d_from_2d_0modes.dat"

modes_removed = range(5, 90, 5)
cases = ['pk_%dmodes.dat' % n for n in modes_removed]

sim_k, sim_lim = get_lower_limit(data_file_root + sim_group + sim_case)

plt.figure()
plt.loglog(sim_k, sim_lim, '.')
for case in cases:
    k, l = get_upper_limit(data_file_root + group_15hr + case)
    plt.loglog(k, l)

plt.figure()
plt.loglog(sim_k, sim_lim, '.')
for case in cases:
    k, l = get_upper_limit(data_file_root + group_1hr + case)
    plt.loglog(k, l)

plt.figure()
for ii, this_k in enumerate(k[[1, 3, 5, 7, 9]]):
    this_P, this_std_err, this_gauss, this_ulim_95 = \
            get_k_data_for_cases(data_file_root + group_1hr, cases, this_k)
    plt.semilogy(modes_removed, this_ulim_95, colors[ii % len(colors)])
    plt.semilogy(modes_removed, this_ulim_95 - this_P,
                 '+' + colors[ii % len(colors)])

plt.figure()
for ii, this_k in enumerate(k[[1, 3, 5, 7, 9]]):
    this_P, this_std_err, this_gauss, this_ulim_95 = \
            get_k_data_for_cases(data_file_root + group_15hr, cases, this_k)
    plt.semilogy(modes_removed, this_ulim_95, colors[ii % len(colors)])
    plt.semilogy(modes_removed, this_ulim_95 - this_P,
                 '+' + colors[ii % len(colors)])



def get_k_data_for_cases(file_root, cases, k):

    P = []
    std_err = []
    gauss = []
    ulim_95 = []

    for case in cases:
        data = np.loadtxt(file_root + case)
        mask = np.isfinite(data[:,1])

        this_k = data[mask,0]

        this_P = data[mask,1]
        this_std = data[mask,2]
        this_gauss = data[mask,3]

        ind = np.where(this_k == k)[0][0]
        this_P = this_P[ind]
        this_std = this_std[ind]
        this_gauss = this_gauss[ind]

        # n - 1 is the correct factor to use here.
        this_std_err = this_std / np.sqrt(5)
        # Factor from the student's-t distribution with 5 degrees of freedom.
        # For many degrees of freedom, this factor would be approximately 2
        # (Gaussian erf for 95% confidance).  The excess comes from marginallizing
        # over the uncertainty in the variance.
        this_err_95 = 2.571 * this_std_err
        this_ulim_95 = this_P + this_err_95
        
        P.append(this_P)
        std_err.append(this_std_err)
        gauss.append(this_gauss)
        ulim_95.append(this_ulim_95)

    return np.array(P), np.array(std_err), np.array(gauss), np.array(ulim_95)




def get_upper_limit(fname):
    data = np.loadtxt(fname)

    # Get the data curves.
    mask = np.isfinite(data[:,1])

    k = data[mask,0]
    P = data[mask,1]
    std = data[mask,2]
    gauss = data[mask,3]

    # n - 1 is the correct factor to use here.
    std_err = std / np.sqrt(5)
    # Factor from the student's-t distribution with 5 degrees of freedom.
    # For many degrees of freedom, this factor would be approximately 2
    # (Gaussian erf for 95% confidance).  The excess comes from marginallizing
    # over the uncertainty in the variance.
    err_95 = 2.571 * std_err
    #err_95 = 2 * gauss

    # Looks like this is already included.
    #k_fact = k**3 / 2 / np.pi**2
    k_fact = 1.

    return k, k_fact * (P + err_95)

def plot_auto_power(fname):
    data = np.loadtxt(fname)

    # Get the data curves.
    mask = np.isfinite(data[:,1])

    k = data[mask,0]
    P = data[mask,1]
    std = data[mask,2]
    gauss = data[mask,3]

    # n - 1 is the correct factor to use here.
    std_err = std / np.sqrt(5)
    # Factor from the student's-t distribution with 5 degrees of freedom.
    # For many degrees of freedom, this factor would be approximately 2
    # (Gaussian erf for 95% confidance).  The excess comes from marginallizing
    # over the uncertainty in the variance.
    err_95 = 2.571 * std_err

    # Looks like this is already included.
    #k_fact = k**3 / 2 / np.pi**2
    k_fact = 1.

    plt.errorbar(k, k_fact * P, k_fact * err_95)


def get_lower_limit(fname):
    sim_data = np.loadtxt(fname)

    # Get the simulation curves.
    sim_mask = np.isfinite(sim_data[:,4])

    sim_k = sim_data[sim_mask,1]
    sim_P = sim_data[sim_mask,4]

    #sim_k_fact = sim_k**3 / 2 / np.pi**2
    sim_k_fact = 1.
    # Measurement from cross power minus 2 sigma.
    sim_fact = 0.29**2

    return sim_k, sim_P * sim_fact * k_fact


#### OLD ####

group = ('GBT_15hr_map_oldcalpolstack_x_GBT_15hr_map_oldcalpolstack'
         '_onesided_ixi_conv1p4_clipnoise_bfwindow_analysis/')
case = 'pk_20modes.dat'

sim_group = ("GBT_1hr_map_oldcalpolstack_x_GBT_1hr_map_oldcalpolstack"
             "_onesided_iqu_conv1p4_clipnoise_puresignal/sim_nobeam/")
sim_case = "power_1d_from_2d_0modes.dat"


# Read data files.
data = np.loadtxt(data_file_root + group + case)
sim_data = np.loadtxt(data_file_root + sim_group + sim_case)

# Get the simulation curves.
sim_mask = np.isfinite(sim_data[:,4])

sim_k = sim_data[sim_mask,1]
sim_P = sim_data[sim_mask,4]

#sim_k_fact = sim_k**3 / 2 / np.pi**2
sim_k_fact = 1.
sim_fact = 0.46**2

# Get the data curves.
mask = np.isfinite(data[:,1])

k = data[mask,0]
P = data[mask,1]
std = data[mask,2]
gauss = data[mask,3]

# n - 1 is the correct factor to use here.
std_err = std / np.sqrt(5)
# Factor from the student's-t distribution with 5 degrees of freedom.
err_95 = 2.571 * std_err

# Looks like this is already included.
#k_fact = k**3 / 2 / np.pi**2
k_fact = 1.

# Quick look at the PS
plt.figure()
plt.loglog(sim_k, sim_fact * sim_k_fact * sim_P)
plt.loglog(k, (P) * k_fact, '.')
plt.loglog(k, (P + err_95) * k_fact)
plt.loglog(k, std_err * k_fact)

# Compare the 2 types of errorbar.
plt.figure()
plt.loglog(k, 2 * gauss)
plt.loglog(k, err_95)

