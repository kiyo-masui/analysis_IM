""""This script post processes auto-correlations into final results."""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


plt.close('all')
plt.ion()

file_root = os.getenv("GBT_YL")

power_home = file_root + "ps_result/"
transfer_home = file_root + "ps_result/bias/"

#case = "auto_GBT_15hr_map_oldcal_legendre_modes_0gwj_50"
#case = "auto_GBT_15hr_map_oldcal_legendre_modes_5gwj_25"
case = "auto_1hr_IE_legendre_modes_0gwj_conv_20"

power_root = power_home + 'power_' + case + '/' + case
transfer_root = transfer_home + case + '/'
#simulation_root = power_home + 'simulation_auto_sim_15hr_oldmap_str_25/'
simulation_root = power_home + ('reference_auto_1hr_IE_legendre_'
                                'modes_0gwj_conv_15/')

power_file = power_root + '_p2_combined.npy'
power_error_file = power_root + '_p2_var_combined.npy'
power_bins_file = power_root + '_k2_combined.npy'
transfer_file = transfer_root + 'b2_bias.npy'
simulation_file = simulation_root + 'simmaps_p2_combined.npy'
sim_bins_file = simulation_root + 'simmaps_k2_combined.npy'

power = np.load(power_file)
bins = np.load(power_bins_file)
transfer = np.load(transfer_file)
error = np.load(power_error_file)
simulation = np.load(simulation_file)
sim_bins = np.load(sim_bins_file)

power *= transfer
error *= transfer

k_t, k_p = np.meshgrid(bins[0], bins[1])

bad = simulation==0

#Z = simulation.copy()
#Z[bad] = 1e-10
#Z = np.log10(Z)

A = power/simulation
A[bad] = 0
Amin = -4
Amax = 4

plt.figure()
plt.pcolor(k_t, k_p, A, vmin=Amin, vmax=Amax)
plt.xscale('log')
plt.yscale('log')
plt.colorbar()


B = 1/transfer
Bmin=0
Bmax=1

plt.figure()
plt.pcolor(k_t, k_p, B, vmin=Bmin, vmax=Bmax)
plt.xscale('log')
plt.yscale('log')
plt.colorbar()


C = np.log10(error/simulation)
Cmin=-1
Cmax=2

plt.figure()
plt.pcolor(k_t, k_p, C, vmin=Cmin, vmax=Cmax)
plt.xscale('log')
plt.yscale('log')
plt.colorbar()


D = np.log10(error)
Dmin=-10
Dmax=-5

plt.figure()
plt.pcolor(k_t, k_p, D, vmin=Dmin, vmax=Dmax)
plt.xscale('log')
plt.yscale('log')
plt.colorbar()

