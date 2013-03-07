"""This script makes plots for putting limits on the HI powerspectrum from
both the cross power and the auto power."""

import numpy as np
from scipy import stats, interpolate, integrate, optimize
import matplotlib.pyplot as plt

#### Parameters ####

# Number of degrees of freedom for pair variance error bars.
NDEF = 5
# Measurements from the cross power paper.
CMEAN = 0.43
CERR = 0.07
LOW = 1e-10  # Small number in units K**2

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

case_15hr = 'pk_10modes.dat'
case_1hr = 'pk_45modes.dat'


#### Functions ####

def get_likelihood_omega(fname1, fname2, fname_sim):
    """Given two auto power measurements and a simulation, get likelihood
    function for b * Omega.
    
    In this funtion by omega, I mean Omega * b.
    """

    # Get x-axis.
    omega = np.arange(0.002, 2, 0.01)  # In units of 1e-3.
    omega_sq = omega**2
    n_omega = len(omega)
    # Load the data.
    k, auto_pow1, auto_std1, auto_gauss1 = get_auto(fname1)
    k2, auto_pow2, auto_std2, auto_gauss2 = get_auto(fname2)
    if not np.allclose(k, k2):
        raise RuntimeError('1 hour and 15 hour k bins not the same.')
    n_k = len(k)
    # Load the simulation (for omega=1e-3, which is 1 in the units we are
    # working in).
    k_sim, sim_pow = get_sim(fname_sim)
    sim_interp = interpolate.interp1d(k_sim, sim_pow)
    sim_pow = sim_interp(k)
    # Calculate the auto power probabilites. The is the probabilty of obtaining
    # the full 'k'-set of auto-power measurements at the given omegas (just the
    # product of the probabilities at each 'k'), with an integral over a
    # foreground power in each bin.
    p_auto_int_f_1 = 1.
    p_auto_int_f_2 = 1.
    for ii in range(n_k):
        p_auto_int_f_1 *= p_auto_given_signal_int_f(omega_sq, \
                auto_pow1[ii] / sim_pow[ii], auto_std1[ii] / sim_pow[ii])
        p_auto_int_f_2 *= p_auto_given_signal_int_f(omega_sq, \
                auto_pow2[ii] / sim_pow[ii], auto_std2[ii] / sim_pow[ii])
    # Calculate cross power pdf.
    p_cross_int_r = p_cross_given_auto_int_r(omega_sq, 1.)
    # Calculate prior.
    #prior = omega  # Flat prior on power.
    prior = np.ones(len(omega))  # Flat prior on Omega, 1/sqrt on power.
    # Construct total likelihood.
    likelihood = p_auto_int_f_1 * p_auto_int_f_2 * p_cross_int_r * prior
    likelihood_auto_1 = p_auto_int_f_1 * prior
    likelihood_auto_2 = p_auto_int_f_2 * prior
    likelihood_cross = p_cross_int_r * prior
    return omega, likelihood, likelihood_auto_1, \
            likelihood_auto_2, likelihood_cross

def get_conf_interval(lam, pdf, conf=0.68):
    """Get upper and lower confidance intervals given likelihood function."""
    # Get and normalize the cdf.
    cdf = np.zeros(len(lam))
    cdf[1:] = integrate.cumtrapz(pdf, lam)
    cdf /= cdf[-1]
    # Get the survival.
    surv = (1. - conf) / 2
    # Interpolate the power to the survival and 1 - survival.
    lam_int = interpolate.interp1d(cdf, lam)
    llim = lam_int(surv)
    ulim = lam_int(1. - surv)
    return llim, ulim

def signal_likelihood_k_bin(sim_pow, auto_pow1, auto_err1, auto_pow2,
                            auto_err2):
    """Get signal likelihood function in a k bin."""
    # Where to calculate the probability densities.
    signal_pow_min = sim_pow * CERR * 0.01
    signal_pow_max = min(auto_pow1 + 5 * auto_err1,
                         auto_pow2 + 5 * auto_err2)
    signal_pow = np.arange(signal_pow_min, signal_pow_max, signal_pow_min)
    # Get the pdf from the cross-power.
    p_cross_int_r = p_cross_given_auto_int_r(signal_pow, sim_pow)
    # Get the pdfs, from the auto-powers.
    p_auto_int_f_1 = p_auto_given_signal_int_f(signal_pow, auto_pow1, auto_err1)
    p_auto_int_f_2 = p_auto_given_signal_int_f(signal_pow, auto_pow2, auto_err2)
    # Calculate prior.
    #prior = np.ones(len(signal_pow))  # Flat prior on power, favors high Omega.
    prior = 1. /  np.sqrt(signal_pow)  # Flat prior on Omega.
    # Combine.
    likelihood = p_cross_int_r * p_auto_int_f_1 * p_auto_int_f_2 * prior
    likelihood = normalize_pdf(signal_pow, likelihood)
    return signal_pow, likelihood

def get_sim(fname):
    # Load the simulation data.
    data = np.loadtxt(fname)
    # Get the simulation curves.
    mask = np.isfinite(data[:,4])
    k = data[mask,1]
    P = data[mask,4]
    return k, P

def get_auto(fname):
    data = np.loadtxt(fname)
    # Get the data curves.
    mask = np.isfinite(data[:,1])
    k = data[mask,0]
    P = data[mask,1]
    std = data[mask,2] / np.sqrt(NDEF)
    gauss = data[mask,3]
    return k, P, std, gauss

def p_cross_given_auto_int_r(signal_pow, sim_pow):
    """Get probability of obtaining cross-power measurement, given auto-power.
    
    Gives the probability of obtaining the cross-power measurement in Masui,
    Switzer, et. al. 2012, as a function of underlying auto-power amplitude,
    then integrates out the correlation coefficient 'r'.
    Also requires a normalizing simulation with Omega * b = 1e-3.

    Function is unnormalized but has the right dependance on `signal_pow`.
    """
    # Scale simulation to measured lower limit from cross-power.
    cross_mean = CMEAN * np.sqrt(sim_pow)
    cross_err = CERR * np.sqrt(sim_pow)
    # See Feb 26, 2013 of Kiyo's notes for derivation.
    amp = np.sqrt(signal_pow)
    p_cross = (stats.norm.cdf((amp - cross_mean)/cross_err)
                - stats.norm.cdf(-cross_mean/cross_err))
    p_cross *= cross_err / amp
    return p_cross

def p_auto_given_signal_int_f(signal_pow, auto_mean, auto_err):
    """Get probability of obtaining auto-power measurement, given auto-power.
    
    Gives the probability of obtaining an auto-power measurement
    as a function of underlying signal-power amplitude, then integrates out
    the unknown foreground power.

    Function is unnormalized but has the right dependance on `signal_pow`.
    """
    # See Feb 26, 2013 of Kiyo's notes for derivation.
    p_auto = stats.t.sf((signal_pow - auto_mean) / auto_err, NDEF)
    return p_auto

def normalize_pdf(lam, pdf):
    """Normalize a function to integrate to unity."""
    norm = integrate.simps(pdf, lam)
    pdf = pdf / norm
    return pdf


#### Making the main plots ####

## Main plot - need to turn into a function. ##
#intervals = [.99, .95, .68]
#interval_cols = ['1.0', '0.8', '0.6', '0.4']
intervals = [.95, .68]
interval_cols = ['1.0', '0.7', '0.5']
imK2 = 1e6
# Set up plot.
f = plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel(r'$\Delta^2 ({\rm mK}^2)$')
plt.xlabel(r'$k ({\rm h/MPc})$')
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
# Load all the data.
k_sim, sim_pow = get_sim(data_file_root + sim_group + sim_case)
sim_interp = interpolate.interp1d(k_sim, sim_pow)
k, auto_pow15, auto_std15, auto_gauss15 = get_auto(data_file_root
                                                   + group_15hr + case_15hr)
k1, auto_pow1, auto_std1, auto_gauss1 = get_auto(data_file_root
                                                 + group_1hr + case_1hr)
if not np.allclose(k, k1):
    raise RuntimeError()
# Allocate memory of the confidance intervals.
n_intervals = len(intervals)
n_k = len(k)
ulim = np.zeros((n_intervals, n_k))
llim = np.zeros((n_intervals, n_k))
# Loop over k and calculate the confidance for each bin.
for ii in range(n_k):
#for ii in [6]:
    # Data points for this bin.
    this_k = k[ii]
    this_sim_pow = sim_interp(k[ii])
    this_auto_pow15 = auto_pow15[ii]
    this_auto_std15 = auto_std15[ii]
    this_auto_pow1 = auto_pow1[ii]
    this_auto_std1 = auto_std1[ii]
    pow, likelihood = signal_likelihood_k_bin(this_sim_pow, this_auto_pow15,
                                              this_auto_std15, this_auto_pow1,
                                              this_auto_std1)
    for jj in range(n_intervals):
        llim[jj,ii], ulim[jj,ii] = get_conf_interval(pow, likelihood,
                                                     intervals[jj])
# Now plot the limits.
for jj in range(n_intervals):
    plt.fill_between(k, imK2*ulim[jj,:], imK2*LOW, color=interval_cols[jj + 1])
for jj in range(n_intervals - 1, -1, -1):
    plt.fill_between(k, imK2*llim[jj,:], imK2*LOW, color=interval_cols[jj])
# Plot the auto-power upper limits.
cross_line = plt.plot(k_sim, imK2*CMEAN**2*sim_pow, 'k', lw=1.)
# Tidy up plot
xticks = np.arange(1, 10, 1)
xticks = list(0.01 * xticks) + list(0.1 * xticks) + list(xticks)
plt.xticks(xticks, xticks)
plt.ylim([imK2*1e-9, imK2*2e-6])
#plt.xlim([1.01*min(k), 0.99*max(k)])
plt.xlim([1.01*min(k), 0.7])



# Full Omega*b likelihood plot.
likelihoods = get_likelihood_omega(data_file_root + group_15hr + case_15hr, 
                                   data_file_root + group_1hr + case_1hr, 
                                   data_file_root + sim_group + sim_case)
omega, likelihood, likelihood_auto_1, likelihood_auto_2, likelihood_cross \
        = likelihoods
# Put omega in full units.
omega *= 1e-3
# For logarithmic axes.
#m = 1
#lomega = omega
#xscale = 'linear'
#xlim = [0.0, 1.5]
m = omega
lomega = np.log(omega)
xscale = 'log'
xlim = [0.00005, 0.002]
# Normalize all.
likelihood = normalize_pdf(lomega, m * likelihood)
likelihood_auto_1 = normalize_pdf(lomega, m * likelihood_auto_1)
likelihood_auto_2 = normalize_pdf(lomega, m * likelihood_auto_2)
likelihood_cross = normalize_pdf(lomega, m * likelihood_cross)
# Normalize the others to have same max.
#m_like = np.amax(likelihood)
#likelihood_auto_1 *= m_like / np.amax(likelihood_auto_1)
#likelihood_auto_2 *= m_like / np.amax(likelihood_auto_2)
#likelihood_cross *= m_like / np.amax(likelihood_cross)
# Plot.
plt.figure()
ax = plt.gca()
ax.set_xscale(xscale)
plt.ylabel(r'$\mathcal{L}\left[\rm{ln}(\Omega_{\rm HI} b_{\rm HI})\right]$')
plt.xlabel(r'$\Omega_{\rm HI} b_{\rm HI}$')
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
line_cross = plt.plot(omega, likelihood_cross, 'b--', lw=2)
line_auto1 = plt.plot(omega, likelihood_auto_1, 'g-.', lw=2)
line_auto2 = plt.plot(omega, likelihood_auto_2, 'r:', lw=2)
line = plt.plot(omega, likelihood, 'k', lw=2)
plt.legend(('cross-power', '15 hr auto-power', '1 hr auto-power', 'combined'),
           'upper left', prop={'size' : 14, 'family' : 'serif'})
plt.xlim(xlim)

plt.show()
