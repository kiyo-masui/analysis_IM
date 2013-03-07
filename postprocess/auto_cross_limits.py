"""This script makes plots for putting limits on the HI powerspectrum from
both the cross power and the auto power."""

import numpy as np
from scipy import stats, interpolate, integrate, optimize
import matplotlib.pyplot as plt

plt.ion()
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
line_styles = ['-', '--', '-.', ':']

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

plt.rcParams.update({'font.size': 16})


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
line_cross = plt.plot(omega, likelihood_cross, 'b--', lw=2)
line_auto1 = plt.plot(omega, likelihood_auto_1, 'g-.', lw=2)
line_auto2 = plt.plot(omega, likelihood_auto_2, 'r:', lw=2)
line = plt.plot(omega, likelihood, 'k', lw=2)
plt.legend(('cross-power', '15 hr auto-power', '1 hr auto-power', 'combined'),
           'upper left', prop={'size' : 14, 'family' : 'serif'})
plt.xlim(xlim)


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


#### OLD ####

## Main plot ##
intervals = [.95, .68]
#intervals = [.68]
interval_cols = ['1.0', '0.7', '0.5']
imK2 = 1e6
# Set up plot.
f = plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
plt.ylabel(r'$\Delta^2 ({\rm mK}^2)$')
plt.xlabel(r'$k ({\rm h/MPc})$')
# Plot the auto-power upper limits.
for ii, inter in enumerate(intervals):
    k, ulim = get_mutual_uppper_limit(data_file_root + group_15hr
             + case_15hr, data_file_root + group_1hr + case_1hr, inter, pw=1)
    plt.fill_between(k, imK2*ulim, imK2*LOW, color=interval_cols[ii + 1])
    # Put markers on upper limits to indicate which is coming from 1hr and
    # 15hr.
    #mask_1 = np.logical_not(mask15)
    #plt.plot(k[mask15], ulim[mask15], 'kx', markersize=7)
    #plt.plot(k[mask1], ulim[mask1], 'k+', markersize=7)
# Plot the cross-power lower limits.
intervals.reverse()
interval_cols.reverse()
for ii, inter in enumerate(intervals):
    k_cross, llim = get_cross_lower_limit(data_file_root + sim_group
                                                + sim_case, inter)
    plt.fill_between(k_cross, imK2*llim, imK2*LOW, color=interval_cols[ii + 1])
# Plot the r=1 cross power line.
k_cross, P_cross  = get_cross_lower_limit(data_file_root + sim_group
                                                + sim_case, 0.)
cross_line = plt.plot(k_cross, imK2*P_cross, 'k', lw=1.)
# Plot the best possible upper limit line.
k, err68 = get_mutual_uppper_limit(data_file_root + group_15hr
             + case_15hr, data_file_root + group_1hr + case_1hr, 0.68, pw=0)
P_cross_interp = interpolate.interp1d(k_cross, P_cross)
P_cross_at_k = P_cross_interp(k)
best_possible_lim = P_cross_at_k + err68
plt.plot(k, imK2*best_possible_lim, 'k--', lw=1.)
#mask_1 = np.logical_not(mask15)
#plt.plot(k[mask15], best_possible_lim[mask15], 'kx', markersize=7)
#plt.plot(k[mask1], best_possible_lim[mask1], 'k+', markersize=7)
# Tidy up plot
xticks = np.arange(1, 10, 1)
xticks = list(0.01 * xticks) + list(0.1 * xticks) + list(xticks)
plt.xticks(xticks, xticks)
plt.ylim([imK2*1e-9, imK2*2e-6])
#plt.xlim([1.01*min(k), 0.99*max(k)])
plt.xlim([1.01*min(k), 0.7])

def get_cross_lower_limit(fname, confidance=.95):
    
    # Get the simulation curves.
    k, P = get_sim(fname)
    
    # Simulation has Omega_HI=1e-3.  Scale to measured lower limit.
    cross_mean = 0.43
    cross_err = 0.07

    # Get the number of sigmas for the confidance interval.
    one_sided_conf = 1 - (1. - confidance) / 2
    z = stats.norm.ppf(one_sided_conf)

    # Measurement from cross power minus 2 sigma.
    sim_fact = (cross_mean - z * cross_err)**2

    return sim_k, P * sim_fact

def get_auto_upper_limit(fname, confidance=0.95, pw=1):
    # Get the data.
    k, P, std, gauss = get_auto(fname)

    # Figure out the student-t factor.
    one_sided_conf = 1 - (1. - confidance) / 2
    t = stats.t.ppf(one_sided_conf, NDEF)

    return k, pw * P + t * std
    
def get_mutual_uppper_limit(fname1, fname2, confidance=0.95, pw=1):
    # Old way
    ## We need slightly higher confidance to deal with selection bias of 2
    ## choices.
    #confidance = 1 - (1. - confidance) / 2
    #k1, l1, = get_auto_upper_limit(fname1, confidance, pw=pw)
    #k2, l2, = get_auto_upper_limit(fname2, confidance, pw=pw)
    #if not np.allclose(k1, k2):
    #    raise NotImplementedError('k values do not agree between 2 spectra.')
    ## Figure out which k bins the first is better.
    #mask1 = l1 < l2
    #l2[mask1] = l1[mask1]
    #return k1, l2, mask1
    # New way
    k1, P1, std1, gauss1 = get_auto(fname1)
    k2, P2, std2, gauss2 = get_auto(fname2)
    if not np.allclose(k1, k2):
        raise NotImplementedError('k values do not agree between 2 spectra.')
    lim = np.empty(len(k1))
    for ii in range(len(lim)):
        lim[ii] = joint_upper_limit(pw*P1[ii], std1[ii], pw*P2[ii], std2[ii],
                                    confidance)
    return k, lim


def joint_upper_limit_old(m1, s1, m2, s2, confidance=0.95):
    # Convert confidance to a desired one sided survival.
    survival = (1. - confidance) / 2
    ## First get the integration step sizes and limits.
    #step = min(s1, s2) / 10.
    #int_lim = max(m1 + 5 * s1, m2 + 5 * s2)
    #lam = np.arange(0, int_lim, step)
    # From To combine constraints, we multiply the survival function.
    # See Feb 216, 2013 of Kiyo's notes.
    sft = lambda lam : (stats.t.sf((lam - m1) / s1, NDEF)
                        * stats.t.sf((lam - m2) / s2, NDEF))
    # Applying prior that both signal and foregrounds are positive is
    # equivalent to normalizing sft(0) to be unity.
    sft_prior = lambda lam : sft(lam) / sft(0)
    # Equation to invert.
    eqn = lambda lam : sft_prior(lam) - survival
    # Use Newton Ralphson to invert the function.
    lam0 = optimize.brentq(eqn, min(m1 - 5*s1, m2 - 5*s2),
                           min(m1 + 5*s1, m2 + 5*s2), rtol=1e-4)

    #sf1 = stats.t.sf((lam - m1) / s1, NDEF)
    #sf2 = stats.t.sf((lam - m2) / s2, NDEF)
    ## Joint distribution.
    #sft = sf1 * sf2
    #print len(p_new), len(pdft)
    return lam0
    


sim_k, sim_lim = get_lower_limit(data_file_root + sim_group + sim_case)

plt.figure()
plt.loglog(sim_k, sim_lim, '.')
for case in cases[1:2]:
    k, l = get_upper_limit(data_file_root + group_15hr + case)
    plt.loglog(k, l)
plt.ylabel(r'$\Delta^2 ({\rm K}^2)$')
plt.xlabel(r'$k ({\rm h/MPc})')

#plt.figure()
plt.loglog(sim_k, sim_lim, '.')
for case in cases[8:9]:
    k, l = get_upper_limit(data_file_root + group_1hr + case)
    plt.loglog(k, l)
plt.ylabel(r'$\Delta^2 ({\rm K}^2)$')
plt.xlabel(r'$k ({\rm h/MPc})')

k_to_plot = k[[2, 4, 6]]
n_k_to_plot = len(k_to_plot)

plt.figure()
#plt.subplot(n_k_to_plot, 1, 0)
for ii, this_k in enumerate(k_to_plot):
    #plt.subplot(n_k_to_plot, 1, ii)
    this_P, this_std_err, this_gauss, this_ulim_95 = \
            get_k_data_for_cases(data_file_root + group_1hr, cases, this_k)
    thus_ulim_95 = this_P + 2 * this_gauss
    plt.semilogy(modes_removed, this_ulim_95, lw=3.,
                 color=colors[ii % len(colors)],
                 linestyle=line_styles[ii % len(line_styles)])
    plt.semilogy(modes_removed, this_ulim_95 - this_P, lw=0.7,
                 color=colors[ii % len(colors)],
                 linestyle=line_styles[ii % len(line_styles)])
plt.ylabel(r'$\Delta^2 ({\rm K}^2)$')
plt.xlabel('modes removed')

plt.figure()
for ii, this_k in enumerate(k_to_plot):
    this_P, this_std_err, this_gauss, this_ulim_95 = \
            get_k_data_for_cases(data_file_root + group_15hr, cases, this_k)
    thus_ulim_95 = this_P + 2 * this_gauss
    plt.semilogy(modes_removed, this_ulim_95, lw=2.,
                 color=colors[ii % len(colors)],
                 linestyle=line_styles[ii % len(line_styles)])
    plt.semilogy(modes_removed, this_ulim_95 - this_P, lw=0.7,
                 color=colors[ii % len(colors)],
                 linestyle=line_styles[ii % len(line_styles)])
plt.ylabel(r'$\Delta^2 ({\rm K}^2)$')
plt.xlabel('modes removed')



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



def get_mutual_upper_limit(fname1, fname2):
    data1 = np.loadtxt(fname1)
    k1 = data1[:,0]
    P1 = data1[:,1]
    std1 = data1[:,2]
    gauss1 = data1[:,3]
    
    data2 = np.loadtxt(fname2)
    k1 = data1[:,0]
    P1 = data1[:,1]
    std1 = data1[:,2]
    gauss1 = data1[:,3]

    


def get_upper_limit(fname, ):
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




#### VERY OLD ####

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

