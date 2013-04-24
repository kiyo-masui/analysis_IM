"""This script makes plots for putting limits on the HI powerspectrum from
both the cross power and the auto power."""

import numpy as np
from scipy import stats, interpolate, integrate, optimize
import matplotlib
import matplotlib.pyplot as plt

# Need this got get rid of 'Type 3 fonts' that MNRAS didn't like.
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

#### Parameters ####

# Number of degrees of freedom for pair variance error bars.
NDEF = 5
# Measurements from the cross power paper.
CMEAN = 0.43
CERR = 0.07
LOW = 1e-10  # Small number in units K**2
# k ranges for joint constraints.
KLOW = 0.075
KLOW = 0.12
KHIGH = 0.30
GAUSS_ERROR = True
DECORR = True

data_file_root = '/cita/h/home-2/eswitzer/code/analysis_IM/'

#root_15hr = ('GBT_15hr_map_oldcalpolstack_x_GBT_15hr_map_oldcalpolstack'
#                '_onesided_ixi_conv1p4_clipnoise_bfwindow')
root_15hr = ('pwrspec_plots/GBT_15hr_map_autopaper_x_GBT_15hr_map_autopaper'
                '_onesided_ixi_conv1p4_clipnoise_bfwindow_svdweighted')

#root_1hr = ('GBT_1hr_map_oldcalpolstack_x_GBT_1hr_map_oldcalpolstack'
#               '_onesided_ixi_conv1p4_clipnoise')
root_1hr = ('pwrspec_plots/GBT_1hr_map_oldcalpolstack_x_GBT_1hr_map'
            '_oldcalpolstack_onesided_ixi_conv1p4_clipnoise_bfwindow'
            '_svdweighted')

# Need to divide by root 6.
#group_1hr = root_1hr + '_noiseweight/'
#group_15hr = root_15hr + '_noiseweight/'

group_15hr = root_15hr + '_analysis/'
group_1hr = root_1hr + '_analysis/'
noise_group_15hr = root_15hr + '_self/'
noise_group_1hr = root_1hr + '_self/'

sim_group_15 = ("pwrspec_plots/GBT_15hr_map_autopaper_x_GBT_15hr_map_"
                "autopaper_onesided_ixi_conv1p4_clipnoise_bfwindow_"
                "svdweighted_puresignal/sim_nobeam/")
sim_group_1 =("pwrspec_plots/GBT_1hr_map_oldcalpolstack_x_GBT_1hr_map"
              "_oldcalpolstack_onesided_ixi_conv1p4_clipnoise_bfwindow"
              "_svdweighted_puresignal/sim_nobeam/")
sim_case = "power_1d_from_2d_0modes.dat"

modes_15hr = 10
modes_1hr = 30

#case_15hr = 'pk_15modes.dat'
#case_1hr = 'pk_45modes.dat'
file_type = 'pk_decorr_%dmodes.dat'
case_15hr = file_type % modes_15hr
case_1hr = file_type % modes_1hr
case_15hr_noise = 'pk_%dmodes.dat' % modes_15hr
case_1hr_noise = 'pk_%dmodes.dat' % modes_1hr
case_15hr_0 = 'pk_0modes.dat'
case_1hr_0 = 'pk_0modes.dat'


#### Functions ####

def get_likelihood_omega(fname1, fname2, fname_sim1, fname_sim2=None):
    """Given two auto power measurements and a simulation, get likelihood
    function for b * Omega.
    
    In this funtion by omega, I mean Omega * b.
    """
    
    if not fname_sim2:
        fname_sim2 = fname_sim1
    # Get x-axis.
    omega = np.arange(0.002, 2, 0.01)  # In units of 1e-3.
    omega_sq = omega**2
    n_omega = len(omega)
    # Load the data.
    data1, data2 = get_auto_match_k(fname1, fname2)
    k, auto_pow1, auto_err1 = data1
    k2, auto_pow2, auto_err2 = data2
    if not np.allclose(k, k2):
        raise RuntimeError('1 hour and 15 hour k bins not the same.')
    n_k = len(k)
    # Load the simulation (for omega=1e-3, which is 1 in the units we are
    # working in).
    k_sim1, sim_pow1 = get_sim(fname_sim1)
    sim_interp1 = interpolate.interp1d(k_sim1, sim_pow1)
    sim_pow1 = sim_interp(k)
    k_sim2, sim_pow2 = get_sim(fname_sim2)
    sim_interp2 = interpolate.interp1d(k_sim2, sim_pow2)
    sim_pow2 = sim_interp(k)
    # Calculate the auto power probabilites. The is the probabilty of obtaining
    # the full 'k'-set of auto-power measurements at the given omegas (just the
    # product of the probabilities at each 'k'), with an integral over a
    # foreground power in each bin.
    p_auto_int_f_1 = 1.
    p_auto_int_f_2 = 1.
    for ii in range(n_k):
        if k[ii] < KLOW or k[ii] > KHIGH:
            continue
        p_auto_int_f_1 *= p_auto_given_signal_int_f(omega_sq, \
                auto_pow1[ii] / sim_pow1[ii], auto_err1[ii] / sim_pow1[ii])
        p_auto_int_f_2 *= p_auto_given_signal_int_f(omega_sq, \
                auto_pow2[ii] / sim_pow2[ii], auto_err2[ii] / sim_pow2[ii])
        # Zero foreground, error bar only case.
        #p_auto_int_f_1 *= p_auto_given_signal_int_f(omega_sq, \
        #        1., auto_err1[ii] / sim_pow1[ii])
        #p_auto_int_f_2 *= p_auto_given_signal_int_f(omega_sq, \
        #        1., auto_err2[ii] / sim_pow2[ii])
    # Calculate cross power pdf.
    p_cross_int_r = p_cross_given_auto_int_r(omega_sq, 1.)
    # Calculate prior.
    prior = omega  # Flat prior on power.
    #prior = np.ones(len(omega))  # Flat prior on Omega, 1/sqrt on power.
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
    return float(llim), float(ulim)

def signal_likelihood_k_bin(sim_pow, auto_pow1, auto_err1, auto_pow2,
                            auto_err2):
    """Get signal likelihood function in a k bin."""
    # Where to calculate the probability densities.
    signal_pow_min = sim_pow * CERR * 0.05
    # XXX
    # signal_pow_min = sim_pow * CERR * 0.01
    signal_pow_max = min(auto_pow1 + 5 * auto_err1,
                         auto_pow2 + 5 * auto_err2)
    signal_pow = np.arange(signal_pow_min, signal_pow_max, signal_pow_min)
    # Get the pdf from the cross-power.
    p_cross_int_r = p_cross_given_auto_int_r(signal_pow, sim_pow)
    # Get the pdfs, from the auto-powers.
    p_auto_int_f_1 = p_auto_given_signal_int_f(signal_pow, auto_pow1, auto_err1)
    p_auto_int_f_2 = p_auto_given_signal_int_f(signal_pow, auto_pow2, auto_err2)
    # Calculate prior.
    prior = np.ones(len(signal_pow))  # Flat prior on power, favors high Omega.
    #prior = 1. /  np.sqrt(signal_pow)  # Flat prior on Omega.
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
    if DECORR:
        err = data[mask,2]
        return k, P, err
    else:
        #std_err = data[mask,2] / np.sqrt(NDEF + 1)
        std_err = data[mask,2]
        # XXX
        gauss_err = data[mask,3]
        if GAUSS_ERROR:
            return k, P, gauss_err
        else:
            return k, P, std_err

def get_noise(fname):
    data = np.loadtxt(fname)
    # Get the data curves.
    mask = np.isfinite(data[:,1])
    k = data[mask,0]
    P = data[mask,1]
    err = data[mask,2] 
    return k, P, err

def get_auto_match_k(fname1, fname2, fun=None):
    if fun is None:
        fun = get_auto
    data1 = list(fun(fname1))
    data2 = list(fun(fname2))
    k1 = data1[0]
    k2 = data2[0]
    mask1 = np.zeros(len(k1), dtype=bool)
    mask2 = np.zeros(len(k2), dtype=bool)
    for ii in range(len(k1)):
        for jj in range(len(k2)):
            if np.allclose(k1[ii], k2[jj]):
                mask1[ii] = True
                mask2[jj] = True
                break
    for ii in range(len(data1)):
        data1[ii] = data1[ii][mask1]
    for ii in range(len(data2)):
        data2[ii] = data2[ii][mask2]
    return tuple(data1), tuple(data2)

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
    if GAUSS_ERROR:
        p_auto = stats.norm.sf((signal_pow - auto_mean) / auto_err)
    else:
        p_auto = stats.t.sf((signal_pow - auto_mean) / auto_err, NDEF)
    return p_auto

def normalize_pdf(lam, pdf):
    """Normalize a function to integrate to unity."""
    norm = integrate.simps(pdf, lam)
    pdf = pdf / norm
    return pdf

def plot_log_pos_neg_errorbar(x, y, error, **kwargs):
    mask = y > 0
    if np.any(mask):
        err_lower = np.maximum(1e-10, y[mask] - error[mask])
        err_lower = y[mask] - err_lower
        line = plt.errorbar(x[mask], y[mask], [err_lower, error[mask]], lw=2,
                     **kwargs)
    mask = y < 0
    if np.any(mask):
        err_lower = np.maximum(1e-10, -y[mask] - error[mask])
        err_lower = -y[mask] - err_lower
        plt.errorbar(x[mask], -y[mask], [err_lower, error[mask]], lw=1, 
                     mfc='None', **kwargs)
    return line
  

#### Making the main plots ####

## Auto-power plot.
imK2 = 1e6
# Set up plot.
f = plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
ylab = plt.ylabel(r'$\Delta^2 ({\rm mK}^2)$')
xlab = plt.xlabel(r'$k ({\rm h/Mpc})$')
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# Load all the data.
data15, data1 = get_auto_match_k(data_file_root + group_15hr + case_15hr,
                                 data_file_root + group_1hr + case_1hr)
k, auto_pow15, auto_err15 = data15
k1, auto_pow1, auto_err1 = data1
# Plot the 15 hr field.
if GAUSS_ERROR:
    error_fact = 1
else:
    error_fact = 1.1  # 68% students t factor.
line15 = plot_log_pos_neg_errorbar(k, imK2 * auto_pow15, 
                                   error_fact * imK2 * auto_err15,
                                   color='g', linestyle='', marker='o')
line1 = plot_log_pos_neg_errorbar(1.02*k, imK2 * auto_pow1, 
                                  error_fact * imK2 * auto_err1,
                                  color='b', linestyle='', marker='s')
# Plot the cross-power r=1 line.
k_sim, sim_pow = get_sim(data_file_root + sim_group_1 + sim_case)
cross_line = plt.plot(k_sim, imK2*CMEAN**2*sim_pow, 'r--', lw=2.)
# Foreground PS.
data15, data1 = get_auto_match_k(data_file_root + group_15hr + case_15hr_0,
                                 data_file_root + group_1hr + case_1hr_0)
k_0, auto_pow15_0, tmp = data15
k1_0, auto_pow1_0, tmp = data1
#fg_line = plt.step(k, imK2*auto_pow1_0, 'b', where='mid', lw=2)
fg_line = plt.step(k, imK2*auto_pow15_0, 'g', where='mid', lw=2)
#plt.step(k, imK2*auto_pow1_0, 'b', where='mid', lw=2)
# Noise PS.
data15, data1 = get_auto_match_k(data_file_root + noise_group_15hr +
                                 case_15hr_noise,
                                 data_file_root + noise_group_1hr +
                                 case_1hr_noise,
                                 fun=get_noise)
k_n, noise_pow_15hr, tmp = data15
k_n1, noise_pow_1hr, tmp = data1
if not np.allclose(k, k_n) or not np.allclose(k, k_n1):
    raise RuntimeError()
#plt.step(k, imK2*noise_pow/np.sqrt(4), 'k:', where='mid', lw=2)
therm15_line = plt.plot(k, imK2*noise_pow_15hr/np.sqrt(4), 'g-.', lw=2)
therm1_line = plt.plot(k, imK2*noise_pow_1hr/np.sqrt(4), 'b:', lw=2)
# Formatting.
xticks = np.arange(1, 10, 1)
xticks = list(0.01 * xticks) + list(0.1 * xticks) + list(xticks)
plt.xticks(xticks, xticks)
plt.ylim([1e-3, 5e3])
plt.xlim([0.97*min(k), 0.7])
plt.savefig('auto_spectrum_no_leg.eps', bbox_extra_artists=(xlab, ylab),
            bbox_inches='tight')
# Legend.
leg = plt.legend([line15,
                  therm15_line[0],
                  fg_line[0],  
                  line1,
                  therm1_line[0],
                  cross_line[0]],
                 [r'deep auto-power',  
                  r'deep noise',
                  r'deep foreground', 
                  r'wide auto-power', 
                  r'wide noise',
                  r'$\Omega_{\rm HI}b_{\rm HI}=0.43 \times 10^{-3}$',], 
                 prop={'size' : 16, 'family' : 'serif'}, 
                 bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                 ncol=2, mode="expand", borderaxespad=0.)
plt.savefig('auto_spectrum.eps', bbox_extra_artists=(xlab, ylab, leg),
            bbox_inches='tight')



## Main plot - need to turn into a function. ##
#intervals = [.99, .95, .68]
#interval_cols = ['1.0', '0.8', '0.6', '0.4']
intervals = [.95, .68]
interval_cols = ['1.0', '0.8', '0.65']
imK2 = 1e6
# Set up plot.
f = plt.figure()
ax = plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')
ylab = plt.ylabel(r'$\Delta^2 ({\rm mK}^2)$')
xlab = plt.xlabel(r'$k ({\rm h/Mpc})$')
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# Load all the data.
k_sim, sim_pow = get_sim(data_file_root + sim_group_1 + sim_case)
sim_interp = interpolate.interp1d(k_sim, sim_pow)
data15, data1 = get_auto_match_k(data_file_root + group_15hr + case_15hr,
                                 data_file_root + group_1hr + case_1hr)
k, auto_pow15, auto_err15 = data15
k1, auto_pow1, auto_err1 = data1
if not np.allclose(k, k1):
    raise RuntimeError()
# Allocate memory of the confidance intervals.
n_intervals = len(intervals)
n_k = len(k)
ulim = np.zeros((n_intervals, n_k))
llim = np.zeros((n_intervals, n_k))
ulim_no_f = np.zeros((n_intervals, n_k))
# Loop over k and calculate the confidance for each bin.
for ii in range(n_k):
#for ii in [6]:
    # Data points for this bin.
    this_k = k[ii]
    this_sim_pow = sim_interp(k[ii])
    
    this_auto_pow15 = auto_pow15[ii]
    this_auto_pow1 = auto_pow1[ii]
    #this_auto_pow15 = this_sim_pow
    #this_auto_pow1 = this_sim_pow
    
    this_auto_err15 = auto_err15[ii]
    this_auto_err1 = auto_err1[ii]
    pow, likelihood = signal_likelihood_k_bin(this_sim_pow, this_auto_pow15,
                                              this_auto_err15, this_auto_pow1,
                                              this_auto_err1)
    for jj in range(n_intervals):
        llim[jj,ii], ulim[jj,ii] = get_conf_interval(pow, likelihood,
                                                     intervals[jj])
    # Limit if there where no foreground contamination.
    pow_no_f, likelihood_no_f = signal_likelihood_k_bin(this_sim_pow, 
                            CMEAN**2 * this_sim_pow, this_auto_err15,
                            CMEAN**2 * this_sim_pow, this_auto_err1)
    for jj in range(n_intervals):
        tmp, ulim_no_f[jj,ii] = get_conf_interval(pow_no_f, likelihood_no_f,
                                                  intervals[jj])
# Get k, left sides and widths for plotting.
k_left = np.empty_like(k)
k_left[0] = k[0] - (k[1] - k[0]) / 2
for ii in range(1,n_k):
    k_left[ii] = 2*k[ii - 1] - k_left[ii - 1]
dk = np.empty_like(k)
dk[:-1] = np.diff(k_left)
dk[-1] = dk[-2]
# Now plot the limits.
bottom = llim[0,:]
for jj in range(n_intervals):
    #plt.fill_between(k, imK2*ulim[jj,:], imK2*LOW, color=interval_cols[jj + 1])
    plt.bar(k_left, imK2*(ulim[jj,:] - bottom), width=dk, bottom=imK2*bottom,
            color=interval_cols[jj + 1], linewidth=0.5)
for jj in range(n_intervals - 1, 0, -1):
    #plt.fill_between(k, imK2*llim[jj,:], imK2*LOW, color=interval_cols[jj])
    plt.bar(k_left, imK2*(llim[jj,:] - bottom), width=dk,
            bottom=imK2*bottom, color=interval_cols[jj], linewidth=0.5)
    #plt.bar(k_left, imK2*llim[jj,:], width=dk,
    #        color=interval_cols[jj], linewidth=None)
# Plot the auto-power upper limits.
cross_line = plt.plot(k_sim, imK2*CMEAN**2*sim_pow, 'r--', lw=2.)
# Plot 95% upper limit if no foregrounds.
no_f_line = plt.step(k_left, imK2*ulim_no_f[0,:], lw=3, where='post',
                     color=(.4,0,0))
# Tidy up plot
xticks = np.arange(1, 10, 1)
xticks = list(0.01 * xticks) + list(0.1 * xticks) + list(xticks)
plt.xticks(xticks, xticks)
plt.ylim([imK2*2e-9, imK2*4e-7])
plt.xlim([0.97*min(k), 0.7])
f.subplots_adjust(top=0.7)
leg = plt.legend([no_f_line[0], cross_line[0]], 
                 [r'statistical limit',
                  r'$\Omega_{\rm HI}b_{\rm HI}=0.43 \times 10^{-3}$'],
                 'lower right', prop={'size' : 16, 'family' : 'serif'})
plt.savefig('allowed_signal.eps', bbox_extra_artists=(xlab, ylab),
            bbox_inches='tight')


# Full Omega*b likelihood plot.
likelihoods = get_likelihood_omega(data_file_root + group_15hr + case_15hr, 
                                   data_file_root + group_1hr + case_1hr, 
                                   data_file_root + sim_group_15 + sim_case,
                                   data_file_root + sim_group_1 + sim_case)
omega, likelihood, likelihood_auto_1, likelihood_auto_2, likelihood_cross \
        = likelihoods
# Put omega in full units.
omega *= 1e-3
# For logarithmic axes.
#m = 1
#lomega = omega
#xscale = 'linear'
#xlim = [0.0, 1.5]
median = get_conf_interval(omega, likelihood, 0.)[0]
interval = get_conf_interval(omega, likelihood, .68)
print "68% confidence interval:", interval
print "or:", (interval[0] + interval[1])/2, "\pm", (interval[1] - interval[0])/2
print "or:", median, "+", interval[1] - median, "-", median - interval[0]
# Log scale.
m = omega
lomega = np.log(omega)
xscale = 'log'
# Linear Scale
#m = 1
#lomega = omega
#xscale = 'linear'
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
ylab = plt.ylabel(r'$p\left[\rm{ln}'
                  r'(\Omega_{\rm HI} b_{\rm HI})\right]$')
xlab = plt.xlabel(r'$\Omega_{\rm HI} b_{\rm HI}$')
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
line_cross = plt.plot(omega, likelihood_cross, 'r--', lw=3)
line_auto1 = plt.plot(omega, likelihood_auto_1, 'g-.', lw=3)
line_auto2 = plt.plot(omega, likelihood_auto_2, 'b:', lw=3)
line = plt.plot(omega, likelihood, 'k', lw=2)
plt.legend(('cross-power', 'deep auto-power', 'wide auto-power', 'combined'),
           'upper left', prop={'size' : 18, 'family' : 'serif'})
plt.xlim(xlim)
plt.savefig('likelihood.eps', bbox_extra_artists=(xlab, ylab),
            bbox_inches='tight')


plt.show()

# Full likelihood table vs n_modes.
modes = range(5, 90, 5)
print "n_modes, 15 hr, 1 hr"
for mode in modes:
    likelihoods = get_likelihood_omega(data_file_root + group_15hr
                                         + file_type % mode, 
                                         data_file_root + group_1hr
                                         + file_type % 0,
                                         data_file_root + sim_group_15 + sim_case,
                                         data_file_root + sim_group_1 + sim_case)
    omega, likelihood, likelihood_auto_1, likelihood_auto_2, likelihood_cross \
            = likelihoods
    ubound15 = get_conf_interval(omega, likelihood, .68)[1]
    likelihoods = get_likelihood_omega(data_file_root + group_15hr
                                         + file_type % 0, 
                                         data_file_root + group_1hr
                                         + file_type % mode,
                                         data_file_root + sim_group_15 + sim_case,
                                         data_file_root + sim_group_1 + sim_case)
    omega, likelihood, likelihood_auto_1, likelihood_auto_2, likelihood_cross \
            = likelihoods
    ubound1 = get_conf_interval(omega, likelihood, .68)[1]
    print "%02d" % mode, "    %4.2f" % ubound15, "    %4.2f" % ubound1


