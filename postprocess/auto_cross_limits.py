"""This script makes plots for putting limits on the HI powerspectrum from
both the cross power and the auto power. It is also adapted to put an upper
limit on shot noise, from the auto power only. It starts in 2D and combines
posteriors down to 1D.

See the functions 'analyze_2d' and 'main_plot' for the bulk of the 
computation. Variables 'auto_root' and 'cross_root' point to the 
directories in which 2D auto and cross power data can be found. Parkes uses
full 2D power data for the analysis, while GBT only uses the Omega_HI*b_HI
result to scale the simulation which bounds the intervals from below;
specify that result in 'CMEAN' and 'CERR'.
"""

import numpy as np
from scipy import stats, interpolate, integrate, optimize
import matplotlib
import matplotlib.pyplot as plt
from mkpower import check_ps as th
import sys
import h5py
from core import algebra as al
from mkpower import ps_summary as pss
# Need this got get rid of 'Type 3 fonts' that MNRAS didn't like.
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

main_plot = False

#### Parameters ####
T_b = 0.0000521424
b_HI = 0.6
b_opt = 1.
# Number of degrees of freedom for pair variance error bars.
NDEF = 5
# Measurements from the cross power paper.
CMEAN = 0.43
CERR = 0.07
LOW = 1e-10  # Small number in units K**2
# k ranges for joint constraints.
KLOW = 0.04
KHIGH = 10.0
GAUSS_ERROR = True
DECORR = True
# Flags for which prior to use.
FORE_FLAT =  True
OMEGA_FLAT = False

data_file_root = '/scratch2/p/pen/nluciw/parkes/analysis_IM/'

#root_15hr = ('GBT_15hr_map_oldcalpolstack_x_GBT_15hr_map_oldcalpolstack'
#                '_onesided_ixi_conv1p4_clipnoise_bfwindow')
#auto_root = data_file_root + '256_thru_260_transfer_select/'
auto_root = data_file_root + 'field_ra165_transfer_select/'
cross_root = data_file_root + 'full446_cros_ps_10mode/'

# Need to divide by root 6.
#group_1hr = root_1hr + '_noiseweight/'
#group_15hr = root_15hr + '_noiseweight/'

auto_file = "comp_varweight_auto_ps_1-20mode_1dpow.txt"
sim_group_1 =("full261_cros_ps_15mode/")
sim_case = "cros_si_15mode_1dpow.txt"

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
    if OMEGA_FLAT:
        prior = np.ones(len(omega))  # Flat prior on Omega, 1/sqrt on power.
    else:
        prior = omega  # Flat prior on power.
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

def signal_likelihood_k_bin(auto_sim, cross_sim, auto_pow1, auto_err1,
                            cross_pow, cross_err, dm_sim, signal_pow):
    """Get signal likelihood function in a k bin."""
    # Where to calculate the probability densities.
    #signal_pow_min = auto_sim * cross_err 
    # XXX
    # signal_pow_min = sim_pow * CERR * 0.01
    #signal_pow_max = np.abs(auto_pow1) + 7 * np.abs(auto_err1)
    #signal_pow = np.arange(signal_pow_min, signal_pow_max, signal_pow_min)
    # Get the pdf from the cross-power.
    p_cross_int_r = p_cross_given_auto_int_r(signal_pow, dm_sim, cross_pow, cross_err)
    # Get the pdfs, from the auto-powers.
    p_auto_int_f_1 = p_auto_given_signal_int_f(signal_pow, auto_pow1, auto_err1)
    #p_auto_int_f_2 = p_auto_given_signal_int_f(signal_pow, auto_pow2, auto_err2)
    # Calculate prior.
    if OMEGA_FLAT:
        prior = 1. /  np.sqrt(signal_pow)  # Flat prior on Omega.
    else:
        prior = np.ones(len(signal_pow))  # Flat prior on power.
    # Combine.
    likelihood = p_cross_int_r * p_auto_int_f_1 * prior
    #likelihood = p_auto_int_f_1 * prior
    likelihood = normalize_pdf(signal_pow, likelihood)
#    plt.plot(signal_pow, normalize_pdf(signal_pow,p_cross_int_r), 'k--', color='r')
#    plt.plot(signal_pow, normalize_pdf(signal_pow,p_auto_int_f_1), 'k-.', color='b')
#    plt.plot(signal_pow, likelihood, 'k:', color='g')
 #   plt.yscale('log')
 #   plt.ylim([1.e1, 1.e3])
#    plt.savefig('check.eps')
#    print 'Saved check, cross=', cross_pow
    return signal_pow, likelihood

def get_sim(fname):
    # Load the simulation data.
    data = np.loadtxt(fname)
    # Get the simulation curves.
    #mask = np.isfinite(data[:,4])
    k = data[:,0]
    P = data[:,1]
    return k, P

def get_cross(fname):
    data = np.loadtxt(fname)
    # Get the data curves.
    k = data[:,0]
    P = data[:,1]
    if DECORR:
        err = data[:,2]
        return k, P, err
    else:
        #std_err = data[mask,2] / np.sqrt(NDEF + 1)
        std_err = data[:,2]
        # XXX
        gauss_err = data[:,2]
       	if GAUSS_ERROR:
            return k, P, gauss_err
        else:
            return k, P, std_err

def get_auto(fname):
    data = np.loadtxt(fname)
    # Get the data curves.
    k = data[:,0]
    P = data[:,1]
    if DECORR:
        err = data[:,2]
        return k, P, err
    else:
        #std_err = data[mask,2] / np.sqrt(NDEF + 1)
        std_err = data[:,2]
        # XXX
        gauss_err = data[:,2]
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

def p_cross_given_auto_int_r(signal_pow, sim_pow, cross_pow, cross_err):
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
#    cross_pow = np.sqrt(cross_pow)
#    cross_err = (1/2.)*cross_pow**(-1./2.)*(cross_err)

    # If using actual cross power data
  #  amp = np.sqrt(signal_pow)*np.sqrt(sim_pow)
    # Else
    amp = np.sqrt(signal_pow)
    p_cross = (stats.norm.cdf((amp - cross_mean)/cross_err)
                - stats.norm.cdf(-cross_mean/cross_err))
#    np.save('amp', amp)
#    np.save('signal_pow', signal_pow)
#    np.save('p_cross', p_cross)
#    np.save('cross_err', cross_err)
    p_cross *= cross_err / amp
    return p_cross

def p_auto_given_signal_int_f(signal_pow, auto_mean, auto_err):
    """Get probability of obtaining auto-power measurement, given auto-power.
    
    Gives the probability of obtaining an auto-power measurement
    as a function of underlying signal-power amplitude, then integrates out
    the unknown foreground power.

    Function is unnormalized but has the right dependance on `signal_pow`.
    """

    if not FORE_FLAT:
        # See Apr 29 of Kiyo's notes.
        # This prior doesn't work since the integral diverges at f->0.
        # Get the integration domain.
        min_f = auto_err / 10.
        max_f = auto_mean + 4 * auto_err
        forground_pow = np.arange(min_f, max_f, min_f)
        # Simpson's rule integration requires an even number of samples.
        if not (len(forground_pow) % 2):
            foreground_pow = forground_pow[:-1]
        # Data pdf depends on s + f, so start with this sum (for every signal
        # input).
        integrand = signal_pow[...,None] + forground_pow
        # Now form the data pdf.
        if GAUSS_ERROR:
            integrand = stats.norm.pdf((integrand - auto_mean) / auto_err)
        else:
            integrand = stats.t.pdf((integrand - auto_mean) / auto_err, NDEF)
        # Add the prior.
        print auto_mean, auto_err
        #integrand[forground_pow > 10] = 0
        integrand *= 1. / (1 + (forground_pow / auto_err / 2)**2)
        #integrand *= 1. / (forground_pow)**0.9
        #integrand[forground_pow < CMEAN**2 /50] =0
        # Integrate.
        p_auto = integrate.simps(integrand, dx=min_f, axis=-1)
    else:
        # See Feb 26, 2013 of Kiyo's notes for derivation.
        if GAUSS_ERROR:
            p_auto = stats.norm.sf((signal_pow - auto_mean) / auto_err)
        else:
            p_auto = stats.t.sf((signal_pow - auto_mean) / auto_err, NDEF)
    return p_auto

def normalize_pdf(lam, pdf):
    """Normalize a function to integrate to unity."""
    norm = integrate.simps(pdf, lam)
    pdf /= norm 
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

def auto_pow_plot():
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
    #  Load all the data.
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

def analyze_2d(shell=True):

    #Load all data
    dm_file = h5py.File(data_file_root + 
                  'full267_cros_ps_15mode/ps_result.hd5', 'r')
    dm_pow, k_num_dm, k_p_edges, k_v_edges =\
        pss.load_power_spectrum('cros_si_15mode_2dpow', dm_file)
    cros_file = h5py.File(data_file_root + 
                  'full446_cros_ps_10mode/avg_result.hd5', 'r')
    cros_pow = al.make_vect(al.load_h5(cros_file, 'cros_ps_10mode_2davg'))
    cros_std = al.make_vect(al.load_h5(cros_file, 'cros_std_10mode_2davg'))

    auto_file = h5py.File(data_file_root + 
                  'field_ra165_transfer_select/avg_result.hd5', 'r')
    auto_pow = al.make_vect(al.load_h5(auto_file, 'ps_2d'))
    auto_std = al.make_vect(al.load_h5(auto_file, 'std_2d'))

    #GBT
#    auto_pow = al.make_vect(al.load(data_file_root + 'tzuching_gbt/power_2d.npy'))
#    auto_std = al.make_vect(al.load(data_file_root + 'tzuching_gbt/err_2d.npy'))

    #Get centres of k bins along each dimension in the 2D pow
    k_p_centre, k_v_centre = pss.get_2d_k_bin_centre(dm_pow)
    k_edges, k_v_edges = pss.get_2d_k_bin_edges(dm_pow)

    #GBT
#    dm_pow = np.ones_like(auto_pow)
#    cros_pow = np.ones_like(auto_pow)
#    cros_std = np.ones_like(auto_pow)
#    k_p_centre = np.load(data_file_root + 'tzuching_gbt/k_cent.npy')
#    k_v_centre = np.load(data_file_root + 'tzuching_gbt/k_cent.npy')
#    k_edges = np.load(data_file_root + 'tzuching_gbt/k_left.npy')
#    k_edges = np.logspace(np.log10(0.0045), np.log10(4.0), 44)

    #Compute the centre of each pixel (ie |k|)
    k_centre = np.sqrt(k_v_centre[:, None]**2 + k_p_centre[None, :]**2)
    #Compute the angle to each pixel, from k_perp
    angle = (np.arccos(k_p_centre[:,None] / k_centre))
    angle = angle.flatten()

    #Convert from dimensionless \Delta to P(k), if desired
    auto_pow *= k_centre**(-3) * 2. * np.pi**2
    auto_std *= k_centre**(-3) * 2. * np.pi**2
    cros_pow *= k_centre**(-3) * 2. * np.pi**2
    cros_std *= k_centre**(-3) * 2. * np.pi**2
    dm_pow *= k_centre**(-3) * 2. * np.pi**2

    #PARKES SPECIFIC: use only first half of data
    #auto_pow = auto_pow[:,:-16]
    #auto_std = auto_std[:,:-16]
    #cros_pow = cros_pow[:,:-16]
    #cros_std = cros_std[:,:-16]
    #dm_pow = dm_pow[:,:-16]
    #k_centre = k_centre[:,:-16]

    #Compute array of only unique |k|
    k_centre = k_centre.flatten()
    k_unique, indices = np.unique(k_centre, return_inverse=True)
    
    #Init lists of power in the 2D arcs
    dm_pow_arcs = []
    cros_pow_arcs = []
    cros_std_arcs = []
    auto_pow_arcs = []
    auto_std_arcs = []
    auto_sim_arcs = []
    k_vals = []
    angles = []

    # Each iteration of the loop selects the powers corresponding to k
    # values in the bin with right edge k_edges[ii+1], appending to the
    # appropriate 'arcs' list
    for ii in range(len(k_edges)-1):

        #mask to select pixels with |k| of bin 
        mask = (k_centre < k_edges[ii+1]) & (k_centre > k_edges[ii])
        #order the data from lowest angle to k_perp to highest angle
        order = np.argsort(angle[mask])

        #use mask to select pixels in a given arc, and order them
        auto_rad = auto_pow.flatten()[mask][order]
        auto_std_rad = auto_std.flatten()[mask][order][auto_rad!=0]

        dm_at_rad = dm_pow.flatten()[mask][order][auto_rad!=0]
        dm_pow_arcs.append(dm_at_rad)

        angles.append(angle[mask][order][auto_rad!=0])
        k = k_centre[mask][order][auto_rad!=0]
        k_vals.append(k)

        cros_rad = cros_pow.flatten()[mask][order][auto_rad!=0]

        cros_std_rad = cros_std.flatten()[mask][order][auto_rad!=0]
        cros_pow_arcs.append(cros_rad)
        cros_std_arcs.append(cros_std_rad)

        auto_sim= th.pk_th_mono(k, T_b**2, b_HI, b_HI)
        auto_sim_arcs.append(auto_sim)
        auto_pow_arcs.append(auto_rad[auto_rad!=0])
        auto_std_arcs.append(auto_std_rad)
        if ii==len(k_edges)-4:
            print auto_std, auto_std_rad

    #Use the plotting routine to make plot. PARKES: take [6:] subset of 
    # all arrays except shell.
    main_plot(angles, auto_pow_arcs, auto_std_arcs,
                  auto_sim_arcs, cros_pow_arcs, 
                  cros_std_arcs, dm_pow_arcs, 
                  k_edges, shell)

def main_plot(angles, auto_pow, auto_err, auto_sim, cross_pow, cross_err,
              dm_pow, k, shell, shellnum=0):
    ## Main plot - need to turn into a function. ##
    intervals = [.95, .68]
    interval_cols = ['1.0', '0.8', '0.65']
    imK2 = 1.
    # Set up plot. Commented lines are for different plot styles; script isn't
    # currently portable.
    f = plt.figure()
    ax = plt.gca()
    #ax.set_yscale('log')
    ax.set_xscale('log')
    ylab = plt.ylabel(r'$p[ln(\Delta^2)]$')
#    xlab = plt.xlabel(r'$cos^{-1}(\frac{\| k \|}{k_\perp})$')
    xlab = plt.xlabel(r'$\Delta^2 ({\rm K}^2)$')
    for item in ([ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    #Initialize lists for quantities of each arc. Names may be incorrect:
    # name choice carried over from old code.
    lhoods = []
    sks = []
    kplot = []
    shot_post = 1

    #PARKES has issues with the last few bins, so use len(angles)-5.
    #Loop over arcs
    arcs  = range(len(angles)-5)
    for jj in arcs:
      if angles[jj].size:
        # Allocate memory of the confidence intervals.
        pows = []
        likelihoods = []
        pows_no_f = []
        likelihoods_no_f = []
        n_intervals = len(intervals)
        n_k = len(angles)
        ulim = np.zeros((n_intervals, n_k))
        llim = np.zeros((n_intervals, n_k))
        ulim_no_f = np.zeros((n_intervals, n_k))

        #These are tests for the choice of the range of the independent variable
        # in the inference analysis, some dynamically changing based on scale of
        # power in that bin. Currently not working well.
        #sig_pow_min = np.min(auto_sim[jj]) * np.mean(cross_err[jj])
        #sig_pow_min = np.max(1e-12, np.min(np.abs(auto_sim[jj])) * np.mean(cross_err[jj]))
        #sig_pow_max = np.mean(auto_pow[jj]) + 4*np.mean(auto_err[jj])

        #Compute array of signal values over which posteriors are computed
        sig_pow_max = 4e-2
        sig_pow_min = 1e-12
        signal_pow = np.logspace(np.log10(sig_pow_min), np.log10(sig_pow_max), num=1e6)
#        signal_pow = np.arange((sig_pow_min), (sig_pow_max), sig_pow_min)

        #This code is also used for calculating an upper limit on shot noise. 
        shot_n_min = 1e-11
        shot_n_max = 1e-3
        shot_noi = np.logspace(np.log10(shot_n_min), np.log10(shot_n_max), num=1e6)
        nshot_noi = np.flipud(shot_noi*(-1))
        shot_noi = np.append(nshot_noi, shot_noi)

        #Select the arc over which we compute the aggregate 
        ang = angles[jj]
        # Loop over angle and calculate the confidence for each bin. 
        for ii in range(len(ang)):
            print 'Angle ', ang[ii]
            # Data points for this angle within arc
            this_k = ang[ii]
            this_auto_sim = auto_sim[jj][ii]
            this_dm_sim = dm_pow[jj][ii]    
            this_auto_pow = auto_pow[jj][ii]
            this_cross_pow = cross_pow[jj][ii]
            this_auto_err = auto_err[jj][ii]
            this_cross_err = cross_err[jj][ii]

            #Not sure if next two lines are necessary
            pow = [0]
            likelihood = [0]

            #Calculate posterior for the angle. pow is the signal over which
            # posterior is computed
         #   pow, likelihood = signal_likelihood_k_bin(this_auto_sim, this_cross_sim,
         #                                         this_auto_pow, this_auto_err, 
         #                                         this_cross_pow, this_cross_err,
         #                                         this_dm_sim, signal_pow)

            # If only upper lim from auto
            likelihood = p_auto_given_signal_int_f(signal_pow, this_auto_pow, this_auto_err)
            pow = signal_pow

            #Posterior for shot noise
            hipass_shot = 5.9e-8
            hipass_err = 0.7e-8
            shot_post = (stats.norm.pdf(signal_pow,
                         loc=hipass_shot, scale=hipass_err))

            #Somewhat hacky way to avoid numerical issues. Function is normalized,
            # so should be fine.
          #  if np.median(shot_post) < 1e-50:
          #      shot_post *= 1e50

            #Avoid numerical issues, and append each angle to list
            if np.sum(np.isnan(likelihood)) < 1:
                pows.append(pow)
                likelihoods.append(likelihood)

        #One temp for each arc
        temp_sk = np.mean(np.array(pows), axis=0)
        print 'CHECK', temp_sk.shape, np.prod(np.array(likelihoods), axis=0).shape
        temp_l = normalize_pdf(temp_sk, np.prod(np.array(likelihoods), axis=0))

        #List of arcs' posteriors
        lhoods.append(temp_l)
        sks.append(temp_sk)
        kplot.append(k[jj])

        #Used to plot posteriors of angles within single arc. Set shell=True
        print jj
        if jj==1000:
            print 'plotted'
            plt.cla()
            plt.plot(temp_sk, temp_l, 'r--', lw=3)
            plt.xlim([1e-12, 5e-2])
            #plt.ylim([0, 0.7])
            #llm, ulm = get_conf_interval(temp_sk, temp_l)
            #plt.axvline(x=llm, ymin=0., ymax=1.0,
            #    color='r', ls='--', lw=2)
            #plt.axvline(x=ulm, ymin=0., ymax=1.0,
            #    color='g', ls='--', lw=2)
            #cdf = np.zeros(len(temp_sk))
            #cdf[1:] = integrate.cumtrapz(temp_l, temp_sk)
            #cdf /= cdf[-1]
            #plt.plot(temp_sk, cdf, 'b--', lw=3)
            plt.xscale('log')
            plt.yscale('log')
            np.save('temp_l', temp_l)
            np.save('temp_sk', temp_sk)
            plt.savefig('likelihood.eps')

            #Plot confidence intervals for signal in angles over single arc
            if shell:
                k = ang 
                n_k = len(k)
                ulim = np.zeros((n_intervals, n_k))
                llim = np.zeros((n_intervals, n_k))
                ulim_no_f = np.zeros((n_intervals, n_k))
                for ll in range(len(likelihoods)):         
                    for mm in range(n_intervals):
                        llim[mm,ll], ulim[mm,ll] = get_conf_interval(pows[ll], likelihoods[ll],
                                                         intervals[mm])
                ax.set_yscale('log')
                ax.set_xscale('linear')
                print k_rad
                ax.text(0.2, 1e-4, r'$\| k \| \simeq %.2f$'%(k_rad))
                ylab = plt.ylabel(r'$\Delta^2 ({\rm K}^2)$')
                xlab = plt.xlabel(r'$cos^{-1}(\frac{k_\perp}{\| k \|})$')
                for item in ([ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(20)
                plt.xlim([0.0, 1.6])
                plt.ylim([6e-12, 2e-3])

    #After all arcs complete, clear plot (not sure if necessary)
    plt.cla()

    #If plotting 1D, set shell=False. Plotting stuff from here on.
    if not shell:
        #Get confidence intervals for each arc
        for ii in range(len(lhoods)):         
            for jj in range(n_intervals):
                llim[jj,ii], ulim[jj,ii] = get_conf_interval(sks[ii], lhoods[ii],
                                                         intervals[jj])
        #Plot settings
        ax.set_yscale('log')
        ax.set_xscale('log')
        ylab = plt.ylabel(r'$P(k) ({\rm K}^2)$')
        xlab = plt.xlabel(r'$k\, (h\, {\rm Mpc}^{-1})$')
        for item in ([ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
        plt.ylim([9e-10, 2e-3]) 
        plt.xlim([4e-2, 2e0]) 

    # Get left edges for k-bins
   # k_left = np.empty_like(k[:len(lhoods)])
   # k_left[0] = k[0] - (k[1] - k[0]) / 2
   # for ii in range(1,len(k[:len(lhoods)])):
   #     k_left[ii] = 2*k[ii - 1] - k_left[ii - 1]
    k_left = kplot
    print 'k_left', k_left

    #Get widths for k-bins
    dk = np.empty_like(k[:len(lhoods)])
    dk[:-1] = np.diff(k_left)
    dk[-1] = dk[-2]

    #Plot the limits.
    bottom = llim[0,:len(lhoods)]
    for jj in range(n_intervals):
        plt.bar(k_left, imK2*(ulim[jj,:len(lhoods)] - bottom), width=dk, bottom=imK2*bottom,
                color=interval_cols[jj + 1], linewidth=0.5)
    for jj in range(n_intervals - 1, 0, -1):
        plt.bar(k_left, imK2*(llim[jj,:len(lhoods)] - bottom), width=dk,
                bottom=imK2*bottom, color=interval_cols[jj], linewidth=0.5)

    #Plot any desired lines here.
#    cross_line = plt.plot(k, th.pk_th_mono(k, T_b**2, b_HI, b_HI)*k**(-3)*2*np.pi**2., 'g--', lw=2.)
#    cross_line = plt.plot(k[cross_pow > 0], cross_pow[cross_pow > 0] * T_b,
#                      'r--', lw=2.)
    auto_line = plt.plot(k, np.ones_like(k)*5.9e-8, 'b--', lw=2.)

    #Tidy up plot
    #xticks = np.arange(1, 10, 2)
    #xticks = list(0.01 * xticks) + list(0.1 * xticks) + list(xticks)
    #plt.xticks(xticks, xticks)
    f.subplots_adjust(top=0.7)
    #leg = plt.legend([cross_line[0]], 
    #                 [r'$\Omega_{\rm HI}b_{\rm HI}=0.43 \times 10^{-3}$'],
    #                 'lower right', prop={'size' : 16, 'family' : 'serif'})

    #Save fig of allowed signal in 1D 
    print 'Saving fig  ', data_file_root + 'tzuching_gbt/allowed_signal.eps'
    plt.savefig('allowed_signal.eps', bbox_inches='tight')
#bbox_extra_artists=(xlab, ylab)
    plt.cla()
    return

    #Shot noise stuff from here on. Irrelevant otherwise.
    confidence = 0.95
    vol_2df = 3e7
    gal_n = 2.5e5
    Tb = 5.2e-5
#    shot_post = shot_post[shot_noi>0]
    shot_noi = shot_noi[shot_noi>0]
    shot_post = normalize_pdf(shot_noi, shot_post)
    cd = np.zeros_like(shot_noi)
    cd[1:] = integrate.cumtrapz(shot_post, shot_noi)
    cd /= cd[-1]
    sn_int = interpolate.interp1d(cd, shot_noi)
    uplim = sn_int(confidence)
    lolim = sn_int(0.05)
    plt.plot(shot_noi, shot_post*shot_noi, 'k', lw=2)
    plt.ylabel(r'$p[ln({T_\mathrm{b}}^2 / \overline{n})]$')
    plt.xlabel(r'${T_\mathrm{b}}^2 / \overline{n}\; [\mathrm{K}^2\, \mathrm{Mpc}^3\, h^{-3}]$')
    plt.xlim([1e-11, 1e-3])
    plt.xscale('log')
    plt.ylim([0, 0.5])
    plt.axvline(x=uplim, ymin=0., ymax=np.max(shot_post*shot_noi),
                color='r', ls='--', lw=2)
    #plt.annotate('${T_\mathrm{b}}^2 / \overline{n} =$%.2e \n'%uplim+
    #             '$\mathrm{HIPASS}\Rightarrow \overline{n} =$%.2e \n'%(Tb**2/uplim)+
    #             '$\mathrm{2dF}\Rightarrow {T_\mathrm{b}} =$%.2e'%(np.sqrt(uplim*gal_n/vol_2df)),
    #             xy=(uplim,np.max(shot_post*shot_noi)), xytext=(uplim*5, 0.35))
    plt.savefig(auto_root + 'shot_noise.eps')
    plt.cla()
    return
    cd = np.zeros_like(shot_noi)
    cd[1:] = integrate.cumtrapz(shot_post*(shot_noi**2), shot_noi**(-1))
    cd /= cd[-1]
    sn_int = interpolate.interp1d(cd, shot_noi**(-1))
    uplim = sn_int(confidence)
    lolim = sn_int(0.05)
    plt.plot(shot_noi**(-1), shot_post*(shot_noi**(2)), 'k', lw=2)
    plt.ylabel(r'$p[(\overline{n} / {T_\mathrm{b}}^2 )]$')
    plt.xlabel(r'$\overline{n} / {T_\mathrm{b}}^2\; [\mathrm{K}^{-2}\, \mathrm{Mpc}^{-3}\, h^{3}]$')
    #plt.xlim([1e-12, 1e-4])
    plt.ylim([0, 0.0001])
    plt.xscale('log')
    #plt.yscale('log')
    plt.axvline(x=lolim, ymin=0., ymax=np.max(shot_post*shot_noi),
                color='r', ls='--', lw=2)
    #plt.annotate('${T_\mathrm{b}}^2 / \overline{n} =$%.2e \n'%uplim+
    #             '$\mathrm{HIPASS}\Rightarrow \overline{n} =$%.2e \n'%(Tb**2/uplim)+
    #             '$\mathrm{2dF}\Rightarrow {T_\mathrm{b}} =$%.2e'%(np.sqrt(uplim*gal_n/vol_2df)),
    #             xy=(uplim,np.max(shot_post*shot_noi)), xytext=(uplim*5, 0.35))
    plt.savefig('shot_noise_inv.eps')

    print 'PLOT SAVED', 'uplowlim', uplim, lolim

if __name__=='__main__':
    
    analyze_2d(shell=False)
