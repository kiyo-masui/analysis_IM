import numpy as np
import scipy as sp
import scipy.interpolate as interp
from simulations import gaussianfield as gf
from utils import cosmology
import astropy.units as u
from astropy import constants as const
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib import ticker
from matplotlib.ticker import LogFormatterMathtext
import pylab
import core.algebra as al
from quadratic_products import pwrspec_estimator as ps
import h5py
import copy
import scipy.integrate as integrate
from mpi4py import MPI
import gc
rcdef = plt.rcParams.copy()

#ps_fname = '/scratch2/p/pen/andersoc/simulated_ps/class_planck_z2point05.npy'
#sim_out = '/scratch2/p/pen/andersoc/simulated_ps/gaussian_real/'

ps_fname = '/scratch2/p/pen/andersoc/simulated_ps/deltaz_point5/class_planck_z2point25.npy'
sim_out = '/scratch2/p/pen/andersoc/simulated_ps/deltaz_point5/gaussian_real/'

#freq in Hz
f_rest = (3.*10**8)/(157.7*10**-6)

def make_gaussian_realizations(sim_out = sim_out, subtract_slice_mean = False, sim_num=100, make_cont_too = False, amp = 0.48, side = 500):
    ps = np.load(ps_fname)
    ps_int = interp.interp1d(ps[:,0],ps[:,1], fill_value = 'extrapolate')

    #Making side Mpc/h cube wth 1 Mpc/h pixel size.
    rf = gf.RandomField([side,side,side])
    rf.powerspectrum = lambda x: ps_int((x[:,:,:,0]**2 + x[:,:,:,1]**2 + x[:,:,:,2]**2)**0.5)

    for num in range(sim_num):
        map = rf.getfield()
        if subtract_slice_mean:
            map -= np.mean(np.mean(map, axis=1),axis=1)[:,None,None]
        np.save(sim_out + 'sim' + str(num), map)
        if make_cont_too:
            cont_map = make_continuum_cube(map, a_c = amp)
            np.save(sim_out + 'sim' + str(num) + 'cont_ac48', cont_map)

def make_gaussian_realizations_parallel(sim_out = sim_out, subtract_slice_mean = False, sim_num=100, make_cont_too = False, amp = 0.48, side = 500, process_n = 0):
    ps = np.load(ps_fname)
    ps_int = interp.interp1d(ps[:,0],ps[:,1], fill_value = 'extrapolate')
    #Making side Mpc/h cube wth 1 Mpc/h pixel size.
    rf = gf.RandomField([side,side,side])
    rf.powerspectrum = lambda x: ps_int((x[:,:,:,0]**2 + x[:,:,:,1]**2 + x[:,:,:,2]**2)**0.5)

    map = rf.getfield()
    if subtract_slice_mean:
        map -= np.mean(np.mean(map, axis=1),axis=1)[:,None,None]
    np.save(sim_out + 'sim' + str(process_n), map)
    cont_map = make_continuum_cube(map, a_c = amp)
    np.save(sim_out + 'sim' + str(process_n) + 'cont_ac48', cont_map)

def dust_intensity(freq, temp=27, amp=1, beta = 1.5, f_crit = 3.*10**12):
    #Use SI units,freq in Hz and temp in K.
    h = np.double(const.h)
    k_B = np.double(const.k_B)
    c = np.double(const.c)
    #exp = lambda x: np.exp(h * x/(k_B*temp)) - 1
    #pref = lambda x: (amp*2*h*x**(3+beta))/(c**2)
    dust = lambda x: ((amp*2*h*x**(3+beta))/(c**2))/(np.exp(h * x/(k_B*temp)) - 1)
    #int = pref/exp
    #ans = np.piecewise(freq, [freq<=f_crit, freq> f_crit], [dust(freq), dust(f_crit)*(freq**(-2)/f_crit**(-2))])
    ans = dust(freq)
    return ans

def dust_intensity_piecewise(freq, temp=27, amp=1, beta = 1.5, f_crit = 3.*10**12):
    #Use SI units,freq in Hz and temp in K.
    h = np.double(const.h)
    k_B = np.double(const.k_B)
    c = np.double(const.c)
    #exp = lambda x: np.exp(h * x/(k_B*temp)) - 1
    #pref = lambda x: (amp*2*h*x**(3+beta))/(c**2)
    dust = lambda x: ((amp*2*h*x**(3+beta))/(c**2))/(np.exp(h * x/(k_B*temp)) - 1)
    #int = pref/exp
    ans = np.piecewise(freq, [freq<=f_crit, freq> f_crit], [dust(freq), dust(f_crit)*(freq**(-2)/f_crit**(-2))])
    #ans = dust(freq)
    return ans

def displacement_to_z(x, z_center):
    cosmo = cosmology.Cosmology()
    #Displacement, x, from center of map, should be in Mpc/h
    z = z_center + 100*x*cosmo.H(z_center)/(cosmo.H()*np.double(const.c/1000.))
    return z

def z_to_freq(z, f_cent = f_rest):
    f = f_cent/(1 + z)
    return f

def dust_map_space(x_source, x_map, z_center, f_rest_line = f_rest, dust_temp = 27):
    #x_source is the displacement of the source from the map center in Mpc/h
    #x_map is the displacement from the center for the point where you're measuring the intensity
    z_source = displacement_to_z(x_source, z_center)
    f_map = z_to_freq(displacement_to_z(x_map, z_center), f_rest)
    #print z_source
    #print f_map
    dust = dust_intensity((1 + z_source)*f_map)*((1+z_source)**(-3))
    return dust

def dust_func_conv_approx(x, z_center, amp):
    #x is displacement in Mpc/h
    cosmo = cosmology.Cosmology()
    z = z_center
    dust = dust_intensity(f_rest - (f_rest/(1+z))*100*x*cosmo.H(z_center)/(cosmo.H()*np.double(const.c/1000.)))
    dust /= dust_intensity(f_rest)
    dust *= amp
    return dust

def make_continuum_cube(int_map, z_center=2.25, pix_size=1, a_c=0.48, f_rest_line = f_rest, dust_temp=27, beta=1.5):
    #a_c is the ratio of the specific intensity of the intensity map line, in the bandwidth centered on the map pixel, to the dust continuum specific intensity in that same bandwidth.
    #int_map is a 3d map, with the first dimension assumed to be the line-of-sight.
    #pixel_size is in Mpc/h.
    shape = int_map.shape
    cont_map = np.zeros(shape, np.double)
    x_inds = np.arange(shape[0])
    disp = (x_inds - shape[0]//2)*pix_size
    #print disp
    for x in x_inds:
        sed_x = dust_map_space(disp[x],disp, z_center=z_center, f_rest_line = f_rest_line, dust_temp = dust_temp)
        cont_map += int_map[x,:,:][None,:,:]*(1/a_c)*sed_x[:,None,None]*(1/sed_x[x])
    return cont_map

def make_2dps_parallel(root, n_sims, subtract_slice_mean = True, reverse_order = False, bin_zeros = False, n = 0):
    int = np.load(root + 'sim' + str(n) + '.npy')
    int = al.make_vect(int)
    int.axes = ('x','y','z')
    int.info['x_centre']=0
    int.info['y_centre']=0
    int.info['z_centre']=0
    int.info['x_delta']=1.
    int.info['y_delta']=1.
    int.info['z_delta']=1.
    cont = np.load(root + 'sim' + str(n) + 'cont_ac48.npy')
    cont = al.make_vect(cont)
    cont.axes = int.axes
    cont.info = int.info
    if subtract_slice_mean:
        int -= np.mean(np.mean(int,axis=1),axis=1)[:,None,None]
        cont -= np.mean(np.mean(cont,axis=1),axis=1)[:,None,None]
    int_power = ps.calculate_xspec(copy.deepcopy(int), copy.deepcopy(int), np.ones(int.shape), np.ones(int.shape), unitless=False)
    gc.collect()
    if not reverse_order:
        int_cont_power = ps.calculate_xspec(copy.deepcopy(int), copy.deepcopy(cont), np.ones(int.shape), np.ones(cont.shape), unitless=False)
        gc.collect()
        int_copy = copy.deepcopy(int)
        cont_and_int = int_copy + copy.deepcopy(cont)
        cross_power = ps.calculate_xspec(int_copy, cont_and_int, np.ones(int.shape), np.ones(cont.shape), unitless=False)
        del(int_copy)
        del(cont_and_int)
        gc.collect()
    else:
        int_cont_power = ps.calculate_xspec(copy.deepcopy(cont), copy.deepcopy(int), np.ones(int.shape), np.ones(cont.shape), unitless=False)
        gc.collect()
        int_copy = copy.deepcopy(int)
        cont_and_int = int_copy + copy.deepcopy(cont)
        cross_power = ps.calculate_xspec(cont_and_int, int_copy, np.ones(int.shape), np.ones(cont.shape), unitless=False)
        del(int_copy)
        del(cont_and_int)
        gc.collect()
    return [int[0], int_cont_power[0], cross_power[0]]

def make_2dps_avgs(root, n_sims, subtract_slice_mean = True, reverse_order = False, bin_zeros = False):
    avg_int = {}
    avg_int_cont = {}
    avg_cross = {}
    avg_int['binavg'] = np.zeros((40,40),np.double)
    avg_int['counts_histo'] = np.zeros((40,40))
    avg_int['bin_cent'] = np.zeros((40,))
    avg_int_cont['binavg'] = np.zeros((40,40),np.double)
    avg_int_cont['counts_histo'] = np.zeros((40,40))
    avg_int_cont['bin_cent'] = np.zeros((40,))
    avg_cross['binavg'] = np.zeros((40,40),np.double)
    avg_cross['counts_histo'] = np.zeros((40,40))
    avg_cross['bin_cent'] = np.zeros((40,))
    for n in xrange(n_sims):
        int = np.load(root + 'sim' + str(n) + '.npy')
        int = al.make_vect(int)
        int.axes = ('x','y','z')
        int.info['x_centre']=0
        int.info['y_centre']=0
        int.info['z_centre']=0
        int.info['x_delta']=1.
        int.info['y_delta']=1.
        int.info['z_delta']=1.
        cont = np.load(root + 'sim' + str(n) + 'cont_ac48.npy')
        cont = al.make_vect(cont)
        cont.axes = int.axes
        cont.info = int.info
        if subtract_slice_mean:
            int -= np.mean(np.mean(int,axis=1),axis=1)[:,None,None]
            cont -= np.mean(np.mean(cont,axis=1),axis=1)[:,None,None]
        int_power = ps.calculate_xspec(copy.deepcopy(int), copy.deepcopy(int), np.ones(int.shape), np.ones(int.shape), unitless=False)
        if not reverse_order:
            int_cont_power = ps.calculate_xspec(copy.deepcopy(int), copy.deepcopy(cont), np.ones(int.shape), np.ones(cont.shape), unitless=False)
            cross_power = ps.calculate_xspec(copy.deepcopy(int), copy.deepcopy(int) + copy.deepcopy(cont), np.ones(int.shape), np.ones(cont.shape), unitless=False)
        else:
            int_cont_power = ps.calculate_xspec(copy.deepcopy(cont), copy.deepcopy(int), np.ones(int.shape), np.ones(cont.shape), unitless=False)
            cross_power = ps.calculate_xspec(copy.deepcopy(int) + copy.deepcopy(cont), copy.deepcopy(int), np.ones(int.shape), np.ones(cont.shape), unitless=False)
        avg_int['binavg'] += int_power[0]['binavg']
        avg_int['counts_histo'] += int_power[0]['counts_histo']
        avg_int_cont['binavg'] += int_cont_power[0]['binavg']
        avg_int_cont['counts_histo'] += int_cont_power[0]['counts_histo']
        avg_cross['binavg'] += cross_power[0]['binavg']
        avg_cross['counts_histo'] += cross_power[0]['counts_histo']
        print str(n)
    avg_int['binavg'] /= n_sims
    avg_int['counts_histo'] /= n_sims
    avg_int_cont['binavg'] /= n_sims
    avg_int_cont['counts_histo'] /= n_sims
    avg_cross['binavg'] /= n_sims
    avg_cross['counts_histo'] /= n_sims
    avg_int['bin_cent'] = int_power[0]['bin_x_center']
    avg_cross['bin_cent'] = cross_power[0]['bin_x_center']
    avg_int_cont['bin_cent'] = int_cont_power[0]['bin_x_center']
    return [avg_int, avg_int_cont, avg_cross]

def plot2Dpowerspectrumlognoneg(pwr, bins_par , bins_perp, output, normalize=False):
    pwrp = np.ma.array(pwr, mask=np.logical_or(pwr<=0, np.isnan(pwr)))
    pwrp = np.abs(pwrp)
    print np.min(pwrp)
    print np.max(pwrp)
    fig = plt.figure()
    imp = pylab.imshow(np.flipud(np.transpose(np.abs(pwrp))), extent=[bins_par[0], bins_par[-1], bins_perp[0], bins_perp[-1]], norm=LogNorm(vmin=np.min(pwrp), vmax=np.max(pwrp)))
    pylab.colorbar(imp, orientation='horizontal', label = r'$(Mpc/h$)^3')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$k_{\perp}$', fontsize=14, color='black')
    plt.ylabel(r'$k_{\parallel}$', fontsize=14, color='black')
    pylab.savefig(output+'.png')
    pylab.clf()

def plot_ps_kperp_avg(file_path, output, abs = False):
    with h5py.File(file_path) as dat:
        mean = np.nanmean(np.array(dat['binavg']), axis = 0)
        k = np.array(dat['bin_cent'])
        mask = np.isnan(mean)
        pylab.clf()
        pos_mask = np.logical_and(mean>0,~mask)
        neg_mask = np.logical_and(mean<0,~mask)
        if abs:
            mean = np.abs(mean)
        #pylab.loglog(k[~mask], mean[~mask], label = 'PS avged over k_perp vs k_par')
        print mean
        pylab.loglog(k[~mask], mean[~mask], label = 'Positive PS avged over k_perp vs k_par', linewidth=2.0)
        pylab.loglog(k[~mask], -mean[~mask], label = 'Negative PS avged over k_perp vs k_par', linewidth=2.0)
        pylab.legend(loc='upper right',prop={'size':8})
        pylab.savefig(output)

def plot_ps_3pairs_kperp_avg(data_root, output, abs = False, lower_limit = False):
    int = h5py.File(data_root + 'avg_intxint.hdf5', 'r')
    int_cont = h5py.File(data_root + 'avg_intxcont.hdf5','r')
    cross = h5py.File(data_root + 'avg_intxintpluscont.hdf5','r')
    labels = ['Galaxy cross line intensity', 'Galaxy cross correlated continuum', 'Galaxy cross microwave map']
    pylab.clf()
    for dat in enumerate([int,int_cont, cross]):
        mean = np.nanmean(np.array(dat[1]['binavg']), axis = 0)
        if dat[0] == 0:
            int_min = np.nanmin(mean)
        k = np.array(dat[1]['bin_cent'])
        if abs:
            mean = np.abs(mean)
        mask = np.isnan(mean)
        pylab.loglog(k[~mask], mean[~mask], label = labels[dat[0]], linewidth = 2.0 )
        pylab.legend(loc='upper right',prop={'size':8})
    if lower_limit:
        pylab.set_ylim(bottom = int_min/10.)
    pylab.savefig(output)
    int.close()
    int_cont.close()
    cross.close()

def plot_2dk(kx, ky, image, filename, logscale=False, cbar=True):
     plt.rcParams.update(rcdef)
     fig, ax1 = plt.subplots()

     plt.rc('font', family='sans-serif')
     plt.rc('font', serif='Helvetica Neue')
     plt.rc('text', usetex='false')

     plt.rcParams.update({'font.size': 12})
     # ax1.get_xaxis() and ax1.xaxis seem to be the same
     ax1.xaxis.set_tick_params(which = 'both', direction='out')
     ax1.yaxis.set_tick_params(which = 'both', direction='out')
     #ax1.xaxis.set_tick_params(which = 'both', width=2)
     #ax1.yaxis.set_tick_params(which = 'both', width=2)
     ax1.tick_params('both', length=4, width=1.5, which='major')
     ax1.tick_params('both', length=3, width=1, which='minor')
     plt.rcParams['xtick.direction'] = 'out'
     plt.rcParams['ytick.direction'] = 'out'

     for axis in ['top','bottom','left','right']:
         ax1.spines[axis].set_linewidth(2)

     plt.rcParams.update({'legend.fontsize': 12})

     if logscale:
         norm = LogNorm()
         format = LogFormatterMathtext()
     else:
         norm = None
         format = None

     im = plt.pcolor(kx, ky, image.T, linewidth=1, cmap='CMRmap', norm=norm)
     im.set_edgecolor('face')
     #plt.axvline(100, color='red')
     plt.xscale('log')
     plt.yscale('log')
     ax1.set_xlim([np.min(kx), np.max(kx)])
     ax1.set_ylim([np.min(ky), np.max(ky)])
     plt.xlabel(r"$k_\perp\,[h{\rm Mpc}^{-1}]$", fontsize=16)
     plt.ylabel(r"$k_\parallel\,[h{\rm Mpc}^{-1}]$", fontsize=16)

     # format=LogFormatterMathtext(),
     if cbar:
         cbar = plt.colorbar(im, orientation='vertical',
                             ax=ax1, shrink=.9, pad=0.02, aspect=10,
                             label=r"", format=format)

         cbar.ax.tick_params(labelsize=10)
         cbar.solids.set_edgecolor("face")

     plt.savefig(filename, format="eps", bbox_inches='tight')
     plt.show()

def plot_mode_num(original, cuts, bins, output):
    pylab.clf()
    pylab.loglog(bins, original**0.5, label = '$\sqrt(N) modes with no cuts$')
    pylab.loglog(bins, cuts**0.5, label = '$\sqrt(N) modes with no cuts$')
    pylab.savefig(output)
    pylab.clf()

def calc_mode_numbers(k_cuts, bins, radius_arr):
    #k_cuts is in shape of 3d power spectrum. It is 0 on cuts and 1 elsewhere.
    from utils import binning
    #arr_flat = k_cuts.flat
    #radius_flat = radius_arr.flat
    modes, frac = binning.bin_an_array(k_cuts, bins, radius_arr)
    return modes, frac

if __name__ == '__main__':
    fir_pow = integrate.quad(dust_intensity_piecewise, 2.449*10**12, 7.059*10**12)[0]
    planck_pow = integrate.quad(dust_intensity_piecewise, 1.268*10**12, 2.586*10**12)[0]
    #df = 3.05*(z_to_freq(displacement_to_z(0,2.05)) - z_to_freq(displacement_to_z(1,2.05)))
    df = 3.25*(z_to_freq(displacement_to_z(0,2.25)) - z_to_freq(displacement_to_z(1,2.25)))
    print df
    voxel_pow = integrate.quad(dust_intensity_piecewise, f_rest - df/2.,f_rest + df/2.)[0]
    #data1000 = '/scratch2/p/pen/andersoc/simulated_ps/gaussian_real/thousand_sims/'
    #make_gaussian_realizations(sim_out = data1000, sim_num=1000, make_cont_too = True)
    #ps_avgs = make_2dps_avgs(data1000, 1000, reverse_order = False)
    #int = h5py.File(data1000 + 'avg_intxint.hdf5')
    #int_cont = h5py.File(data1000 + 'avg_intxcont.hdf5')
    #cross = h5py.File(data1000 + 'avg_intxintpluscont.hdf5')
    #for key in ps_avgs[2].keys():
    #    int.create_dataset(key, data=ps_avgs[0][key])
    #    int_cont.create_dataset(key, data=ps_avgs[1][key])
    #    cross.create_dataset(key, data=ps_avgs[2][key])
    #int.close()
    #int_cont.close()
    #cross.close()

    #data_48planck = '/scratch2/p/pen/andersoc/simulated_ps/gaussian_real/planck48/'
    #data_fir_ratio_10000 = '/scratch2/p/pen/andersoc/simulated_ps/gaussian_real/fir_ratio_10000/'
    data_48planck = '/scratch2/p/pen/andersoc/simulated_ps/deltaz_point5/gaussian_real/planck48/'
    data_48planck_linedimmerby50 = '/scratch2/p/pen/andersoc/simulated_ps/deltaz_point5/gaussian_real/planck48_linedimmerby50/'
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    #make_gaussian_realizations_parallel(sim_out = data_48planck, sim_num=100, make_cont_too = True, amp = 0.48*planck_pow/voxel_pow, process_n = rank)
    #make_gaussian_realizations_parallel(sim_out = data_48planck_linedimmerby50, sim_num=100, make_cont_too = True, amp = 0.48*0.01*planck_pow/voxel_pow, process_n = rank)
    #make_gaussian_realizations(sim_out = data_fir_ratio_10000, sim_num=100, make_cont_too = True, amp = 0.0001*fir_pow/voxel_pow)
    #dir_array = [data_48planck, data_fir_ratio_10000]
    dir_array = [data_48planck, data_48planck_linedimmerby50]

    orig_data = sim_out


    for sim_out in dir_array:
        #ps_avgs = make_2dps_avgs(sim_out, 100, reverse_order = False)
        ps_n = make_2dps_parallel(sim_out, 100, reverse_order = False, n = rank)
        comm.Barrier()
        for ind in np.arange(100)[1:]:
            if rank == ind:
                comm.send(ps_n, dest=0, tag=13)
            if rank == 0:
                rec_data = comm.recv(source=ind, tag=13)
                for dic_ind in np.arange(len(ps_n)):
                    ps_n[dic_ind]['binavg'] += rec_data[dic_ind]['binavg']
                    ps_n[dic_ind]['counts_histo'] += rec_data[dic_ind]['counts_histo']
        comm.Barrier()
        ps_n[dic_ind]['binavg'] /= 100.
        ps_n[dic_ind]['counts_histo'] /= 100.
        ps_avgs = ps_n
        if rank == 0: 
            int = h5py.File(sim_out + 'avg_intxint.hdf5')
            #int_cont = h5py.File(sim_out + 'avg_intxcont.hdf5')
            cross = h5py.File(sim_out + 'avg_intxintpluscont.hdf5')
            int_cont = h5py.File(sim_out + 'avg_intxcont.hdf5')
            #int_cont = h5py.File(sim_out + 'avg_contxint.hdf5')
            #cross = h5py.File(sim_out + 'avg_intpluscontxint.hdf5')
            for key in ps_avgs[2].keys():
                if key not in int:
                    int.create_dataset(key, data=ps_avgs[0][key])
                if key not in int_cont:
                    int_cont.create_dataset(key, data=ps_avgs[1][key])
                if key not in cross:
                    cross.create_dataset(key, data=ps_avgs[2][key])
            int.close()
            int_cont.close()
            cross.close()
            int = h5py.File(sim_out + 'avg_intxint.hdf5', 'r')
            int_cont = h5py.File(sim_out + 'avg_intxcont.hdf5','r')
            cross = h5py.File(sim_out + 'avg_intxintpluscont.hdf5','r')
            cont_over_int = h5py.File(sim_out + 'avg_intxcont_divby_avg_intxint.hdf5')
            for key in int.keys():
                if key == 'binavg':
                    cont_over_int.create_dataset(key, data=np.array(int_cont[key])/np.array(int[key]))
                else:
                    cont_over_int.create_dataset(key, data=int[key])
            ps_avgs=[int, int_cont, cross]
            plot2Dpowerspectrumlognoneg(ps_avgs[0]['binavg'], ps_avgs[0]['bin_cent'], ps_avgs[0]['bin_cent'],sim_out + 'plots/' + 'avg_intxint')
            plot2Dpowerspectrumlognoneg(ps_avgs[1]['binavg'], ps_avgs[1]['bin_cent'], ps_avgs[1]['bin_cent'],sim_out + 'plots/' + 'avg_intxcont')
            plot2Dpowerspectrumlognoneg(ps_avgs[2]['binavg'], ps_avgs[2]['bin_cent'], ps_avgs[2]['bin_cent'],sim_out + 'plots/' + 'avg_intxintpluscont')
            plot2Dpowerspectrumlognoneg(np.array(ps_avgs[2]['binavg'])/np.array(ps_avgs[0]['binavg']), ps_avgs[2]['bin_cent'], ps_avgs[2]['bin_cent'],sim_out + 'plots/' + 'ratio_mean_cross_to_auto')
            plot2Dpowerspectrumlognoneg(ps_avgs[1]['binavg'], ps_avgs[1]['bin_cent'], ps_avgs[1]['bin_cent'],sim_out + 'plots/' + 'avg_contxint')
            plot2Dpowerspectrumlognoneg(ps_avgs[2]['binavg'], ps_avgs[2]['bin_cent'], ps_avgs[2]['bin_cent'],sim_out + 'plots/' + 'avg_intpluscontxint')
            plot2Dpowerspectrumlognoneg(np.array(ps_avgs[2]['binavg'])/np.array(ps_avgs[0]['binavg']), ps_avgs[2]['bin_cent'], ps_avgs[2]['bin_cent'],sim_out + 'plots/' + 'ratio_mean_crossrev_to_auto')
            plot2Dpowerspectrumlognoneg(np.array(ps_avgs[0]['binavg'])/(np.array(ps_avgs[0]['binavg'])+np.array(ps_avgs[1]['binavg'])), ps_avgs[2]['bin_cent'], ps_avgs[2]['bin_cent'],sim_out + 'plots/' + 'ratio_mean_auto_to_auto_plus_intxcont')
            plot2Dpowerspectrumlognoneg(np.array(ps_avgs[2]['binavg'])/(np.array(ps_avgs[0]['binavg'])+np.array(ps_avgs[1]['binavg'])), ps_avgs[2]['bin_cent'], ps_avgs[2]['bin_cent'],sim_out + 'plots/' + 'ratio_mean_cross_to_auto_plus_intxcont')
            plot2Dpowerspectrumlognoneg(cont_over_int['binavg'], cont_over_int['bin_cent'], cont_over_int['bin_cent'],sim_out + 'plots/' + 'ratio_mean_intxcont_to_auto')
            cont_over_int.close()
            plot_ps_kperp_avg(sim_out + 'avg_intxcont_divby_avg_intxint.hdf5', sim_out + 'plots/' + 'contxint_divby_intxint_avg_over_k_perp')
            #modes_file = h5py.File(orig_data + 'modes_lost_cut_2lowkpar.hdf5')
            #print modes_file.keys()
            #modes_orig = np.array(modes_file['modes_orig'])
            #modes_cut = modes_orig*np.array(modes_file['percent_modes_aftercuts'])
            #bins = np.array(ps_avgs[2]['bin_cent'])
            #plot_mode_num(modes_orig, modes_cut, bins, sim_out + 'plots/' + 'num_modes_kept')
            int.close()
            int_cont.close()
            cross.close()
            #modes_file.close()
            plot_ps_kperp_avg(sim_out + 'avg_intxint.hdf5', sim_out + 'plots/' + 'intxint_avg_over_k_perp')
            plot_ps_kperp_avg(sim_out + 'avg_intxcont.hdf5', sim_out + 'plots/' + 'intxcont_avg_over_k_perp')
            plot_ps_kperp_avg(sim_out + 'avg_intxintpluscont.hdf5', sim_out + 'plots/' + 'intxintpluscont_avg_over_k_perp')
            plot_ps_3pairs_kperp_avg(sim_out, sim_out + 'plots/' + 'avg_over_k_powerspectra', abs = True, lower_limit = False)
    #for ind in np.arange(100):
    #    int_map = np.load(sim_out + 'sim' + str(ind) + '.npy')
    #    cont_map = make_continuum_cube(int_map)
    #    np.save(sim_out + 'sim' + str(ind) + 'cont_ac48.npy', cont_map)

