import core.algebra as al
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
from itertools import combinations
import pylab
from matplotlib.colors import LogNorm
from matplotlib import ticker
import math

def ploteigenvalues(file_name, output, key = '21cm_only', x_max = 0):
    file = h5py.File(file_name, 'r')
    vals = file['svd_vals'][key]
    vals = np.log10((vals/vals[0])**0.5)
    if x_max == 0:
        pylab.plot(range(len(vals)), vals, 'b',label='Sqrt Eigenvalue')
    else: 
	pylab.plot(range(x_max), vals[0:x_max], 'b', label='Sqrt Eigenvalue')
    pylab.xlabel('Mode number')
    pylab.ylabel('log(sqrt(eigenvalue))')
    pylab.savefig(output+'.png')
    pylab.clf()

#Switzer_replication_plot

def ploteigenvalues(file_name, output, key = '21cm_only', x_max = 0):
    file = h5py.File(file_name, 'r')
    vals = file['svd_vals'][key]
    vals = (vals/vals[0])**0.5
    if x_max == 0:
        pylab.loglog(np.arange(len(vals)) + 1, vals, 'b',label='Sqrt Eigenvalue', basex=10)
    else:
        pylab.plot(range(x_max), vals[0:x_max], 'b', label='Sqrt Eigenvalue')
    pylab.xlabel('Mode number')
    pylab.ylabel('Map-space amplitude')
    pylab.savefig(output+'.png')
    pylab.clf()

def ploteigenvalues_multiple(file_name1,file_name2,  output, key = '21cm_only', x_max = 0):
    file1 = h5py.File(file_name1, 'r')
    file2 = h5py.File(file_name2, 'r')
    vals1 = file1['svd_vals'][key]
    vals2 = file2['svd_vals'][key]
    vals1 = (vals1/vals1[0])**0.5
    vals2 = (vals2/vals2[0])**0.5
    if x_max == 0:
#        pylab.loglog(np.arange(len(vals1)) + 1, vals1, 'b',label='Gaussian Beam', basex=10)
        pylab.loglog(np.arange(len(vals1)) + 1, vals1, 'b',label='Airy Beam ($\sigma_{BP}$ = 15%)', basex=10)
        pylab.loglog(np.arange(len(vals2)) + 1, vals2, 'r',label='Airy Beam', basex=10)
    else:
        pylab.plot(range(x_max), vals[0:x_max], 'b', label='Sqrt Eigenvalue')
    pylab.xlabel('Mode number')
    pylab.ylabel('Map-space amplitude')
    pylab.legend(loc='lower right',prop={'size':8})
    pylab.savefig(output+'.png')
    pylab.clf()

def ploteigenvaluesnoise(file_name_noise_map, file_name_noise_only, output, key = '21cm_only', x_max = 0):
    file1 = h5py.File(file_name_noise_map, 'r')
    file2 = h5py.File(file_name_noise_only, 'r')
    vals1 = file1['svd_vals'][key]
    vals2 = file2['svd_vals'][key]
    norm = vals1[0]
    vals1 = np.log10((vals1/norm)**0.5)
    vals2 = np.log10((vals2/norm)**0.5)
    if x_max == 0:
        pylab.plot(range(256), vals1, 'b',label='Sqrt Eigenvalue')
        pylab.plot(range(256), vals2, 'g--',label='Sqrt Eigenvalue')
    else:
        pylab.plot(range(x_max), vals1[0:x_max], 'b', label='Sqrt Eigenvalue')
        pylab.plot(range(x_max), vals2[0:x_max], 'g--', label='Sqrt Eigenvalue')
    pylab.xlabel('Mode number')
    pylab.ylabel('log(sqrt(eigenvalue))')
    pylab.savefig(output+'.png')
    pylab.clf()

def ploteigenvalueslists(file_names, output, key = '21cm_only', x_max = 0):
    for file_name in file_names:
        file = h5py.File(file_name, 'r')
        vals = file['svd_vals'][key]
        vals = np.log10((vals/vals[0])**0.5)
        if x_max == 0:
            pylab.plot(range(256), vals, 'b',label='Sqrt Eigenvalue')
        else:
            pylab.plot(range(x_max), vals[0:x_max], 'b', label='Sqrt Eigenvalue')
    pylab.xlabel('Mode number')
    pylab.ylabel('log(sqrt(eigenvalue))')
    pylab.savefig(output+'.png')
    pylab.clf()

def ploteigenvectors(file_name, output, modes_key = 'svd_modes1',key = '21cm_only', first_n = 10, n_freqs = 256, f_limits = [700,900], gbt_cuts =  False, chris_cuts = False):
    file = h5py.File(file_name, 'r')
    vects = file[modes_key][key]
    #freqs = range(256)
    freqs=np.linspace(f_limits[0], f_limits[1], n_freqs)
    if gbt_cuts:
        #range 64-102, 109-129, 135-253
        freqs = np.linspace(900, 700, 256)
        #freqs = np.delete(freqs, range(0,63) + [103,104,105,106,107,108,130,131,132,133,134,254,255])
    for num in range(first_n):
        #freqs = range(len(vects[num]))
        the_vals = vects[num] + num + 1
        y_vals = the_vals
        if chris_cuts:
            y_vals = np.full(256,np.nan)
            y_vals[63:103] = the_vals[0:40]
            y_vals[109:130] = the_vals[40:61]
            y_vals[135:179] = the_vals[61:105]
            y_vals[218:246] = the_vals[105:133]
        elif gbt_cuts:
            y_vals = np.full(256,np.nan)
            y_vals[63:103] = the_vals[0:40]
            y_vals[109:130] = the_vals[40:61]
            y_vals[135:254] = the_vals[61:180]
        pylab.plot(freqs, y_vals)
        #plt.xticks([0,50,100,150,200,250])
        #plt.xticks([700,750,800,850,900])
        plt.yticks(range(first_n))
        #pylab.xlim([0,256])
        #pylab.xlim(f_limits)
        #pylab.ylim([-1,3])                        
    #pylab.xlabel('Frequency Channel')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Mode number')
    pylab.savefig(output + '.png')
    pylab.clf()

def ploteigenvectors_multiple(file_name1, file_name2, output, key = '21cm_only', first_n = 10, n_freqs = 256):
    file1 = h5py.File(file_name1, 'r')
    file2 = h5py.File(file_name2, 'r')
    vects1 = file1['svd_modes1'][key]
    vects2 = file2['svd_modes1'][key]
    #freqs = range(256)
    freqs=np.linspace(700, 900, n_freqs)
    for num in range(first_n):
        #freqs = range(len(vects[num]))
        pylab.plot(freqs, 3*vects1[num]+num)
        pylab.plot(freqs, 3*vects2[num]+num)
        #plt.xticks([0,50,100,150,200,250])
        plt.xticks([700,750,800,850,900])
        plt.yticks(range(1))
        #pylab.xlim([0,256])
        pylab.xlim([700,900])
        #pylab.ylim([-0.5,2.5])
    #pylab.xlabel('Frequency Channel')
    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Mode number')
    pylab.savefig(output + '.png')
    pylab.clf()

# Test below
def ploteigenvectors_test(file_name, output, key = '21cm_only',first_n = 5):
    fig, ax = plt.subplots()   
    file = h5py.File(file_name, 'r')
    vects = file['svd_modes1'][key]
#    freqs = range(256)
    for num in range(first_n):
        freqs = range(len(vects[num]))
        pylab.plot(freqs, vects[num]+num, label = 'Mode ' + str(num))

    pylab.xlabel('Frequency (MHz)')
    pylab.ylabel('Mode number')
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[0] = '700'    
    labels[1] = '850'
    ax.set_xticklabels(labels)
    ax.set_xticks([0, 180])
    ax.set_yticks([-1, 0, 1, 2, 3, 4, 5, 6, 7])   
    pylab.savefig(output + '.png')
    pylab.clf()

def plotmapslice(file_name, output, freq_slice = 100, max = False, avg = False, div = False, div_by=16. , map_by_hand = False, data=np.zeros((0,))):
    if not map_by_hand:
        map = al.load(file_name)
    else:
        map = data
    if not avg:
        a = map[freq_slice,:,:]
    if avg == True:
        a = np.mean(map,axis=0)
    if div:
        a /= div_by
    if max:
        im = pylab.imshow(a, cmap = 'hot', extent=[0, 30, 30, 0], vmin=-5, vmax=5)
    else:
        im = pylab.imshow(a, cmap = 'hot')
    plt.gca().invert_yaxis()
    pylab.colorbar(im, orientation='horizontal', label = 'Temperature (mK)')
    plt.xlabel('RA (degrees)', fontsize=14, color='black')
    plt.ylabel('Dec (degrees)', fontsize=14, color='black')
    pylab.savefig(output+'.png')
    pylab.clf()

def plot_map_mode(mode_path, mode_ind, output):
    modes = al.load(mode_path)
    #Normalize by freq
    modes /= modes.shape[0]**0.5
    info = modes.info
    middle = [info['ra_centre'], info['dec_centre']]
    delta = [info['ra_delta'], info['dec_delta']]
    pix_half = [modes.shape[1]/2,modes.shape[2]/2]
    extent = [middle[0] - pix_half[0]*delta[0], middle[0] + pix_half[0]*delta[0], middle[1] - pix_half[1]*delta[1], middle[1] + pix_half[1]*delta[1]]
    im = pylab.imshow(np.transpose(modes[mode_ind]), cmap = 'hot', origin = 'lower',  extent = extent)
    pylab.colorbar(im, orientation='horizontal', label = 'Temperature (K)')
    plt.xlabel('RA (degrees)', fontsize=14, color='black')
    plt.ylabel('Dec (degrees)', fontsize=14, color='black')
    pylab.savefig(output+'.png')
    pylab.clf()

def plot_map_modes(mode_path, mode_inds, output):
    modes = al.load(mode_path)
    #Normalize by freq
    modes /= modes.shape[0]**0.5
    info = modes.info
    middle = [info['ra_centre'], info['dec_centre']]
    delta = [info['ra_delta'], info['dec_delta']]
    pix_half = [modes.shape[1]/2,modes.shape[2]/2]
    extent = [middle[0] - pix_half[0]*delta[0], middle[0] + pix_half[0]*delta[0], middle[1] - pix_half[1]*delta[1], middle[1] + pix_half[1]*delta[1]]
    for mode_ind in mode_inds:
        pylab.subplot(math.ceil(np.float(len(mode_inds))/2), 2, mode_ind +1)
        im = pylab.imshow(np.transpose(modes[mode_ind]), cmap = 'hot', origin = 'lower',  extent = extent)
        cbar = pylab.colorbar(im, orientation='vertical')
        cbar.ax.tick_params(labelsize=8)
    #plt.xlabel('RA (degrees)', fontsize=14, color='black')
    #plt.ylabel('Dec (degrees)', fontsize=14, color='black')
    pylab.savefig(output+'.png')
    pylab.clf()

def plot1Dpowerspectrum(file_namep, file_namek, output):
    filep = np.load(file_namep)
    filek = np.load(file_namek)
    pylab.plot(np.log10(filek), np.log10(filep), label = '1-D Power Spectrum')
    pylab.xlabel('log(k)')
    pylab.ylabel('log(P(k))(K^2)')
    pylab.savefig(output + '.png')
    pylab.clf()

def plot1Dpowerspectrum_multiple(file_namep1, file_namek1, file_namep2, file_namek2, file_namep3, file_namek3, file_namep4, file_namek4, file_namep7, file_namek7, file_namep8, file_namek8, output):
    filep1 = np.load(file_namep1)
    filek1 = np.load(file_namek1)
    filep2 = np.load(file_namep2)
    filek2 = np.load(file_namek2)
    filep3 = np.load(file_namep3)
    filek3 = np.load(file_namek3) 
    filep4 = np.load(file_namep4)
    filek4 = np.load(file_namek4)
#    filep5 = np.load(file_namep5)
#    filek5 = np.load(file_namek5)
#    filep6 = np.load(file_namep6)
#    filek6 = np.load(file_namek6)
    filep7 = np.load(file_namep7)
    filek7 = np.load(file_namek7)
    filep8 = np.load(file_namep8)
    filek8 = np.load(file_namek8)
    pylab.plot(np.log10(filek1), np.log10(filep1), label = '0 modes removed')
    pylab.plot(np.log10(filek2), np.log10(filep2), label = '1 mode removed')
    pylab.plot(np.log10(filek3), np.log10(filep3), label = '2 modes removed')
    pylab.plot(np.log10(filek4), np.log10(filep4), label = '3 modes removed')
#    pylab.plot(np.log10(filek5), np.log10(filep5), label = '4 modes removed')
#    pylab.plot(np.log10(filek6), np.log10(filep6), label = '5 modes removed')
    pylab.plot(np.log10(filek7), np.log10(filep7), label = '10 modes removed')
    pylab.plot(np.log10(filek8), np.log10(filep8), label = '21 cm signal')
    pylab.xlabel('log(k)')
    pylab.ylabel(r'log(P(k)($mK^2$))')
    pylab.legend(loc='lower right',prop={'size':8})
    pylab.title(r'Power Spectra (21-cm Signal + Foregrounds + Noise)')
    pylab.savefig(output + '.png')
    pylab.clf()


def plot2Dpowerspectrum(file_namep, output, p_max=2):
    pwr = al.load(file_namep)
    pwr = np.ma.array(pwr, mask=pwr==0)
    if p_max == 0:
        im = pylab.imshow(np.flipud(np.transpose(pwr)), extent=[0.001, 5, 0.001, 5])
    else:
        im = pylab.imshow(np.flipud(np.transpose(pwr)), extent=[0.001, 5, 0.001, 5], vmax = p_max)
    pylab.colorbar(im, orientation='horizontal', label = r'$mK^2$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$k_{\perp}$', fontsize=14, color='black')
    plt.ylabel(r'$k_{\parallel}$', fontsize=14, color='black')
    pylab.savefig(output+'.png')
    pylab.clf()

def plot2Dpowerspectrumlognoneg(file_namep, output, normalize=False):
    pwr = al.load(file_namep)
    pwrp = np.ma.array(pwr, mask=pwr<=0)
#    mask = np.isfinite(pwrp)
#    if normalize:
#        print np.unravel_index(np.argmax(pwrp),pwrp.shape)
        #print(pwrp[:,9])
        #pwrp /= pwrp[:,9]
#        print(pwrp[33,:])
#        pwrp /= pwrp[33,:]
#        print(pwrp)
    #pwrp = np.ma.array(pwrp, mask != mask)
#    pwrp = np.ma.array(pwrp, mask=np.logical_or(pwrp<=0, np.isnan(pwrp)))
    fig = plt.figure()
    imp = pylab.imshow(np.flipud(np.transpose(pwrp)), extent=[0.001, 5, 0.001, 5], norm=LogNorm(vmin=np.min(pwrp), vmax=np.max(pwrp)))
    pylab.colorbar(imp, orientation='horizontal', label = r'$mK^2$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$k_{\perp}$', fontsize=14, color='black')
    plt.ylabel(r'$k_{\parallel}$', fontsize=14, color='black')
    pylab.savefig(output+'.png')
    pylab.clf()

def plot2Dpowerspectrumlog(file_namep, output):
    pwr = al.load(file_namep)
    pwrp = np.ma.array(pwr, mask=pwr<=0)
    pwrn = -1*np.ma.array(pwr, mask=pwr>=0)
    fig = plt.figure()
    a = fig.add_subplot(1,2,1)
    imp = pylab.imshow(np.flipud(np.transpose(pwrp)), extent=[0.001, 5, 0.001, 5], norm=LogNorm(vmin=np.min(pwrp), vmax=np.max(pwrp)))
    pylab.colorbar(imp, orientation='horizontal', label = r'$mK^2$')
#    pylab.colorbar(imp, orientation='horizontal', label = r'$mK^2$', ticks=[10**-4,10**0,10**3,10**7])
#    tick_locator = ticker.MaxNLocator(nbins=4)
#    cb.locator = tick_locator
#    cb.ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
#    cb.update_ticks()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$k_{\perp}$', fontsize=14, color='black')
    plt.ylabel(r'$k_{\parallel}$', fontsize=14, color='black')
    a = fig.add_subplot(1,2,2)
    imn = pylab.imshow(np.flipud(np.transpose(pwrn)), extent=[0.001, 5, 0.001, 5], norm=LogNorm(vmin=np.min(pwrn), vmax=np.max(pwrn)))
    pylab.colorbar(imn, orientation='horizontal', label = r'$mK^2$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$k_{\perp}$', fontsize=14, color='black')
    plt.ylabel(r'$k_{\parallel}$', fontsize=14, color='black')
    fig.tight_layout()
    pylab.savefig(output+'.png')
    pylab.clf()


def plotbeam(ang_range, freqs, output, y_max = 0.02, pix_size=1./15, achrom_airy=False):
    from map import beam 
    if not achrom_airy:
        airy_beam = beam.AiryBeam(100)
    else:
        airy_beam = beam.AchromaticAiryBeam(100)
    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                     0.281176247549, 0.270856788455, 0.26745856078,
                     0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                                 dtype=float)
    freq_data *= 1.0e6
    gaussian_beam = beam.GaussianBeam(beam_data, freq_data)
    angles = np.arange(ang_range[0], ang_range[1], pix_size) 
    for freq in freqs: 
        amp1 = airy_beam.beam_function(angles, 1000000 * np.array([freq]))
        amp2 = gaussian_beam.beam_function(angles, 1000000 * np.array([freq]))
        amp2 /= np.max(amp2)
        if not achrom_airy:
            pylab.plot(angles, amp1, color= 'r', label = 'Airy Beam ' + str(freq) + '(MHz)')
            pylab.plot(angles, amp2, color = 'b', label = 'Gaussian Beam ' + str(freq) + '(MHz)')
        else:
            pylab.plot(angles, amp1, label = 'Achromatic Airy Beam ' +  str(freq) + '(MHz)')
    #pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlabel('Angle (degrees)')
    pylab.ylabel('Amplitude')
    pylab.ylim((10**-7, 1))
    pylab.legend(loc='lower right',prop={'size':8})
    #pylab.title(r'Airy Beam Function')
    pylab.savefig(output + '.png')
    pylab.clf()


def makebeam(ang_range, freqs, y_max = 0.02, pix_size=1./15):
    from map import beam
    airy_beam = beam.AiryBeam(100)
    beam_data = sp.array([0.316148488246, 0.306805630985, 0.293729620792,
                     0.281176247549, 0.270856788455, 0.26745856078,
                     0.258910010848, 0.249188429031])
    freq_data = sp.array([695, 725, 755, 785, 815, 845, 875, 905],
                                 dtype=float)
    freq_data *= 1.0e6
    gaussian_beam = beam.GaussianBeam(beam_data, freq_data)
    angles = np.arange(ang_range[0], ang_range[1], pix_size)
    for freq in freqs:
        amp1 = airy_beam.beam_function(angles, 1000000 * np.array([freq]))
        amp2 = gaussian_beam.beam_function(angles, 1000000 * np.array([freq]))
    print(np.trapz(amp1))
    print(np.trapz(amp2))


def plotbp(file_name, output):
    bp = np.load(file_name)
    prylab.plot(bp)
    pylab.savefig(output + '.png')
    pylab.clf()


def plot_eigenvector_dot_products(file_name, output, pair_key = 'pair', evals_too = False, max_dot = False, eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='', common_var = False):
    file = h5py.File(file_name, 'r')
    vects1 = file['svd_modes1'][pair_key]
    vects2 = file['svd_modes2'][pair_key]
    dot_products = []
    for num in range(vects1.shape[1]):
        dot_products.append(np.dot(vects1[num],vects2[num]))
    dot_products = np.array(dot_products)
    fig, ax1 = pylab.subplots()
    if not (spatial_lr_dot or common_var):
        ax1.semilogx(1+np.arange(vects1.shape[1]), dot_products, 'b')
        ax1.set_xlabel('Mode number')
        ax1.set_ylabel('Dot Product of L and R Frequency Vector', color = 'b')
    if max_dot:
        ax1.set_ylim([-0.3,1.1])
    ax1.tick_params('y', colors='b')
    if common_var:
        vars = file[eval_key][pair_key]
        vars *= dot_products
        vars /= len(vars)
        tot_vars = np.array([np.sum(vars[i:]) for i in range(len(vars))])
        ax1.loglog(1+np.arange(vects1.shape[1]), tot_vars, 'b')
        ax1.set_xlabel('Mode number')
        ax1.set_ylabel('Variance with n-1 SVD modes removed (K^2), total is ' + str(tot_vars[0]), color = 'b')
    if evals_too:
        vals = file[eval_key][pair_key]
        vals = (vals/vals[0])**0.5
        ax2 = ax1.twinx()
        ax2.loglog(1+np.arange(vects1.shape[1]),vals, 'r')
        ax2.set_ylabel('Sqrt Mode Amplitude', color = 'r')
        ax2.tick_params('y', colors='r')
    if spatial_lr_dot:
        map1 = np.load(spatial_file1)
        map2 = np.load(spatial_file2)
        #normalize
        norm1 = np.transpose(np.transpose(map1)/np.transpose(np.sum(np.sum(map1**2,axis=2),axis=1)**0.5)) 
        norm2 = np.transpose(np.transpose(map2)/np.transpose(np.sum(np.sum(map2**2,axis=2),axis=1)**0.5))
        dots = np.sum(np.sum(norm1*norm2, axis=2),axis=1)
        ax1.semilogx(1+np.arange(dots.shape[0]), dots, 'g')
        ax1.set_xlabel('Mode number')
        ax1.set_ylabel('Dot Product of L and R Spatial Vector', color = 'g')
    pylab.savefig(output + '.png')       
    pylab.clf()

def left_deproject(mat, vec):
    #First dimension of mat must equal dim of vec.
    #Normalize vec.
    vec /= np.sqrt(np.dot(vec, vec))
    proj = np.transpose(np.transpose(np.ones(mat.shape)*np.tensordot(mat,vec, axes=[0,0]))*vec)
    ans = mat - proj
    return ans

if __name__ == '__main__':
    gbt_diff_path = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/map_diff_autos_1hr_80-68_azrm_short_svd_flagged5/orig_cuts_noweights/'
    diff_out = gbt_diff_path + 'plots/'
    diff_path_cuts = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/map_diff_autos_1hr_80-68_azrm_short_svd_flagged5/chris_cuts_noweights/'
    diff_cut_plots = diff_path_cuts + 'plots/'
    chris_cuts_dir = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/1hr_80-68_azrm_short_svd_flagged5/chris_cuts_131chan/'
    chris_cuts_plots = chris_cuts_dir + 'plots/'
    orig_cuts_dir = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/1hr_80-68_azrm_short_svd_flagged5/orig_cuts_180chan/'
    orig_cuts_plots = orig_cuts_dir + 'plots/'
    diff_cuts_deproj = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/map_diff_autos_1hr_80-68_azrm_short_svd_flagged5/chris_cuts_deproj_mapmodes/'
    diff_cuts_deproj_plots = diff_cuts_deproj + 'plots/'

    parkes_path = '/scratch2/p/pen/andersoc/second_parkes_pipe/cleaned_maps_bp_divide/hitconv_sync27_mbcal/I/ra165/'
    parkes_plots = '/scratch2/p/pen/andersoc/second_parkes_pipe/plots/SVD/first_round_ra165/'

    keys = ['A_minus_B', 'A_minus_C', 'A_minus_D', 'B_minus_C', 'B_minus_D', 'C_minus_D']
    key_names = ['A_with_B', 'A_with_C', 'A_with_D', 'B_with_C', 'B_with_D', 'C_with_D']
    keys_svd = ['beam1_with_beam2', 'beam1_with_beam3', 'beam1_with_beam4', 'beam2_with_beam3', 'beam2_with_beam4', 'beam3_with_beam4']
    
    parkes_keys = ['123_with_10to13','123_with_456','123_with_789',
 '456_with_10to13','456_with_789','789_with_10to13']

    gbt_deproj = '/scratch2/p/pen/andersoc/GBT_cleaning_results/1hr/1hr_80-68_azrm_short_svd_flagged5/deproj_30noise_chris_cuts/'
    gbt_deproj_plots = gbt_deproj + 'plots/'

    plot_add_on = '_freq_dot'
    plot_add_on2 = '_common_variance'

    for keys in zip(key_names[1:], keys_svd[1:]):
        key = keys[0]
        mode_key = keys[1]
        plot_eigenvector_dot_products(gbt_deproj + 'I/' + 'SVD.hd5', gbt_deproj_plots + key + plot_add_on, pair_key = key, evals_too = True, max_dot = False,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='')
        plot_eigenvector_dot_products(gbt_deproj + 'I/' + 'SVD.hd5', gbt_deproj_plots + key + plot_add_on2, pair_key = key, evals_too = True, max_dot = False,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='',common_var = True)
        #ploteigenvectors(gbt_deproj + 'I/' + 'SVD.hd5', gbt_deproj_plots + 'freq_vectors/' + 'map' + key[0] + '_' + key + '_6freq_modes', modes_key = 'svd_modes1',key = key, first_n = 6, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True, chris_cuts = True)
        #ploteigenvectors(gbt_deproj + 'I/' + 'SVD.hd5', gbt_deproj_plots + 'freq_vectors/' + 'map' + key[-1] + '_' + key + '_6freq_modes', modes_key = 'svd_modes2',key = key, first_n = 6, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True, chris_cuts = True)
        #for index in range(5):
            #plot_map_mode(gbt_deproj + 'I/' + 'sec_' + mode_key[0:5] + '_modes_clean_map_I_with_' + mode_key[-5:] + '_133modes.npy', mode_ind = index , output = gbt_deproj_plots + 'spatial_modes/' + key[0] + '_' + key + '_mode' + str(index))
            #plot_map_mode(chris_cuts_dir + 'I/' + 'sec_' + mode_key[0:5] + '_modes_clean_map_I_with_' + mode_key[-5:] + '_131modes.npy', mode_ind = index , output = chris_cuts_plots + 'spatial_modes/' + key[0] + '_' + key + '_mode' + str(index)) 
        mode_inds = range(6)
        plot_map_modes(diff_cuts_deproj+ 'I/' + 'sec_' + mode_key[0:5] + '_modes_clean_map_I_with_' + mode_key[-5:] + '_133modes.npy',  mode_inds = mode_inds , output = diff_cuts_deproj_plots + 'spatial_modes/' + key + '_6modes')
        plot_map_modes(chris_cuts_dir + 'I/' + 'sec_' + mode_key[0:5] + '_modes_clean_map_I_with_' + mode_key[-5:] + '_131modes.npy', mode_inds = mode_inds , output = chris_cuts_plots + 'spatial_modes/' + key[0] + '_' + key + '_6modes')
        plot_map_modes(chris_cuts_dir + 'I/' + 'sec_' + mode_key[-5:] + '_modes_clean_map_I_with_' + mode_key[0:5] + '_131modes.npy', mode_inds = mode_inds , output = chris_cuts_plots + 'spatial_modes/' + key[-1] + '_' + key + '_6modes')
        plot_map_modes(gbt_deproj + 'I/' + 'sec_' + mode_key[0:5] + '_modes_clean_map_I_with_' + mode_key[-5:] + '_133modes.npy', mode_inds = mode_inds ,output = gbt_deproj_plots + 'spatial_modes/' + key[0] + '_' + key +  '_6modes')
        plot_map_modes(gbt_deproj + 'I/' + 'sec_' + mode_key[-5:] + '_modes_clean_map_I_with_' + mode_key[0:5] + '_133modes.npy', mode_inds = mode_inds ,output = gbt_deproj_plots + 'spatial_modes/' + key[-1] + '_' + key +  '_6modes')

    keys = ['A_minus_B', 'A_minus_C', 'A_minus_D', 'B_minus_C', 'B_minus_D', 'C_minus_D']
    for key in keys[1:]:
        print key
        plot_eigenvector_dot_products(diff_cuts_deproj + 'I/' + 'SVD.hd5', diff_cuts_deproj_plots + key + plot_add_on, pair_key = key, evals_too = True, max_dot = True,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='')
        ploteigenvectors(diff_cuts_deproj + 'I/' + 'SVD.hd5', diff_cuts_deproj_plots + 'freq_vectors/' + key + '_10freq_modes', modes_key = 'svd_modes1',key = key, first_n = 10, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True, chris_cuts = True)
        #plot_eigenvector_dot_products(diff_path_cuts + 'I/' + 'SVD.hd5', parkes_plots + key + plot_add_on, pair_key = key, evals_too = True, max_dot = True,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='')
        #ploteigenvectors(diff_path_cuts + 'I/' + 'SVD.hd5', diff_cut_plots + 'freq_vectors/' + key + '_6freq_modes', modes_key = 'svd_modes1',key = key, first_n = 6, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True, chris_cuts = True)

    #for key in parkes_keys:
        #plot_eigenvector_dot_products(parkes_path + 'SVD.hd5', parkes_plots + key + plot_add_on, pair_key = key, evals_too = True, max_dot = False,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='')
        #plot_eigenvector_dot_products(parkes_path + 'SVD.hd5', parkes_plots + key, pair_key = key, evals_too = True, max_dot = True,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='')
        #ploteigenvectors(parkes_path + 'SVD.hd5', parkes_plots + 'freq_vectors/' + key + '_left_'  + '_6freq_modes', modes_key = 'svd_modes1',key = key, first_n = 6, n_freqs = 64, f_limits = [1283,1347], gbt_cuts =  False, chris_cuts = False)
        #ploteigenvectors(parkes_path + 'SVD.hd5', parkes_plots + 'freq_vectors/' + key + '_right_'  + '_6freq_modes', modes_key = 'svd_modes2',key = key, first_n = 6, n_freqs = 64, f_limits = [1283,1347], gbt_cuts =  False, chris_cuts = False)

    for key_pair in zip(keys_svd, key_names):
        plot_eigenvector_dot_products(chris_cuts_dir + 'I/' + 'SVD.hd5', chris_cuts_plots + key_pair[1] + plot_add_on2, pair_key = key_pair[0], evals_too = True, max_dot = False,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='', common_var = True)
        #plot_eigenvector_dot_products(chris_cuts_dir + 'I/' + 'SVD.hd5', chris_cuts_plots + key_pair[1] + plot_add_on, pair_key = key_pair[0], evals_too = True, max_dot = False,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='')
        #ploteigenvectors(chris_cuts_dir + 'I/' + 'SVD.hd5', chris_cuts_plots + 'freq_vectors/' + 'map' + key_pair[1][0] + '_' + key_pair[1] + '_6freq_modes', modes_key = 'svd_modes1',key = key_pair[0], first_n = 6, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True, chris_cuts = True)
        #ploteigenvectors(chris_cuts_dir + 'I/' + 'SVD.hd5', chris_cuts_plots + 'freq_vectors/' + 'map' + key_pair[1][-1] + '_' + key_pair[1] + '_6freq_modes', modes_key = 'svd_modes2',key = key_pair[0], first_n = 6, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True, chris_cuts = True)
        plot_eigenvector_dot_products(orig_cuts_dir + 'I/' + 'SVD.hd5', orig_cuts_plots + key_pair[1] + plot_add_on2, pair_key = key_pair[0], evals_too = True, max_dot = False,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='', common_var = True)
        #plot_eigenvector_dot_products(orig_cuts_dir + 'I/' + 'SVD.hd5', orig_cuts_plots + key_pair[1] + plot_add_on, pair_key = key_pair[0], evals_too = True, max_dot = False,eval_key = 'svd_vals', spatial_lr_dot = False, spatial_file1='', spatial_file2='')
        #ploteigenvectors(orig_cuts_dir + 'I/' + 'SVD.hd5', orig_cuts_plots + 'freq_vectors/' + 'map' + key_pair[1][0] + '_' + key_pair[1] + '_6freq_modes', modes_key = 'svd_modes1',key = key_pair[0], first_n = 6, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True)
        #ploteigenvectors(orig_cuts_dir + 'I/' + 'SVD.hd5', orig_cuts_plots + 'freq_vectors/' + 'map' + key_pair[1][-1] + '_' + key_pair[1] + '_6freq_modes', modes_key = 'svd_modes2',key = key_pair[0], first_n = 6, n_freqs = 256, f_limits = [700,900], gbt_cuts =  True)
