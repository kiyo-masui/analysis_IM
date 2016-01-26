import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import itertools
from matplotlib.backends.backend_pdf import PdfPages
from scipy.ndimage.filters import convolve
import copy
from mpi4py import MPI

comm = MPI.COMM_WORLD

data_dir = '/scratch2/p/pen/andersoc/cluster_hi_agn_sim_navarro/guest/'

halo = 'Halos_'
hi = 'Region_'
suff = '.txt'

c = 300000
# c in km/s
H = 100
# H in km s-1 Mpc-1
f0 = 1315
beam_rad = 1.4*(14.1/60.0)*(np.pi/180)*(1394.5/1347.5)
z0 = (1420.0/1315.0) - 1
h = 6.626*10**-34
# kg m^2 s-1
k = 1.38*10**-23
# kg m^2 s-2 K-1
f_rest = 1420.4*10**6
# Hz
A10 = 2.85*10**-15
# Hz
Ms = 1.99*10**30
Mhi = 1.67*10**-27
mpc = 3.086*10**19
# km

def data_array(type):
    for num in range(29):
        num += 1
        a = np.loadtxt(data_dir + type + str(num) + suff)
        if num == 1:
            b = a
        else:
            b = np.concatenate((b,a), axis=0)
        print str(num) + ' , now ' + str(b.shape[0]) + ' halos.'
    return b

def data_region(type, region):
    a = np.loadtxt(data_dir + type + str(region) + suff)
    return a

def dist_sq(p1, p2):
    #p2 can be an array of coordinates, axis 0 is the particle, rest is x,y,z
    delta = p1 - p2
    return np.sum(np.square(delta), axis = p2.ndim -1)

def halo_masked(parts, center, radius):
    dist = dist_sq(center, parts[:,0:3])
    mask = dist<=radius**2
    return parts[mask]

def p_hi_mass(parts, p_cen, p_size):
    #p_size must be numpy array
    # returns hi mass in single spatial pixel in Msun/h
    p_min = p_cen - p_size/2.
    p_max = p_cen + p_size/2.
    test_min = sp.greater(parts[:,0:3],p_min)
    test_max = sp.less(parts[:,0:3],p_max)
    mask = np.logical_and(np.all(test_min, 1),np.all(test_max, 1))
    mass = np.sum(parts[:,6][mask])
    return mass

def dop_shift(parts, rad_ax):
    ax = {'x':3,'y':4, 'z':5}
    pos_ax = {'x':0,'y':1, 'z':2}
    #parts is the list of HI particles
    #rad_ax: either x,y, or z.  Which axis to take as radial.
    gamma = (1-(parts[:,3]**2 + parts[:,4]**2 + parts[:,5]**2)/c**2)**-0.5
    beta_ax = parts[:,ax[rad_ax]]/c
    z = gamma*(1 + beta_ax)
    #parts[:,pos_ax[rad_ax]] += (c*dz)/H
    #parts[:,pos_ax[rad_ax]] = f0/z
    return f0/z

def freq_grid(parts, f_ax, center):
    ax = {'x':0,'y':1, 'z':2}
    #Will regrid f_ax to freq coordinates, assuming 1315.5 MHz at Halo center.
    #z = H*d/c
    #f = f0/(z+1), with f0 at 1315.5 MHz
    parts[:,ax[f_ax]] -= center
    parts[:ax[f_ax]] = 1315.5/((H*parts[:ax[f_ax]])/c + 1)

def hi_hist(parts, p_cen, p_size, rad):
    weights = parts[:,6]
    pos = parts[:,0:3]
    range = [[p_cen[0] - rad[0], p_cen[0] + rad[0]], [p_cen[1] - rad[1], p_cen[1] + rad[1]], [p_cen[2] - rad[2], p_cen[2] + rad[2]] ]
    bins = 2*rad/p_size
    hi_mass = np.histogramdd(pos, bins = bins, range = range, weights = weights)
    return hi_mass

def hi_freq_hist(part, center, f_ax, r200, p_size, f_p_size, r_factor, f_rad=0, r200_overwrite=0, brightness='True'):
    # parts is HI particles data, cent is halo center in physical coordinates,
    # f_ax is the axis that will become freq, r200 is the halo r200,
    # p_size is the gridding size in real space, f_p_size is same 
    # thing in freq space, in MHz, r_factor is factor beyond r200 to extend
    # the hitmap.  It also applies to the freq range.
    # If brightness is true, will compute brightness temp difference due to
    # HI.  Assumes optically thin, Ts>>Tcmb limit.
    
    cent = copy.deepcopy(center)
    parts = copy.deepcopy(part)
    pos_ax = {'x':0,'y':1, 'z':2}
    #First, find HI particles within r200 of center.
    halo_hi = halo_masked(parts, cent, r200)
    #print halo_hi.shape
    #Now, find the freq range of these particles
    freqs = dop_shift(halo_hi, f_ax)
    #print freqs
    if halo_hi.shape[0] == 0:
        f_range = [f0,f0]
    else:
        f_range = [np.min(freqs), np.max(freqs)]
    #print f_range
    if halo_hi.shape[0] == 0:
        f_cen = 0
    else:
        f_cen = np.sum(np.multiply(freqs, halo_hi[:,6]))/np.sum(halo_hi[:,6])
    if f_rad ==0:
        f_rad = np.max(np.absolute(f_cen - f_range))
    parts[:,pos_ax[f_ax]] = dop_shift(parts, f_ax)

    p_size = [p_size, p_size, p_size]
    p_size[pos_ax[f_ax]] = f_p_size

    if r200_overwrite !=0:
        r200 = r200_overwrite
    rad = np.array([r200, r200, r200])
    rad[pos_ax[f_ax]] = f_rad
    rad *= r_factor
    #print rad

    f_cent = cent
    f_cent[pos_ax[f_ax]] = f_cen
    #print f_cent
    #hist = hi_hist(parts, f_cent, p_size, rad)

    if brightness:
        coeff = (3*(c**2)/(32*np.pi*p_size[0]**2))*(A10*h)*(1/(f_rest*k))*(1/(f_p_size*10.0**6))*(Ms/Mhi)*(1/mpc**2)
        parts[:,6] *= coeff
    #print parts.shape
    #print f_cent
    #print p_size,
    #print rad
    hist = hi_hist(parts, f_cent, p_size, rad)

    return hist

def beam_kernel(fwhm, p_size, sigs):
    # Sigs is the number of sd's to extend the beam profile to.
    # 
    sig = fwhm / (2. * sp.sqrt(2. * sp.log(2.)))
    size = int(np.ceil(sigs*sig/p_size))
    size *= 2
    size += 1
    dr_sq = np.fromfunction(lambda i, j: ((i-(size-1)/2)*(p_size/sig))**2 + ((j-(size-1)/2)*(p_size/sig))**2, (size,size), dtype=float)
    profile = sp.exp(-dr_sq / (2. * sig ** 2.))
    profile *= 1. / (2. * sp.pi * sig ** 2.)
    profile /= np.sum(profile)
    return profile

def conv(input, beam, pad):
    #Really just a wrapper of scipy function.
    # pad should be zero, for zero padding at borders.
    return convolve(input, beam, mode='constant', cval=pad)

def hist_conv(hist, fwhm, p_size, sigs, pad):
    #Assumes that the last index is the frequency index.
    beam = beam_kernel(fwhm, p_size, sigs)
    for f in range(hist[0].shape[2]):
        #print hist[0][:,:,f].shape
        #print beam.shape
        hist[0][:,:,f] = conv(hist[0][:,:,f], beam, pad)

def hi_pdf(list, name, same_scale='False'):
    pdf = PdfPages(name + '.pdf')
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 9}

    matplotlib.rc('font', **font)
    hist = list[0]
    edges = list[1]
    for d in xrange(hist.shape[2]):
        if not same_scale:
            plt.imshow(hist[:,:,d], extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]])
        else:
            plt.imshow(hist[:,:,d], extent = [edges[0][0], edges[0][-1], edges[1][0], edges[1][-1]],vmin=0,vmax=np.max(hist))
        plt.colorbar()
        pdf.savefig()
        plt.close()
    pdf.close()
    

def scatter_halos(halos, filename):
    fig, ax = plt.subplots()
    ax.set_xlim([np.min(halos[:,3]), np.max(halos[:,3])])
    ax.set_ylim([np.min(halos[:,4]), np.max(halos[:,4])])
    ax.scatter(halos[:,3],halos[:,4])
    fig.savefig(filename)
    plt.clf()

def hi_halo_stack(halo_list, hi_list, f_ax, f_rad, pix, pix_f, r_fac, rad, fwhm, sigs):
    for index in xrange(halo_list.shape[0]):
        center = halo_list[index,0:3]
        #print center
        r200 = halo_list[index,3]
        hist = hi_freq_hist(hi_list, center, f_ax, r200, pix, pix_f, r_fac, f_rad=f_rad, r200_overwrite = rad)
        hist_conv(hist, fwhm, pix, sigs, 0)
        if index ==0:
            tot = hist
        else:
            ref = hist[0]
            ref_tot = tot[0]
            ref_tot += ref
    ref = tot[0]
    #ref /= halo_list.shape[0]
    return tot

def hi_full_stack(f_ax, f_rad, pix, pix_f, r_fac, rad, fwhm, sigs, m_min=0, m_max=10**16, r_min =0, r_max = 10):
    rank = comm.Get_rank()
    size = comm.Get_size()
    per = -(-29//size)
    if rank == size -1:
        region_list = range(rank*per + 1, 30)
    else:
        region_list = range(rank*per + 1, (rank+1)*per+1)
    print size
    print rank
    print region_list
    halo_count = 0
    for index in region_list:
        print str(index)
        halo_list = data_region(halo, index)
        mass_mask = halo_list[:,4] <= m_max 
        halo_list = halo_list[mass_mask]
        mass_mask = halo_list[:,4] > m_min
        halo_list = halo_list[mass_mask]
        r_mask = halo_list[:,3] <= r_max
        halo_list = halo_list[r_mask]  
        r_mask = halo_list[:,3] > r_min
        halo_list = halo_list[r_mask]

        halo_count += halo_list.shape[0]
        hi_list = data_region(hi, index)
        if halo_list.shape[0] > 0: 
            hist = hi_halo_stack(halo_list, hi_list, f_ax, f_rad, pix, pix_f, r_fac, rad, fwhm, sigs)
            try:
                tot
            except NameError:
                tot = hist
            else:
                ref = hist[0]
                ref_tot = tot[0]
                ref_tot += ref
    ref = tot[0]
    #ref /= per
    comm.Barrier()
    temp = np.zeros(ref.shape)
    halo_node = 0
    for node in range(1,size):
        if rank == node:
            comm.Send(ref, dest=0, tag=13)
            comm.Send(halo_count, dest=0, tag=14)
            print 'Process ' + str(node) + ' sent hist to node zero.'
        if rank == 0:
            comm.Recv(temp, source = node, tag=13)
            comm.Recv(halo_node, source = node, tag=13)
            halo_count += halo_node
            ref += temp
            print 'Node zero received hist from node ' + str(node) + '.'
    comm.Barrier()
    ref /= halo_count
    print 'Node ' + str(rank) + ' has ' + str(halo_count) + ' Halos.'
    return tot

def halo_r200_hist(halos, filename):
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 17}
    n, bins, patches = plt.hist(halos[:,3], 50)
    plt.axis([np.min(halos[:,3]),np.max(halos[:,3]),0,np.max(n)])
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    matplotlib.rc('font', **font)
    plt.xlabel('R200 in Mpc/h')
    plt.ylabel('Halo count, total 449')
    plt.savefig(filename)
    return [n,bins,patches]

def halo_mass_hist(halos, filename):
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 17}
    n, bins, patches = plt.hist(halos[:,4], 50)
    plt.axis([np.min(halos[:,3]),np.max(halos[:,4]),0,np.max(n)])
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    matplotlib.rc('font', **font)
    plt.xlabel('Mass in solar M/h')
    plt.ylabel('Halo count, total 449')
    plt.savefig(filename)
    return [n,bins,patches]

if __name__ == '__main__':
    #full_stack = hi_full_stack('z', 10, .1, 1, 1.2, 2.3, beam_rad*c*z0/H, 3, m_min = 10**12.9192, m_max = 10**16)
    stack = hi_full_stack('z', 10, .1, 1, 1.2, 2.3, beam_rad*c*z0/H, 3, r_min = 1.0)
    print 'Next'
    if comm.Get_rank() == 0:
        print 'Making pdf'
        hi_pdf(stack, '/gss01/scratch2/p/pen/andersoc/workcopy/parkes_analysis_IM/hi_r200above1Mpc_stack_rsd', same_scale='True') 
