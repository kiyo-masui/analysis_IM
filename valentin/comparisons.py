import numpy as np
import core.algebra as al
import time

from valentin.random_map import *
from scipy.signal import fftconvolve 

#---------------------------------------------------------------------------------------

#Different versions of convolutions

#convolution so that x=0 is at index (l/2,m/2,n/2)   ('same' considers x=0 at (l-1)/2, etc.)
def convolve(map1, map2):
	(l,m,n) = map1.shape
	full_map = fftconvolve(map1, map2, 'full')
	same_map = np.zeros((l,m,n))

	for k in range(l):
		for i in range(m):
			for j in range(n):
				same_map[k,i,j] = full_map[(l/2 + k), (m/2 + i), (n/2 + j)]
	return same_map


#convolution, the maps are assumed to represent one period of a periodic universe (same assumption as for Fourier transform)
def periodic_convolve(map1,map2):
	(l,m,n) = map1.shape
	full_map = fftconvolve(map1, map2, 'full')
	periodic_map = np.zeros((l,m,n))

	full_map[:,:,:(n-1)] += full_map[:,:,n:]
	full_map[:,:,n:] = full_map[:,:,:(n-1)]
	full_map[:,:(m-1),:] += full_map[:,m:,:]
	full_map[:,m:,:] = full_map[:,:(m-1),:]
	full_map[:(l-1),:,:] += full_map[l:,:,:]
	full_map[l:,:,:] = full_map[:(l-1),:,:]

	for k in range(l):
		for i in range(m):
			for j in range(n):
				periodic_map[k,i,j] = full_map[(l/2 + k), (m/2 + i), (n/2 + j)]
	return periodic_map

#convolution via Fourier space (same assumption as above)
def own_fftconvolve(map1,map2):
	f_map1 = np.fft.ifftshift(map1)
	f_map2 = np.fft.ifftshift(map2)
	f_map1 = np.fft.fftn(f_map1)
	f_map2 = np.fft.fftn(f_map2)
	f_conv = f_map1 * f_map2
	f_conv = np.fft.ifftn(f_conv)
	f_conv = np.fft.fftshift(f_conv)
	return f_conv.real

def conv_compare():
	map1 = map_init(100,100,100)
	map2 = map_init(100,100,100)
	
	t = time.time()
	fftconvolve(map1,map2,'full')
	t = time.time() - t
	print "fftconvolve :",t
	
	t = time.time()
	own_fftconvolve(map1,map2)
	t = time.time() - t
	print "own fftconvolve :",t	
	
	t = time.time()
	convolve(map1,map2)
	t = time.time() - t
	print "convolve :",t
	
	t = time.time()
	periodic_convolve(map1,map2)
	t = time.time() - t
	print "periodic convolve :",t
	
#after comparison, own_fftconvolve is really faster (more than 50x)

#---------------------------------------------------------------------------------------

#Attempt to create a fast version of kernel calculation, avoiding triple for loop
def fast_kernel(k_axes):
	kz = k_axes[0]
	kx = k_axes[1]
	ky = k_axes[2]
	k_xx = (kx*kx)[:,None] * np.ones((1,len(ky)))
	k_yy = np.ones((len(kx),1)) * (ky*ky)[None,:]
	k_xy = kx[:,None] * ky[None,:]
	rad_2d = k_xx + k_yy #calculate kx^2 + ky^2
	k_xx = k_xx/rad_2d - 0.5
	k_yy = k_yy/rad_2d - 0.5
	k_xy = k_xy/rad_2d
#	try:
	k_z_3d = (kz*kz)[:,None,None] * (1/rad_2d)[None,:,:] #a tester!!
#	except Exception:
#		print "Ok"
	k_xx = (1 + k_z_3d)*k_xx[None,:,:]
	k_yy = (1 + k_z_3d)*k_yy[None,:,:]
	k_xy = (1 + k_z_3d)*k_xy[None,:,:]
	k_xx[:,len(kx)/2,len(ky)/2] = 0. #Temporary solution, ok only if x_centre,y_centre = 0
	k_yy[:,len(kx)/2,len(ky)/2] = 0.
	k_xy[:,len(kx)/2,len(ky)/2] = 0.
		#Shift so that the fonction ifftn work
	k_xx = np.fft.ifftshift(k_xx)
	k_yy = np.fft.ifftshift(k_yy)
	k_xy = np.fft.ifftshift(k_xy)
	return k_xx,k_yy,k_xy
	
#Old version
def map_kernel(k_axes):
	t = time.time()
	l,m,n = self.info['shape']
	ker_xx = np.zeros((l,m,n))
	ker_yy = np.zeros((l,m,n))
	ker_xy = np.zeros((l,m,n))
	
	for k in range(l):
		for i in range(m):
			for j in range(n):
				k_z = k_axes[0][k]
				k_x = k_axes[1][i]
				k_y = k_axes[2][j]
				cst2 = k_x*k_x + k_y*k_y + 0.0
				if (cst2 != 0):
					cst = (cst2 + k_z*k_z)/cst2
					ker_xx[k,i,j] = cst*(k_x*k_x/cst2 - 0.5)
					ker_yy[k,i,j] = cst*(k_y*k_y/cst2 - 0.5)
					ker_xy[k,i,j] = cst*(k_x*k_y/cst2)
					
		#Shift so that the fonction ifftn work
	ker_xx = np.fft.ifftshift(ker_xx)
	ker_yy = np.fft.ifftshift(ker_yy)
	ker_xy = np.fft.ifftshift(ker_xy)
	
	t = time.time() - t
	print t
	return ker_xx,ker_yy,ker_xy
	
def kernel_compare():
	k_axis = 2*np.pi/322.*(np.arange(128) - 64)
	t = time.time()
	map_kernel([k_axis,k_axis,k_axis])
	t = time.time() - t
	print '3 for loops version : %.3f sec' %
	t = time.time()
	fast_kernel([k_axis,k_axis,k_axis])
	t = time.time() - t
	print 'fast version : %.3f sec' %t

#after comparison, the fast version is really faster (more than 100x)

#---------------------------------------------------------------------------------------

#Two ideas for implementation of gaussian window calculation
def gaussian_window(map, sigma):
	z = map.get_axis('z')
	x = map.get_axis('x')
	y = map.get_axis('y')
	mesh = np.indices(map.shape)
	sph = np.zeros(mesh.shape)
	sph[0] = z[mesh[0]] - map.info['z_centre']
	sph[1] = x[mesh[1]] - map.info['x_centre']
	sph[2] = y[mesh[2]] - map.info['y_centre']
	rad = sph[0]**2 + sph[1]**2 + sph[2]**2
	return (sigma*np.sqrt(2*np.pi))**(-3) * np.exp(-rad/(2.*sigma*sigma))	
	
def alt_gaussian(map,sigma):
	l,m,n = map.shape
	z = map.get_axis('z')
	x = map.get_axis('x')
	y = map.get_axis('y')
	rad_3d = (z*z)[:,None,None] * np.ones((1,m,n))
	rad_3d += (x*x)[None,:,None] * np.ones((l,1,n))
	rad_3d += (y*y)[None,None,:] * np.ones((l,m,1))
	return (sigma*np.sqrt(2*np.pi))**(-3) * np.exp(-rad_3d/(2.*sigma*sigma))

def gaus_compare():
	map = sim_random_map((512,512,512),(100.,100.,100.))
	t = time.time()
	gaussian_window(map, 8.)
	t = time.time() - t
	print 'v1 : %.3f sec' %t
	t = time.time()
	alt_gaussian(map, 8.)
	t = time.time() - t
	print 'v2 : %.3f sec' %t
	
#after comparison, alt_version is faster (about 2x)

#---------------------------------------------------------------------------------------

#Attempt to create a fast version of noisy modes filter, avoiding triple for loop
def fast_kz_filter(kappa, k_axes):
	l,m,n = kappa.shape
	kz = np.fft.ifftshift(k_axes[0])
	kx = np.fft.ifftshift(k_axes[1])
	ky = np.fft.ifftshift(k_axes[2])
	rad_2d = (kx*kx)[None,:,None] * np.ones((l,1,n))
	rad_2d += (ky*ky)[None,None,:] * np.ones((l,m,1))
	filt = rad_2d - (kz*kz)[:,None,None]*np.ones((1,m,n))
	filt = filt > 0
	f_kappa = np.fft.fftn(kappa)
	f_kappa *= filt
	return np.fft.ifftn(f_kappa).real

#Old version	
def kz_filter(kappa, k_axes):
	kz = np.fft.ifftshift(k_axes[0])
	kx = np.fft.ifftshift(k_axes[1])
	ky = np.fft.ifftshift(k_axes[2])
	l,m,n = kappa.shape
	f_kappa = np.fft.fftn(kappa)
	for k in range(l):
		for i in range(m):
			for j in range(n):
				if (kx[i]*kx[i] + ky[j]*ky[j]) <= kz[k]*kz[k]:
					f_kappa[k,i,j]=0
	return np.fft.ifftn(f_kappa).real

def kz_compare():
	kappa = np.random.normal(size = (128,128,128), scale = 1.)
	k_axis = 2*np.pi/322.*(np.arange(128) - 64)
	t = time.time()
	kz_filter(kappa, [k_axis,k_axis,k_axis])
	t = time.time() - t
	print 'old version : %.3f sec' %t
	t = time.time()
	fast_kz_filter(kappa, [k_axis,k_axis,k_axis])
	t = time.time() - t
	print 'fast version : %.3f sec' %t
	
#after comparison, the fast version is faster
#but could be even faster if kappa was already in Fourier space

#---------------------------------------------------------------------------------------

##Attempt to create a fast version of quadratic estimators, reducing the number of fft
#def quad_est_compare(map_s):
#	try :
#		t = time.time()
#		map_s.fast_quad_est()
#		t2 = time.time() - t
#		
#		t = time.time()
#		map_s.gaussian_window(map_s.info['sigma_recon'])
#		map_s.quadratic_estimators(map_s.info['sigma_recon'],log_smooth=True)
#		t1 = time.time() - t
#		
#		print '---------------------------'
#		print 'old version : %.3f sec' %t1
#		print 'fast version : %.3f sec' %t2
#	except NameError:
#		print "mids must be initialized"
		


