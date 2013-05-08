import numpy as np
import core.algebra as al

import matplotlib.pyplot as plt

from scipy.signal import fftconvolve 
import time


#generate a random map from a Gaussian distribution
def map_init(l_map = 64,m_map = 128,n_map = 128):
	sigma_map = 1. 
	return np.random.normal(size=(l_map,m_map,n_map), scale=sigma_map)


def gaussian_window(map, sigma=1.):
	l,m,n = map.shape
	w = np.zeros((l,m,n))
	z_0 = l/2
	x_0 = m/2
	y_0 = n/2
	for k in range(l):
		z_k = k - z_0
		for i in range(m):
			x_i = i - x_0
			for j in range(n):
				y_j = j - y_0
				r = x_i**2 + y_j**2 + z_k**2
				w[k,i,j] = np.exp(-r/(2.*sigma*sigma))
	w = (sigma*np.sqrt(2*np.pi))**(-3) * w		
	return w
	

def partial_derivatives(w, sigma=1.):
	l,m,n = w.shape
	z_0 = l/2
	x_0 = m/2
	y_0 = n/2
	w_x = np.zeros((l,m,n))
	w_y = np.zeros((l,m,n))		
	for i in range(m):
		x_i = i - x_0
		w_x[:,i,:] = -x_i/(sigma*sigma)*w[:,i,:]
	for j in range(n):
		y_j = j - y_0
		w_y[:,:,j] = -y_j/(sigma*sigma)*w[:,:,j]
	return w_x,w_y

	
def quadratic_estimators(map, sigma=1.):
	t = time.time()
	w = gaussian_window(map, sigma)
	(w_x,w_y) = partial_derivatives(w, sigma)
	del w
	delta_x = convolve(map,w_x)
	delta_y = convolve(map,w_y)
	del w_x, w_y
	t_xx = delta_x*delta_x
	t_yy = delta_y*delta_y
	t_xy = delta_x*delta_y
	t = time.time() - t
	print t
	return t_xx,t_yy,t_xy
	
	
def map_kernel(map):
	t = time.time()
	l,m,n = map.shape
	ker_xx = np.zeros((l,m,n))
	ker_yy = np.zeros((l,m,n))
	ker_xy = np.zeros((l,m,n))
	
	for k in range(l):
		for i in range(m):
			for j in range(n):
				x = (i - m/2) #/((m-1)*d_x)
				y = (j - n/2) #/((n-1)*d_y)
				z = (k - l/2) #/((l-1)*d_z)
				cst2 = x*x + y*y + 0.0
				if (cst2 != 0):
					cst = (cst2 + z*z)/cst2
					ker_xx[k,i,j] = cst*(x*x/cst2 - 0.5)
					ker_yy[k,i,j] = cst*(y*y/cst2 - 0.5)
					ker_xy[k,i,j] = cst*(x*y/cst2)
					
	#Shift so that the fonction ifftn work
	ker_xx = np.fft.ifftshift(ker_xx)
	ker_yy = np.fft.ifftshift(ker_yy)
	ker_xy = np.fft.ifftshift(ker_xy)
	
	t = time.time() - t
	print t
	return ker_xx,ker_yy,ker_xy
	
	
	
def density_field(t_set, k_set):
	t = time.time()
	(t_xx,t_yy,t_xy) = t_set
		#Shift so that x=0 is the first element of the matrix
	t_xx = np.fft.ifftshift(t_xx)
	t_yy = np.fft.ifftshift(t_yy)
	t_xy = np.fft.ifftshift(t_xy)
		#Apply the fft
	t_xx = np.fft.fftn(t_xx)
	t_yy = np.fft.fftn(t_yy)
	t_xy = np.fft.fftn(t_xy)
		#Get kappa in k space
	kappa_xx = t_xx * k_set[0]
	kappa_yy = t_yy * k_set[1]
	kappa_xy = t_xy * k_set[2]
		#Get kappa in real space
	kappa_xx = np.fft.ifftn(kappa_xx)
	kappa_yy = np.fft.ifftn(kappa_yy)
	kappa_xy = np.fft.ifftn(kappa_xy)
		#Keep only the real part (imaginary part should be negligible)
	kappa_xx = kappa_xx.real
	kappa_yy = kappa_yy.real
	kappa_xy = kappa_xy.real
		#Shift so that x=0 is the centre of the matrix
	kappa_xx = np.fft.fftshift(kappa_xx)
	kappa_yy = np.fft.fftshift(kappa_yy)
	kappa_xy = np.fft.fftshift(kappa_xy)
		#Merge contributions to kappa
	kappa = kappa_xx + kappa_yy + 2*kappa_xy
	t = time.time() - t
	print t
	return kappa


def full_process(l_map = 64, m_map = 128, n_map = 128):
	map = map_init(l_map,m_map,n_map)
	sigma_w = 4.
		
	t = time.time()
	print "Computing quadratic estimators in ..."
	t_set = quadratic_estimators(map, sigma_w)
	print "Computing kernel in ..."
	k_set = map_kernel(map)
	print "Compting convolution in ..."
	kappa = density_field(t_set,k_set)
	t = time.time() - t
	print "Total time : ",t
	return map,kappa



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
	
#	for j in range(n-1):
#		full_map[:,:,j] += full_map[:,:,n+j]
#		full_map[:,:,n+j] = full_map[:,:,j]
#	for i in range(m-1):
#		full_map[:,i,:] += full_map[:,m+i,:]
#		full_map[:,m+i,:] = full_map[:,i,:]
#	for k in range(l-1):	
#		full_map[k,:,:] += full_map[l+k,:,:]
#		full_map[l+k,:,:] = full_map[k,:,:]
		
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


	
k_xx = k_set[0]
t_xx = t_set[0]

if_t_xx = np.fft.ifftshift(t_xx)
if_t_xx = np.fft.fftn(if_t_xx)
if_kap_xx = if_t_xx * k_xx
if_kap_xx = np.fft.ifftn(if_kap_xx)
if_kap_xx = if_kap_xx.real
if_kap_xx = np.fft.fftshift(if_kap_xx)

f_k_xx = np.fft.ifftn(k_xx)
f_k_xx = f_k_xx.real
f_k_xx = np.fft.fftshift(f_k_xx)
f_kap_xx = periodic_convolve(t_xx, f_k_xx)

f_k_yy = np.fft.ifftn(k_set[1])
f_k_yy = f_k_yy.real
f_k_yy = np.fft.fftshift(f_k_yy)
f_kap_yy = periodic_convolve(t_set[1], f_k_yy)

f_k_xy = np.fft.ifftn(k_set[2])
f_k_xy = f_k_xy.real
f_k_xy = np.fft.fftshift(f_k_xy)
f_kap_xy = periodic_convolve(t_set[2], f_k_xy)
kappa2 = f_kap_xx + f_kap_yy + 2*f_kap_xy

