import numpy as np
import matplotlib.pyplot as plt

import core.algebra as al

from scipy.signal import fftconvolve 
import time 


file_name = "/mnt/raid-project/gmrt/eswitzer/GBT/cleaned_maps/GBT_15hr_map_autopaper/sec_C_cleaned_clean_map_I_with_B_20modes.npy"

map = al.load(file_name)
map = al.make_vect(map)

#convert nan to 0 in map
l,m,n = map.shape
for i in range(l):
	for j in range(m):
		for k in range(n):
			if np.isnan(map[i,j,k]):
				map[i,j,k]=0
del i,j,k,l,m,n


#axis values
freq = map.get_axis('freq')
ra = map.get_axis('ra')
dec = map.get_axis('dec')

#freq_centre = map.info['freq_centre'] # /!\ the real freq_centre is freq[127] (cf fast_gaussian_window)
#ra_centre = map.info['ra_centre']
#dec_centre = map.info['dec_centre']

freq_delta = map.info['freq_delta']
ra_delta = map.info['ra_delta']
dec_delta = map.info['dec_delta']

#gaussian window parameter 
sigma = 4*dec_delta

freq_centre = freq[127]
ra_centre = ra[38]
dec_centre = dec[21]

#different version of gaussian window
def gaussian_window(map):
	l,m,n = map.shape
	w = np.zeros((l,m,n))
	for k in range(l):
		for i in range(m):
			for j in range(n):
				r = (ra[i] - ra_centre)**2 + (dec[j] - dec_centre)**2 + (freq[k] - freq_centre)**2
#r = (i-n/2.)**2 + (j-l/2.)**2 + (k-m/2.)**2 # probably need to change this to get the exact value of coordinates x y z instead of indices
				w[k,i,j] = np.exp(-r/(2*sigma*sigma))
	return w

#best version
def fast_gaussian_window(map):
	l,m,n = map.shape
	w = np.zeros((l,m,n))
	for i in range(m):
		for j in range(n):
			r = (ra[i] - ra_centre)**2 + (dec[j] - dec_centre)**2
			w[127,i,j] = np.exp(-r/(2*sigma*sigma))
	return w

def full_gaussian_window(map):
	l,m,n = map.shape
	w = np.zeros((2*l-1,2*m-1,2*n-1))
	for k in range(l):
		for i in range(m):
			for j in range(n):
				r = (ra[i] - ra[0])**2 + (dec[j] - dec[0])**2 + (freq[k] - freq[0])**2
				val = np.exp(-r/(2*sigma*sigma))		
				w[l+k-1,m+i-1,n+j-1] = val
				w[l+k-1,m-i-1,n+j-1] = val
				w[l+k-1,m+i-1,n-j-1] = val
				w[l+k-1,m-i-1,n-j-1] = val
				w[l-k-1,m+i-1,n+j-1] = val
				w[l-k-1,m-i-1,n+j-1] = val
				w[l-k-1,m+i-1,n-j-1] = val
				w[l-k-1,m-i-1,n-j-1] = val
	return w

def fast_full_gaussian_window(map):
	l,m,n = map.shape
	w = np.zeros((2*l-1,2*m-1,2*n-1))
	for i in range(m):
		for j in range(n):
			r = (ra[i] - ra[0])**2 + (dec[j] - dec[0])**2
			val = np.exp(-r/(2*sigma*sigma))		
			w[255,m+i-1,n+j-1] = val
			w[255,m-i-1,n+j-1] = val
			w[255,m+i-1,n-j-1] = val
			w[255,m-i-1,n-j-1] = val
	return w

				
def partial_derivatives(w):
	l,m,n = w.shape
	w_x = np.zeros((l,m,n))
	w_y = np.zeros((l,m,n))		
	for i in range(m):
		x_i = ra[i] - ra_centre
		w_x[:,i,:] = -x_i/(sigma*sigma)*w[:,i,:]
	for j in range(n):
		x_j = dec[j] - dec_centre
		w_y[:,:,j] = -x_j/(sigma*sigma)*w[:,:,j]
	return w_x,w_y


def quadratic_estimators(map):
	t = time.time()
	w = fast_gaussian_window(map)
	(w_x,w_y) = partial_derivatives(w)
	del w
	delta_x = fftconvolve(map,w_x,'same')
	delta_y = fftconvolve(map,w_y,'same')
	del w_x, w_y
	t_xx = delta_x*delta_x
	t_yy = delta_y*delta_y
	t_xy = delta_x*delta_y
	t = time.time() - t
	print t
	return t_xx,t_yy,t_xy



	
	