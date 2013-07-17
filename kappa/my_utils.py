import numpy as np
import core.algebra as al
import matplotlib.pyplot as plt
from utils import batch_handler
import time

#----------------------------------------
#Map initialization
#----------------------------------------

#generate a random map from a Gaussian distribution
#map_shape is the number of pixels per dimension, map_size is the physical size of the map (in Mpc/h) 
def sim_random_map(map_shape, map_size, sigma_map=1.):
	map = np.random.normal(size=map_shape, scale=sigma_map)
	map = al.make_vect(map,('z','x','y'))
	map.set_axis_info('z', 0., map_size[0]/map.shape[0])
	map.set_axis_info('x', 0., map_size[1]/map.shape[1])
	map.set_axis_info('y', 0., map_size[2]/map.shape[2])
	return map

#get large map
def large_init(k=4):
	path = '/mnt/raid-project/gmrt/goblot/simulations/'
	map = np.load(path+'recon67_0.054den_512px.npy')
	map = low_resolution(map, k)
	map = al.make_vect(map,('z','x','y'))
	map.set_axis_info('z',0.,k*322./512)
	map.set_axis_info('x',0.,k*322./512)
	map.set_axis_info('y',0.,k*322./512)
	return map

#----------------------------------------
#Useful stuff
#----------------------------------------

#convolution using fft
@batch_handler.log_timing
def own_fftconvolve(map1,map2):
	f_map1 = np.fft.ifftshift(map1)
	f_map2 = np.fft.ifftshift(map2)
	del map1,map2
	t1 = time.time()
	f_map1 = np.fft.fftn(f_map1)
	t2 = time.time()
	print '[TIMING] fft:',t2-t1
	f_map2 = np.fft.fftn(f_map2)
	t1 = time.time()
	print '[TIMING] fft:',t1-t2
	f_conv = f_map1 * f_map2
	t2 = time.time()
	print '[TIMING] multiplication:',t2-t1
	del f_map1,f_map2
	f_conv = np.fft.ifftn(f_conv)
	t1 = time.time()
	print '[TIMING] ifft:',t1-t2
	f_conv = np.fft.fftshift(f_conv)
	return f_conv.real

#animation to show every slices z
def show_seq(mat,step=0.2, ref=0, extent=None):
	n = mat.shape[0]
	fig, ax = plt.subplots()

	im = ax.imshow(mat[ref,:,:].T,extent=extent)
	fig.show()
	for i in range(n):
		im.set_data(mat[i,:,:].T)
		fig.canvas.draw()
		time.sleep(step)

#compact a 3d map by a factor k
def low_resolution(map, k=2):
	l,m,n = map.shape
	map_red = map.reshape(l/k,k,m/k,k,n/k,k)
	map_red = np.sum(map_red, axis=5)
	map_red = np.sum(map_red, axis=3)
	map_red = np.sum(map_red, axis=1)
	return map_red/k**3

#short vesion of density field calculation (not tested yet)	
def short_density_field(t_set,k_set):
	t = time.time()
	f_t = np.fft.fftn(t_set,axes=[1,2,3])
	f_t = f_t * k_set
	kappa = f_t[0] + f_t[1] + 2*f_t[2]
	kappa = np.fft.ifftn(kappa).real
	t = time.time() - t
	print t
	return kappa

#make sure the min value of a delta map is > -1	
def rescale(map, epsilon=1e-5):
	return map*(-1+epsilon)/map.min()

#other solution
def cut(map, epsilon=1e-5):
	valid = map>=(-1+epsilon)
	return valid*map - (1-epsilon)*(1-valid)
